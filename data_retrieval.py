import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path
import boto3

import pandas as pd

from cred import s3, buck


# Fragment Data:
def get_frag_lengths2(frag_data: Path) -> dict:
    '''
    Returns dict of dict connecting
        key: fragment id
        values:
            fragment sequence
            if fragment is an element
    '''
    frag_id_2_frag_len = {}
    with frag_data.open('r') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            is_element = False
            if row['block_quant_file_name'][:2].upper() == 'PA' or row['plate_name'][:2].upper() == 'PA':
                is_element = True

            try:
                conc = str(round(float(row['concentration']), 2))
                if conc == '0.0' or is_element: conc = ''
            except ValueError:
                conc = ''
            frag_id = row['fragment_name']
            plate_name = row['block_quant_file_name'].lower().partition('.')[0]
            if frag_id not in frag_id_2_frag_len:
                frag_id_2_frag_len[frag_id] = {
                    'concentration': [conc],
                    'plate_name': [plate_name],
                    'fragment_position': [row['fragment_position']],
                    'fragment_sequence': [row['fragment_sequence']],
                    'is_element': [is_element]
                }
            else:
                frag_id_2_frag_len[frag_id]['concentration'].append(conc)
                frag_id_2_frag_len[frag_id]['plate_name'].append(plate_name)
                frag_id_2_frag_len[frag_id]['fragment_position'].append(row['fragment_position'])
                frag_id_2_frag_len[frag_id]['fragment_sequence'].append(row['fragment_sequence'])
                frag_id_2_frag_len[frag_id]['is_element'].append(is_element)
    return frag_id_2_frag_len


# Sequence ID + Batch to failed fragment count
def get_number_failed_frags(frag_data: pd.DataFrame, frag_id_2_frag_len) -> dict:
    '''
    Returns dict of sequence id:
        number of passing eBlocks
        Total fragments
        sequence tier
        fragment sequences
        element count
    '''
    seq_2_skip = set() #doesn't have proper frag data or concentrations are different between alphaprod and st1
    seq_id_2_frag_pass = {}
    for i, row in frag_data.iterrows():
        if row['frag id'][:2].isalpha():
            # skips controls that aren't recorded in alphaprod's block quant
            continue
        seq_id = str(row['Elegen Gene ID']) + '_' + str(row['Destination Plate ID'].partition('_')[0])#should be '[gene id]_[Batch_num without date]'
        passing = False
        try:
            if 'Conc' not in row:
                passing = row['eBlock qc flag'] == 'Pass'
            else:
                conc = float(row['Conc'])
                if conc == '0.0' or row['eBlock Plate ID'][:2].upper() == 'PA':
                    passing = row['eBlock qc flag'] == 'Pass'
                else:
                    passing = (conc > 50 and row['eBlock qc flag'] == 'Pass')
        except ValueError:
            passing = row['eBlock qc flag'] == 'Pass'

        # st1 named its fragments with [gene_id].1_frag_[1-9]
        try:
            conc = str(round(float(row['Conc']), 2))
            if conc == '0.0' or row['eBlock Plate ID'][:2] == 'pa': conc = ''
        except:
            conc = ''
        frag_id = row['frag id']
        dot_idx = row['frag id'].find('.')
        if dot_idx != -1:
            frag_id = frag_id[:dot_idx] + frag_id[dot_idx+2:]
            # couple of st1's that have fragments >=10 that are marked as frag_1. Renaming to match alphaprod
            tmp = 0
            while frag_id not in frag_id_2_frag_len:
                frag_id = row['frag id'][:dot_idx] + row['frag id'][dot_idx+2:] + str(tmp) + conc
                tmp += 1
                if tmp > 10:
                    break
        if frag_id not in frag_id_2_frag_len:
                print(frag_id)
                seq_2_skip.add(seq_id)
                continue
        # Trying to connect sequences to their fragments any way possible:
        frag_sequence = ''
        element = 0
        plate_name = row['eBlock Plate ID'].lower()
        if row['eBlock Plate ID'][-7:] == 'control':
            ctrl_char = 'a'
            plate_name = plate_name[:-8] + '_control_' + ctrl_char
            while plate_name not in frag_id_2_frag_len[frag_id]['plate_name']:
                ctrl_char = chr(ord(ctrl_char) + 1)
                plate_name = plate_name[:-1] + ctrl_char
                if ctrl_char == 'h':
                    break
        if conc and conc in frag_id_2_frag_len[frag_id]['concentration']:
            idx = frag_id_2_frag_len[frag_id]['concentration'].index(conc)
            frag_sequence = frag_id_2_frag_len[frag_id]['fragment_sequence'][idx]
            element = int(frag_id_2_frag_len[frag_id]['is_element'][idx])
        elif row['eBlock Plate Well Position'] in frag_id_2_frag_len[frag_id]['fragment_position']:
            idx = frag_id_2_frag_len[frag_id]['fragment_position'].index(row['eBlock Plate Well Position'])
            frag_sequence = frag_id_2_frag_len[frag_id]['fragment_sequence'][idx]
            element = int(frag_id_2_frag_len[frag_id]['is_element'][idx])
        elif plate_name in frag_id_2_frag_len[frag_id]['plate_name']:
            idx = frag_id_2_frag_len[frag_id]['plate_name'].index(plate_name)
            frag_sequence = frag_id_2_frag_len[frag_id]['fragment_sequence'][idx]
            element = int(frag_id_2_frag_len[frag_id]['is_element'][idx])
        else:
            seq_2_skip.add(seq_id)
            continue

        if seq_id in seq_id_2_frag_pass:
            seq_id_2_frag_pass[seq_id]['passing_frags'].append(passing)
            seq_id_2_frag_pass[seq_id]['total_frags'] += 1
            seq_id_2_frag_pass[seq_id]['frag_sequences'].append(frag_sequence)
            seq_id_2_frag_pass[seq_id]['element_count'].append(element)
            seq_id_2_frag_pass[seq_id]['concentrations'].append(conc)
        else:
            if frag_id not in frag_id_2_frag_len or seq_id in seq_2_skip:
                print(frag_id)
                continue
            seq_id_2_frag_pass[seq_id] = {
                'sequence': row['Gene Seq'].upper(),
                'passing_frags': [passing],
                'total_frags': 1,
                'tier': ['Synthesis Tier'],
                'frag_sequences': [frag_sequence],
                'element_count': [element],
                'concentrations': [conc]
            }
    return seq_id_2_frag_pass


# Using st_final files found from sql (bfx) system
# Sequence-ID to Confirmation Data
# connecting 7454 / 10010 sequences with existing fragment data to its confirmation data
def get_conf_data(
        df: pd.DataFrame,
        batch: str,
        date: str
    ):
    '''
    Returns dict of sequence id:
        wells passing (i.e., passing dial-out),
        total wells,
        date
    '''
    seq_id_2_conf_pass = {}
    for i, row in df.iterrows():
        seq_id = str(row['Elegen gene ID']) + '_' + batch
        try:
            passing = row['sequence_err'].lower() == 'pass' and row['ngs_overall'].lower() == 'pass' and row['sp1_pass'].lower() == 'pass' and row['sp2_pass'].lower() == 'pass'
        except:
            passing = row['ngs_overall'].lower() == 'pass' and row['sp1_pass'].lower() == 'pass' and row['sp2_pass'].lower() == 'pass'
        if seq_id in seq_id_2_conf_pass:
            seq_id_2_conf_pass[seq_id]['wells_passing'] += int(passing)
            seq_id_2_conf_pass[seq_id]['total_wells'] += 1
        else:
            seq_id_2_conf_pass[seq_id] = {
                'wells_passing': int(passing),
                'total_wells': 1,
                'date': date,
                'barcode_count': 0,
                'deletions': []#just saving top 5
            }
    return seq_id_2_conf_pass
