import logging
import numpy as np
import secrets
from collections import defaultdict
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import re
from pathlib import Path
import gzip
import shutil
import math
from math import log2
import os
import ast
import boto3

import pandas as pd
import sys
import json
from difflib import SequenceMatcher

import speedy

from cred import s3, buck
from dp_path_finder import find_graph, process_output_path

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


def longest_common_substring(s1, s2):
    m = len(s1)
    n = len(s2)
    # Initialize a 2D array to store lengths of longest common suffixes
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    max_length = 0
    end = 0
    dist5 = 0
    dist3 = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # ensuring it's off the exact alignment
            if i == j:
                continue

            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] > max_length:
                    max_length = dp[i][j]
                    end = i - 1
                    dist5 = i
                    dist3 = j - max_length
            else:
                dp[i][j] = 0

    # Extract the longest common substring
    return s1[end - max_length + 1: end + 1], dist5, dist3


def ATGC_energy(s1, mismatch_penalty=0.5) -> int:
    score = 0
    for c in s1:
        if c == 'A' or c == 'T':
            score += 2
        elif c == 'G' or c == 'C':
            score += 3
        else:
            score -= mismatch_penalty
    return score


def lcs(S, T):
    m = len(S)
    n = len(T)
    counter = [[0] * (n + 1) for _ in range(m + 1)]
    longest = 0
    lcs_set = set()
    
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    # lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])
            else:
                counter[i+1][j+1] = 0
    return lcs_set

def top_5_longest_common_substrings(S, T):
    lcs_set = lcs(S, T)
    tmp_list = []
    for cs in lcs_set:
        tmp_list.append((cs, ATGC_energy(cs)))
    lcs_list = sorted(tmp_list, key=lambda x: x[1], reverse=True)
    return lcs_list[:5]


def reprocess_well(seq_dels: dict, gene_batch_id: str, frag_conc):
    deletion_groups = []#[[[list of strts for group], [list of ends for group], group_percent_of_reads]]
    errorneous_deletions = set()
    # iterate over all exact deletions
    for coordinates,data_fields in seq_dels.items():
        # log2 of percent reads with exact deletion
        data_fields['Log2_exact_deletion'] = max(log2(data_fields['percentage_of_reads']), -15)

        # nullify old values for furthest right/left; this will change after combining wells
        data_fields['furthest_right'] = False
        data_fields['furthest_left'] = False

        # save bases at beginning/end of deletion because they'll be needed but are removed in final csv for saving space
        # bases_5prime = data_fields['bases_at_beginning_of_deletion']
        # bases_3prime = data_fields['bases_at_end_of_deletion']
        fiveprime_end = data_fields['bases_at_beginning_of_deletion'][:50] + data_fields['bases_at_beginning_of_deletion'][51:]
        threeprime_end = data_fields['bases_at_end_of_deletion'][:50] + data_fields['bases_at_end_of_deletion'][51:]

        # increase boundary to put repeat on right side of deletion; doing this for re-defining longest repeat on left-hand side of deletion
        lR = 50
        while lR < min(len(fiveprime_end), len(threeprime_end)) and fiveprime_end[lR] == threeprime_end[lR]:
            lR += 1
        rR = lR
        rep, dist5, dist3 = longest_common_substring(fiveprime_end[:lR], threeprime_end[:rR])
        data_fields['longest_repeat_before_del_end_del'] = rep
        data_fields['dist_before_del-longest_repeat_before_del_end_del'] = dist5 - lR
        data_fields['dist_end_del-longest_repeat_before_del_end_del'] = lR - dist3

        lL = 49
        while lR >= 0 and lR < min(len(fiveprime_end), len(threeprime_end)) and fiveprime_end[lL] == threeprime_end[lL]:
            lL -= 1
        rL = lL
        rep, dist5, dist3 = longest_common_substring(fiveprime_end[lL:], threeprime_end[rL:])
        data_fields['longest_repeat_beginning_del_after_del'] = rep
        data_fields['dist_beg_del-longest_repeat_beginning_del_after_del'] = dist5
        data_fields['dist_after_del-longest_repeat_beginning_del_after_del'] = -dist3

        # New fields for repeat energy scores
        if data_fields['longest_repeat_beginning_del_after_del'] == 5:
            data_fields['longest_repeat_beginning_del_after_del'] = ''
            # ^for a single line error. Figure out what happened there!
        data_fields['energy_longest_repeat_beginning_del_after_del'] = ATGC_energy(data_fields['longest_repeat_beginning_del_after_del'])
        data_fields['energy_longest_repeat_before_del_end_del'] = ATGC_energy(data_fields['longest_repeat_before_del_end_del'])

        # re-create deletion groups
        added_group = False
        for i, (group_strt, group_end, _) in enumerate(deletion_groups):
            avg_strt = sum(group_strt) / len(group_strt)
            avg_end = sum(group_end) / len(group_end)
            # new method for grouping deletions; old method was +/- 100 instead of 50
            if abs(avg_strt - data_fields['del_start']) <= 50 and abs(avg_end - data_fields['del_end']) <= 50:
                deletion_groups[i][0].append(data_fields['del_start'])
                deletion_groups[i][1].append(data_fields['del_end'])
                deletion_groups[i][2] += data_fields['percentage_of_reads']
                added_group = True
                break
        # didn't find deletion group? create a new one...
        if not added_group:
            deletion_groups.append([
                [data_fields['del_start']],
                [data_fields['del_end']],
                data_fields['percentage_of_reads']
            ])

    # calculate log2(percentage of total reads in deletion group)
    for i in range(len(deletion_groups)):
        # append current deletion group percentage for totalling overlapping deletion groups' percentages
        deletion_groups[i].append(deletion_groups[i][2])
        midpt = (sum(deletion_groups[i][0]) / len(deletion_groups[i][0]) + sum(deletion_groups[i][1]) / len(deletion_groups[i][1])) / 2.
        # find other deletion groups overlapping this one's midpoint
        for j in range(len(deletion_groups)):
            if i == j: continue
            avg_strt = sum(deletion_groups[j][0]) / len(deletion_groups[j][0])
            avg_end = sum(deletion_groups[j][1]) / len(deletion_groups[j][1])
            # If new group overlapping old group: sum percentages of deletions spanning this region
            if avg_strt < midpt < avg_end:
                deletion_groups[i][3] += deletion_groups[j][2]
    # have to restart the for-loop otherwise the percentages will get changed to logs and screw up the percent spanning deletions calculations for other groups
    for i in range(len(deletion_groups)):
        # log2 of total fraction of all exact-deletions in this deletion group
        deletion_groups[i][2] = max(-15, log2(deletion_groups[i][2]))
        # log2 of total fraction of reads that are NOT in spanning deletion group
        if deletion_groups[i][3] >= 1 and deletion_groups[i][3] - 1.09 < 0:
            deletion_groups[i][3] = -15
            continue
        deletion_groups[i][3] = max(-15, log2(1 - deletion_groups[i][3]))

    # Going back and saving this information into deletion dataset
    for coordinates,data_fields in seq_dels.items():
        added_group = False
        for i, (group_strt, group_end, group_log2_percent_of_reads, log2_total_reads_not_in_spanning_del_group) in enumerate(deletion_groups):
            if data_fields['del_start'] in group_strt and data_fields['del_end'] in group_end:
                data_fields['Log2_deletion_group_percent_of_reads'] = group_log2_percent_of_reads
                data_fields['Log2_total_reads_not_in_spanning_del_group'] = log2_total_reads_not_in_spanning_del_group

                # recording deletion group coordinates
                data_fields['group_start'] = sum(group_strt) / len(group_strt)
                data_fields['group_end'] = sum(group_end) / len(group_end)

                # Adding fragment concentrations as requested
                data_fields['fragment_concentrations'] = frag_conc

                # number of wells deletion is in
                data_fields['appears_in_x_wells'] = 1

                added_group = True
                break
        if not added_group:
            # exception for json messing up formatting once the line got too long
            if coordinates in errorneous_deletions:
                continue
            # should never raise valueerror unless there's a bug
            raise ValueError('Deletion group not found')

    return seq_dels


# Secondary Structure Code from helminth/element_secondary/sequence_attributes

def ATGC_secondary_structures(s1: str, s2: str, loop_len = 0, mismatch_penalty= 0.5) -> float:
    i1,i2 = 0,0
    score = 0
    while i1 < len(s1) and i2 < len(s2):
        if s1[i1] == s2[i2]:
            if s1[i1] == 'A' or s1[i1] == 'T':
                score += 2
            else:
                score += 3
        else:
            score -= mismatch_penalty
        i1 += 1
        i2 += 1

    if loop_len > 3:
        score -= mismatch_penalty*(loop_len-3)
    return score

def rev_comp(seq: str):
    revs = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    rev_seq = ''.join(x for x in seq[::-1])
    ret = ''
    for c in rev_seq:
        ret += revs[c]
    return ret


def find_repeats(sequence: str) -> list:
    repeat_dict = speedy.find_repeats(sequence)
    repeat_list = []
    for length in repeat_dict:
        if len(repeat_dict[length]) > 0 and (length in [2, 3] or length >= 8):
            for r in repeat_dict[length]:
                repeat_list.append(
                    {
                        "id": secrets.token_hex(4),
                        "positions": sorted(r),
                        "length": length,
                        "sequence": sequence[r[0] : r[0] + length],
                        "ATGC_score": ATGC_secondary_structures(sequence[r[0] : r[0] + length], sequence[r[0] : r[0] + length])
                    }
                )
    return sorted(repeat_list, key=lambda x: x['ATGC_score'], reverse=True)


def find_near_repeats(sequence):
    repeat_init = speedy.find_near_repeats(sequence)
    near_repeat_list = []
    for r in repeat_init:
        near_repeat_list.append(
            {
                "id": secrets.token_hex(4),
                "positions": sorted([r["position1"], r["position2"]]),
                "length": r["length"],
                "match_count": r["match_count"],
                "mismatch_count": r["mismatch_count"],
                "sequence": sequence[r["position1"] : r["position1"] + r["length"]],
                "ATGC_score": ATGC_secondary_structures(sequence[r["position1"] : r["position1"] + r["length"]],
                                                        rev_comp(sequence[r["position2"] : r["position2"] + r["length"]]))
            }
        )
    return sorted(near_repeat_list, key=lambda x: x['ATGC_score'], reverse=True)


def find_hairpins(sequence):
    hairpin_dict = speedy.find_hairpins(sequence)
    hairpin_list = []
    for length in hairpin_dict:
        if len(hairpin_dict[length]) > 0:
            for h in hairpin_dict[length]:
                hairpin_list.append(
                    {
                        "id": secrets.token_hex(4),
                        "positions": sorted(h),
                        "length": length,
                        "loop_length": h[1] - (h[0] + length),
                        "sequence": sequence[h[0] : h[0] + length],
                        "ATGC_score": ATGC_secondary_structures(sequence[h[0] : h[0] + length],
                                                                sequence[h[0] : h[0] + length],
                                                                h[1] - (h[0] + length))
                    }
                )
    return sorted(hairpin_list, key=lambda x: x['ATGC_score'], reverse=True)


def find_near_hairpins(sequence):
    hairpin_init = speedy.find_near_hairpins(sequence)
    near_hairpin_list = []
    for r in hairpin_init:
        near_hairpin_list.append(
            {
                "id": secrets.token_hex(4),
                "positions": sorted([r["position1"], r["position2"]]),
                "length": r["length"],
                "loop_length": r["position2"] - r["position1"] - r["length"],
                "match_count": r["match_count"],
                "mismatch_count": r["mismatch_count"],
                "sequence": sequence[r["position1"] : r["position1"] + r["length"]],
                "ATGC_score": ATGC_secondary_structures(sequence[r["position1"] : r["position1"] + r["length"]],
                                                        rev_comp(sequence[r["position2"] : r["position2"] + r["length"]]),
                                                        r["position2"] - r["position1"] - r["length"])
            }
        )
    return sorted(near_hairpin_list, key=lambda x: x['ATGC_score'], reverse=True)


def fix_energy_score(matching_over_deletion: str, mismatch_penalty=0.5):
    right_energy, right_offset, current_score = 0, 0, 0
    not_found_matching = True
    for c in matching_over_deletion[21:]:
        if c == 'A' or c == 'T':
            current_score += 2
            not_found_matching = False
        elif c == 'G' or c == 'C':
            current_score += 3
            not_found_matching = False
        else:
            current_score -= mismatch_penalty
            if not_found_matching:
                right_offset += 1
        right_energy = max(right_energy, current_score)

    left_offset, left_energy, current_score = 0, 0, 0
    not_found_matching = True
    if right_offset == 0:
        j = 21
        while j < len(matching_over_deletion) and matching_over_deletion[j] != '_':
            left_offset -= 1
            if matching_over_deletion[j] == 'A' or matching_over_deletion[j] == 'T':
                current_score += 2
            elif matching_over_deletion[j] == 'G' or matching_over_deletion[j] == 'C':
                current_score += 3
            j += 1
        not_found_matching = False
    for c in matching_over_deletion[19:0:-1]:
        if c == 'A' or c == 'T':
            current_score += 2
            not_found_matching = False
        elif c == 'G' or c == 'C':
            current_score += 3
            not_found_matching = False
        else:
            current_score -= mismatch_penalty
            if not_found_matching:
                left_offset += 1
        left_energy = max(left_energy, current_score)
    return right_energy, right_offset, left_energy, left_offset



def find_only_hairpins(del_segment: str) -> dict:
    # hairpins = find_hairpins(del_segment)
    # near_hairpins = find_near_hairpins(del_segment)
    try:
        hairpins = find_hairpins(del_segment)
        near_hairpins = find_near_hairpins(del_segment)
    except:
        # print(del_segment, len(del_segment))
        # print("HAIRPIN ERROR. Segment Length: ", len(del_segment))
        return {}

    # repeats, hairpins, near_repeats, near_hairpins = [], [], [], []
    # ToDo: calculate energies and re-order them
    return {
        'hairpins': hairpins,
        'near_hairpins': near_hairpins
    }


def find_only_repeats(del_segment: str) -> list:
    # repeats = find_repeats(del_segment)
    try:
        repeats = find_repeats(del_segment)
    except:
        # print(del_segment, len(del_segment))
        return []
    return repeats


def find_only_near_repeats(del_segment: str) -> list:
    # repeats = find_near_repeats(del_segment)
    try:
        repeats = find_near_repeats(del_segment)
    except:
        # print(del_segment, len(del_segment))
        return []
    return repeats



def change_hairpin_notation(hr: dict, origin=0, invert=1):
    summary = []
    for i, h in enumerate(hr):
        summary.append({
            'start': (h['positions'][0] - origin)*invert,
            'end': (h['positions'][1]+h['length'] - origin)*invert,
            'stem': h['length'],
            'loop': h['loop_length'],
            'energy': h['ATGC_score']
        })
        if 'mismatch_count' in h:
            summary[i]['mismatch_count'] = h['mismatch_count']
    return summary


def longest_substring(s1, s2) -> str:
    # Create a SequenceMatcher object with the two strings
    seq_match = SequenceMatcher(None, s1, s2)

    # Find the longest matching substring
    match = seq_match.find_longest_match(0, len(s1), 0, len(s2))

    # If a match is found, return the substring
    if match.size != 0:
        return s1[match.a: match.a + match.size]
    else:
        return ''


def process_all_wells(cw: dict, well_cnt: int, gene_batch_id: str, seq: str):

    '''
    Input:
        cw (combined-wells): dict of deletion data for multiple wells of the same gene+batch
        well_cnt: number of wells attributed to the sequence
    Returns:
        cw: updated combined-wells
    '''
    # COMBINE EXACT DELETIONS
    # Iterating each well
    consolidated_deletions = {}
    for well in cw:
        # Iterating each exact deletion
        for coordinates,data_fields in well.items():
            # CHECK if deletions coordinates are identical AND del-groups are similar
            if coordinates in consolidated_deletions and (
                abs(consolidated_deletions[coordinates]['group_start'] - data_fields['group_start']) <= 100 and
                abs(consolidated_deletions[coordinates]['group_end'] - data_fields['group_end']) <= 100 and
                consolidated_deletions[coordinates]['group_start'] < data_fields['group_end'] and
                consolidated_deletions[coordinates]['group_end'] < data_fields['group_start']
            ):
                
                consolidated_deletions[coordinates]['appears_in_x_wells'] += 1
                data_fields['appears_in_x_wells'] += 1 #===================================================================================================================== COULD VERY MUCH BE WRONG =====================================================================================================================
                consolidated_deletions[coordinates]['Log2_exact_deletion'] += data_fields['Log2_exact_deletion']
                consolidated_deletions[coordinates]['perfect_sequence_across_deletion'] |= (data_fields['perfect_sequence_across_deletion'].lower() == 'true')
                frac = data_fields['percentage_of_reads']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(frac/(1-frac))'] += 10 + log2(max(2**-15, frac / (1 - frac))) if frac != 1 else 2**63 - 1
                group_frac = data_fields['group_percent']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(group_frac/(1-group_frac))'] += 10 + log2(max(2 **-15, group_frac / (1 - group_frac))) if group_frac != 1 else 2**63 - 1

            else:
                consolidated_deletions[coordinates] = data_fields
                consolidated_deletions[coordinates]['perfect_sequence_across_deletion'] = (data_fields['perfect_sequence_across_deletion'].lower() == 'true')
                frac = data_fields['percentage_of_reads']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(frac/(1-frac))'] = 10 + log2(max(2**-15, frac / (1 - frac))) if frac != 1 else 2**63 - 1
                group_frac = data_fields['group_percent']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(group_frac/(1-group_frac))'] = 10 + log2(max(2**-15, group_frac / (1 - group_frac))) if group_frac != 1 else 2**63 - 1
            if group_frac < 0.05:
                print(' ======================================= CHECK ERROR IN GROUP FRAC RECORDING ==================================================')
    
    # Establish average group deletions across combined wells' deletions
    group_del_groups = []#[[group_strt, group_end, log2(frequency), well_ct_with_group ]]
    for coords,cd_data in consolidated_deletions.items():
        cd_data['Log2_exact_deletion'] /= cd_data['appears_in_x_wells']
        cd_data['Avg_wells_10+log2(frac/(1-frac))'] /= cd_data['appears_in_x_wells']
        cd_data['Avg_wells_10+log2(group_frac/(1-group_frac))'] /= cd_data['appears_in_x_wells']

        added_group = False
        for i, (group_strts, group_ends, log2_group_percent, exact_del_data, ct) in enumerate(group_del_groups):
            avg_group_strt = sum(group_strts) / len(group_strts)
            avg_group_end = sum(group_ends) / len(group_ends)
            if abs(avg_group_strt - cd_data['group_start']) <= 100 and abs(avg_group_end - cd_data['group_end']) <= 100:
                group_del_groups[i][0].add(cd_data['group_start'])
                group_del_groups[i][1].add(cd_data['group_end'])
                group_del_groups[i][2] += cd_data['Log2_deletion_group_percent_of_reads']
                group_del_groups[i][3].append(cd_data)#just store the whole group's individual deletions
                group_del_groups[i][4] += 1
                added_group = True
                break
        if not added_group:
            group_del_groups.append(
                [{cd_data['group_start']},# make into set to not double count group boundaries
                    {cd_data['group_end']},
                    cd_data['Log2_deletion_group_percent_of_reads'],
                    [cd_data],#just store the whole group's individual deletions
                    1]#count
            )
    for i, (group_strts, group_ends, log2_group_percent, exact_del_data, ct) in enumerate(group_del_groups):
        log2_group_percent /= ct

    # # Finding secondary structures in each group
    # group_sec_struct = []
    # avg_group_sec_struct = []
    # for i, (strts, ends, log2_group_percent, exact_del_data) in enumerate(group_del_groups):
    #     group_del_groups[i][2] /= well_cnt

    #     # Get Fragment Bounds:
    #     if exact_del_data[0]['distance_to_fragment_overhang_5prime'] > 0:
    #         frag_strt = exact_del_data[0]['distance_to_fragment_overhang_5prime']
    #     else:
    #         frag_strt = exact_del_data[0]['del_start']
    #     if exact_del_data[0]['distance_to_fragment_overhang_3prime'] > 0:
    #         frag_end = exact_del_data[0]['distance_to_fragment_overhang_3prime']
    #     else:
    #         frag_end = exact_del_data[0]['del_end']

        # # Secondary Structures between narrowest deletion segment
        # rightmost_strt = int(max(strts))
        # leftmost_end = int(min(ends))
        # group_del_groups[i].append(rightmost_strt)
        # group_del_groups[i].append(leftmost_end)


        # group_sec_struct.append(
        #     find_secondary_structures(seq[rightmost_strt-25:leftmost_end+25], seq[frag_strt:frag_end])
        # )

        # # Secondary Structures between avg deletion section
        # avg_strt = int(sum(strts) / len(strts))
        # avg_end = int(sum(ends) / len(ends))
        # avg_group_sec_struct.append(
        #     find_secondary_structures(seq[avg_strt-25:avg_end+25], seq[frag_strt:frag_end])
        # )

        # mid_pt = (avg_strt + avg_end) / 2.
        # current_frequency = log2_group_percent
        # for j, _ in enumerate(group_del_groups):
        #     if i == j:
        #         continue
        #     avg_group_strt = sum(group_del_groups[j][0]) / len(group_del_groups[j][0])
        #     avg_group_end = sum(group_del_groups[j][1]) / len(group_del_groups[j][1])
        #     if avg_group_strt <= mid_pt and avg_group_end >= mid_pt:
        #         current_frequency += group_del_groups[j][2]
        # group_del_groups[i].append(1 - current_frequency)

    # At This Point
    # group_del_groups = [
    #       [0. [group_starts],
    #       1.[group_ends],
    #       2. log2_group_frequency,
    #       3. [exact-del data]
    #       4. rightmost_strt,
    #       5. leftmost_strt,
    #       6. frequency(not in spanning-deletion)],
    #       ...
    #       [...]
    # ]

    # Now re-calculate top10 deletions + rightmost/leftmost exact-deletion for each group
    gene_id, batch_num = gene_batch_id.rsplit('_', 1)
    # dumpable_data = {gene_id: {}}
    # for i, (strts, ends, log2_group_percent, exact_del_data, rightmost_strt, leftmost_end, frequency) in enumerate(group_del_groups):
    for i, (strts, ends, log2_group_percent, exact_del_data, _) in enumerate(group_del_groups):
        exact_del_data = sorted(exact_del_data, key=lambda x: x['Log2_exact_deletion'], reverse=True)

        # TOP 10 DELETIONS
        for exact_del in exact_del_data[:10]:
            dumpable_data = {gene_id: {}}
            # Exact-Del Secondary Structures
            fiveprime_end = exact_del['bases_at_beginning_of_deletion'][:50] + exact_del['bases_at_beginning_of_deletion'][51:]
            threeprime_end = exact_del['bases_at_end_of_deletion'][:50] + exact_del['bases_at_end_of_deletion'][51:]
            
            # lL, lR, rL, rR are all offsets from the deletion boundary marked by minimap
            #   all 4 values are absolute (positive) and represent the range of possible "real" deletion boundaries
            #   between lL to lR and rL to rR should be a repeat that makes the deletion boundary ambiguous
            lR, rR = 0, 0
            idx = 50
            while idx < len(threeprime_end) and fiveprime_end[idx] == threeprime_end[idx]:
                idx += 1
                lR += 1
                rR += 1

            lL, rL = 0, 0
            idx = 49
            while idx >= 0 and fiveprime_end[idx] == threeprime_end[idx]:
                idx -= 1
                lL += 1
                rL += 1
            
            exact_del['LL'] = -lL
            exact_del['LR'] = lR
            exact_del['RL'] = rL
            exact_del['RR'] = -rR
            ds = exact_del['del_start']
            de = exact_del['del_end']

            seq1_dp_left = seq[ds+lR-50:ds+lR]#inboard
            seq2_dp_left = seq[de+rR-50:de+rR]#inside deletion
            dpLeftPath, dpLeftEnergy, dpLeftSequence = find_graph(seq1_dp_left,seq2_dp_left)
            dpLeftPathDict = process_output_path(dpLeftPath)

            seq2_dp_right = seq[ds-lL:ds-lL+50]#inside deletion
            seq1_dp_right = seq[de-rL:de-rL+50]#inboard
            dpRightPath, dpRightEnergy, dpRightSequence = find_graph(seq1_dp_right,seq2_dp_right)
            dpRightPathDict = process_output_path(dpRightPath)
            
            ds = exact_del['del_start']
            de = exact_del['del_end']
            if exact_del['distance_to_fragment_overhang_5prime'] > 0:
                fs = ds - exact_del['distance_to_fragment_overhang_5prime']
            else:
                fs = ds + exact_del['distance_to_fragment_overhang_5prime']
            if exact_del['distance_to_fragment_overhang_3prime'] > 0:
                fe = de + exact_del['distance_to_fragment_overhang_3prime']
            else:
                fe = de - exact_del['distance_to_fragment_overhang_5prime']

            # +/-50bp around LR:
            lr50_hairpins = find_only_hairpins(seq[ds+lR-50:ds+lR+50])
            if len(lr50_hairpins):
                sec_struct_2_keep = {}
                for k,v in lr50_hairpins.items():
                    sec_struct_2_keep[k] = []
                    for strukt in v[:5]:
                        sec_struct_2_keep[k].append(strukt)
                exact_del['5_prime_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'], 50)
                exact_del['5_prime_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'], 50)
            else:
                exact_del['5_prime_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
                exact_del['5_prime_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]

            # +/-50bp around RL:
            rl50_hairpins = find_only_hairpins(seq[de-lR-50:de-lR+50])
            if len(rl50_hairpins):
                sec_struct_2_keep = {}
                for k,v in rl50_hairpins.items():
                    sec_struct_2_keep[k] = []
                    for strukt in v[:5]:
                        sec_struct_2_keep[k].append(strukt)
                exact_del['3_prime_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'], 50, -1)
                exact_del['3_prime_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'], 50, -1)
            else:
                exact_del['3_prime_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
                exact_del['3_prime_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]


            # # IHFL (Interfering with Hybridization from Left): 
            # #   LL-6<=HR<=LR+1 
            # IHFL_sec_struct = find_only_hairpins(seq[fs:ds+lR+30])
            # if len(IHFL_sec_struct):
            #     sec_struct_2_keep = {}
            #     for k,v in IHFL_sec_struct.items():
            #         sec_struct_2_keep[k] = []
            #         top2limit = 2
            #         for strukt in v:
            #             h_end = strukt['positions'][1] + strukt['length']
            #             h_len = strukt['length']
            #             if (
            #                 (ds-fs-lL <= h_end + 6 <= ds - fs + lR + h_len)
            #                 # and (strukt['ATGC_score'] > 0)
            #                 ):
            #                 # subtracting bases overlapping deletion from the ATGC energy score 
            #                 if h_end + 1 > ds - fs + lR:
            #                     bpc = h_end + 1 - (ds - fs + lR)#del-end - hairpin-end (should be positive if overlapping)
            #                     strukt['ATGC_score'] -= ATGC_energy(strukt['sequence'][:bpc])# for reverse complement, the first few bases of the left part of stem is what you'd be subtracting
            #                     # if strukt['ATGC_score'] <= 0: continue
            #                 sec_struct_2_keep[k].append(strukt)
            #                 top2limit -= 1
            #             if top2limit <= 0:
            #                 break
            #     exact_del['IHFL_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'], ds+lR-fs)
            #     exact_del['IHFL_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'], ds+lR-fs)
            # else:
            #     exact_del['IHFL_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
            #     exact_del['IHFL_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]

            # # PDFL (Polymerase Dropoff from Left):
            # #   LL-1<=HL<=LR+9
            # PDFL_sec_struct = find_only_hairpins(seq[ds-lL-1:fe])
            # if len(PDFL_sec_struct):
            #     sec_struct_2_keep = {}
            #     for k,v in PDFL_sec_struct.items():
            #         sec_struct_2_keep[k] = []
            #         top2limit = 2
            #         for strukt in v:
            #             if (
            #                 (1 <= strukt['positions'][0] <= lL+1+lR+10)
            #                 # and (strukt['ATGC_score'] > 0)
            #                 ):
            #                 sec_struct_2_keep[k].append(strukt)
            #                 top2limit -= 1
            #             if top2limit <= 0:
            #                 break
            #     exact_del['PDFL_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'], 1+lL+lR)
            #     exact_del['PDFL_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'], 1+lL+lR)
            # else:
            #     exact_del['PDFL_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
            #     exact_del['PDFL_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]

            # # IHFR (Interfering with Hybridization from Right): 
            # #   RL-1<=HL<=RR+9 
            # IHFR_sec_struct = find_only_hairpins(seq[de-rL-30:fe])
            # if len(IHFR_sec_struct):
            #     sec_struct_2_keep = {}
            #     for k,v in IHFR_sec_struct.items():
            #         sec_struct_2_keep[k] = []
            #         top2limit = 2
            #         for strukt in v:
            #             h_start = strukt['positions'][0]
            #             h_len = strukt['length']
            #             if (
            #                 (30-h_len <= h_start-6 <= rL+30+rR)
            #                 # and (strukt['ATGC_score'] > 0)
            #                 ):
            #                 # subtracting bases overlapping deletion from the ATGC energy score 
            #                 if h_start < 30 + 1:
            #                     bpc = 30 - h_start + 1
            #                     strukt['ATGC_score'] -= ATGC_energy(strukt['sequence'][:bpc])
            #                     # if strukt['ATGC_score'] <= 0: continue
            #                 sec_struct_2_keep[k].append(strukt)
            #                 top2limit -= 1
            #             if top2limit <= 0:
            #                 break
            #     exact_del['IHFR_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'], 30, invert=-1)
            #     exact_del['IHFR_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'], 30, invert=-1)
            # else:
            #     exact_del['IHFR_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
            #     exact_del['IHFR_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]

            # # PDFR (Polymerase Dropoff from Right):
            # #   RL-9<=HR<=RR+1
            # PDFR_sec_struct = find_only_hairpins(seq[fs:de+rR+2])
            # if len(PDFR_sec_struct):
            #     sec_struct_2_keep = {}
            #     for k,v in PDFR_sec_struct.items():
            #         sec_struct_2_keep[k] = []
            #         top2limit = 2
            #         for strukt in v:
            #             h_end = strukt['positions'][1] + strukt['length']
            #             if (
            #                 (de-fs-rL-10 <= h_end <= de-fs+rR)
            #                 # and (strukt['ATGC_score'] > 0)
            #                 ):
            #                 sec_struct_2_keep[k].append(strukt)
            #                 top2limit -= 1
            #             if top2limit <= 0:
            #                 break
            #     exact_del['PDFR_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'], de-fs-rL, invert=-1)
            #     exact_del['PDFR_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'], de-fs-rL, invert=-1)
            # else:
            #     exact_del['PDFR_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
            #     exact_del['PDFR_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]

            # OVERALL SECONDARY STRUCTURES
            overall_hairpins = find_only_hairpins(seq[ds-20:de+20])
            if len(overall_hairpins):
                sec_struct_2_keep = {}
                for k,v in overall_hairpins.items():
                    sec_struct_2_keep[k] = []
                    top5limit = 5
                    for strukt in v:
                        if (strukt['positions'][0] < 20+lR and 20+lR-strukt['positions'][0] < strukt['length']/2) or ((strukt['positions'][1] + strukt['length'] > 20+de-ds-rL) and ((strukt['positions'][1]+strukt['length']) - (20+de-ds-rL) < strukt['length']/2)):
                            # if (strukt['ATGC_score'] > 0):
                            strukt['positions'][0] -= (20+lR)
                            strukt['positions'][1] -= (20+de-ds-rL)
                            strukt['positions'][1] *= -1
                            sec_struct_2_keep[k].append(strukt)
                            top5limit -= 1
                        if top5limit <= 0:
                            break
                exact_del['Overall_hairpins'] = change_hairpin_notation(sec_struct_2_keep['hairpins'])
                exact_del['Overall_near_hairpins'] = change_hairpin_notation(sec_struct_2_keep['near_hairpins'])
            else:
                exact_del['Overall_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''}]
                exact_del['Overall_near_hairpins'] = [{'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''}]
            # overall repeats
            overall_reps = find_only_repeats(seq[ds-lL-35:de+rR+35])
            overall_near_reps = find_only_near_repeats(seq[ds-lL-35:de+rR+35])
            if len(overall_reps):
                sec_struct_2_keep = []
                top5limit = 5
                for strukt in overall_reps:
                    # if (
                    #     (strukt['positions'][0] >= lR-35)
                    #     or (strukt['positions'][1]+strukt['length'] <= de-ds+65-rL)
                    #     # and (strukt['ATGC_score'] > 0)
                    #     ):
                    sec_struct_2_keep.append({
                        'start': strukt['positions'][0] - 35 - lL,
                        'end': 35 + de - ds + rR - strukt['positions'][1]+strukt['length'],
                        'length': strukt['length'],
                        'energy': strukt['ATGC_score']
                    })
                    top5limit -= 1
                    if top5limit <= 0:
                        break
                exact_del['Overall_repeats'] = sec_struct_2_keep
            else:
                exact_del['Overall_repeats'] = [{'start': '', 'end': '', 'length': '', 'energy': ''}]
            if len(overall_near_reps):
                sec_struct_2_keep = []
                top5limit = 5
                for strukt in overall_near_reps:
                    if (
                        (strukt['positions'][0] >= 5+lR) 
                        or (strukt['positions'][1]+strukt['length'] <= de-ds+65-rL)
                        # and (strukt['ATGC_score'] > 0)
                        ):
                        sec_struct_2_keep.append({
                            'start': strukt['positions'][0] - 35 - lR,
                            'end': 35 + de - ds - rL - strukt['positions'][1]+strukt['length'],
                            'length': strukt['length'],
                            'mismatch_count': strukt['mismatch_count'],
                            'energy': strukt['ATGC_score']
                        })
                        top5limit -= 1
                    if top5limit <= 0:
                        break
                exact_del['Overall_near_repeats'] = sec_struct_2_keep
            else:
                exact_del['Overall_near_repeats'] = [{'start': '', 'end': '', 'length': '', 'energy': '', 'mismatch_count': ''}]

            # ================================================== REVERSE REPEATS ==================================================
            overall_reverse_reps = top_5_longest_common_substrings(seq[ds-30+lR:de+30-rL], seq[ds-30+lR:de+30-rL:-1])
            if len(overall_reverse_reps):
                sec_struct_2_keep = []
                for (revrep, atgc_score) in overall_reverse_reps:
                    idx1_1 = seq[ds-30+lR:de+30-rL].find(revrep) - 30 - lR
                    idx1_2 = seq[ds-30+lR:de+30-rL].find(revrep[::-1]) - 30 - lR
                    idx2_1 = seq[ds-30+lR:de+30-rL:-1].find(revrep) - 30 - rL
                    idx2_2 = seq[ds-30+lR:de+30-rL:-1].find(revrep[::-1]) - 30 - rL
                    sec_struct_2_keep.append({
                        'start': min(idx1_1,idx1_2),
                        'end': min(idx2_1,idx2_2),
                        'length': len(revrep),
                        'energy': atgc_score
                    })

                exact_del['Overall_reverse_repeats'] = sec_struct_2_keep
            else:
                exact_del['Overall_reverse_repeats'] = [{'start': '', 'end': '', 'length': '', 'energy': ''}]



            # Fixing right/left energy+offset
            right_energy, right_offset, left_energy, left_offset = fix_energy_score(exact_del['matching_over_deletion'])
            
            # fixing longest_repeat offset
            long_rep = longest_substring(fiveprime_end, threeprime_end)
            exact_del['longest_repeat'] = long_rep
            fiveprime_offset = (len(long_rep) + fiveprime_end.find(long_rep)) - (50 + lR)
            threeprime_offset = (49 - rL) - (threeprime_end.find(long_rep) + len(long_rep))

            after_del = threeprime_end[50-rL:]
            beginning_del = fiveprime_end[50+lR:]
            longest_repeat_beginning_del_after_del = longest_substring(beginning_del, after_del)
            repeat_dist_beg_del = beginning_del.find(longest_repeat_beginning_del_after_del) + len(longest_repeat_beginning_del_after_del)#can be zero-indexed
            repeat_dist_after_del = -(1 + after_del.find(longest_repeat_beginning_del_after_del) + len(longest_repeat_beginning_del_after_del))

            before_del = fiveprime_end[:50+lR]
            end_del = threeprime_end[:50-rL]
            longest_repeat_before_del_end_del = longest_substring(before_del, end_del)
            repeat_dist_before_del = 50 + lR - (len(longest_repeat_before_del_end_del) + before_del.find(longest_repeat_before_del_end_del))
            repeat_dist_end_del = 50 - rL - (1 + len(longest_repeat_before_del_end_del) + end_del.find(longest_repeat_before_del_end_del))#can be zero-indexed

            # if rightmost_strt == exact_del['del_start']:
            #     exact_del['furthest_right'] = True
            # if leftmost_end == exact_del['del_end']:
            #     exact_del['furthest_left'] = True



            tmp = {
                # Sequence meta-data:
                'batch_num': batch_num,
                'total_wells': well_cnt,
                'Log2 Total Reads in Well': exact_del['log2_total_reads_well'],

                # Exact-Deletion Data
                'del_start': exact_del['del_start'],
                'del_end': exact_del['del_end'],
                'del_length': exact_del['del_end'] - exact_del['del_start'],
                
                'Avg_wells_10+log2(frac/(1-frac))': exact_del['Avg_wells_10+log2(frac/(1-frac))'],
                'Avg_wells_10+log2(group_frac/(1-group_frac))': exact_del['Avg_wells_10+log2(group_frac/(1-group_frac))'],
                'Left_End_Score': exact_del['left_end_score'],
                'Right_End_Score': exact_del['right_end_score'],
                'LL': exact_del['LL'],
                'LR': exact_del['LR'],
                'RL': exact_del['RL'],
                'RR': exact_del['RR'],

                '10bp_before_deletion': exact_del['bases_at_beginning_of_deletion'][40:50],
                '10bp_into_deletion_start': exact_del['bases_at_beginning_of_deletion'][51:61],
                '10bp_before_deletion_end': exact_del['bases_at_end_of_deletion'][40:50],
                '10bp_after_deletion': exact_del['bases_at_end_of_deletion'][51:61],
                'del_in_x_wells': exact_del['appears_in_x_wells'],
                # 'similar_to_past_deletion': exact_del['similar_to_past_deletion'],
                'read_ct': exact_del['read_ct'],
                'percentage_of_reads': exact_del['percentage_of_reads'],
                'matching_over_deletion': exact_del['matching_over_deletion'],
                'left_Single_Match': exact_del['left_matching_score_mismatch_1'],
                'left_Tail': exact_del['left_matching_score_match_1'],
                'left-Match_1': exact_del['left_matching_score_mismatch_2'],
                'left-GC_content_1': exact_del['left_match_group_GC_score_content_1'],
                'left-Mismatch_1': exact_del['left_matching_score_match_2'],
                'left-Match_2': exact_del['left_matching_score_mismatch_3'],
                'left-GC_content_2': exact_del['left_match_group_GC_score_content_2'],
                'left-Mismatch_2': exact_del['left_matching_score_match_3'],
                'left-Match_3': exact_del['left_matching_score_mismatch_4'],
                'left-GC_content_3': exact_del['left_match_group_GC_score_content_3'],
                'left-Mismatch_3': exact_del['left_matching_score_match_4'],
                'left-Match_4': exact_del['left_matching_score_mismatch_5'],
                'left-GC_content_4': exact_del['left_match_group_GC_score_content_4'],
                'left-Mismatch_4': exact_del['left_matching_score_match_5'],
                # 'left-Match_5': exact_del['left_matching_score_mismatch_2'],
                # 'left-GC_content_5': exact_del['left_match_group_GC_score_content_5'],
                # 'left-Mismatch_5': exact_del['left_matching_score_mismatch_2'],

                'right_Single_Match': exact_del['right_matching_score_mismatch_1'],
                'right_Tail': exact_del['right_matching_score_match_1'],
                'right-Match_1': exact_del['right_matching_score_mismatch_2'],
                'right-GC_content_1': exact_del['right_match_group_GC_score_content_1'],
                'right-Mismatch_1': exact_del['right_matching_score_match_2'],
                'right-Match_2': exact_del['right_matching_score_mismatch_3'],
                'right-GC_content_2': exact_del['right_match_group_GC_score_content_2'],
                'right-Mismatch_2': exact_del['right_matching_score_match_3'],
                'right-Match_3': exact_del['right_matching_score_mismatch_4'],
                'right-GC_content_3': exact_del['right_match_group_GC_score_content_3'],
                'right-Mismatch_3': exact_del['right_matching_score_match_4'],
                'right-Match_4': exact_del['right_matching_score_mismatch_5'],
                'right-GC_content_4': exact_del['right_match_group_GC_score_content_4'],
                'right-Mismatch_4': exact_del['right_matching_score_match_5'],

                'right_energy': right_energy,
                'right_offset': right_offset,
                'left_energy': left_energy,
                'left_offset': left_offset,
                'gc_content_frag_start_LR': exact_del['gc_content_frag_start_LR'],

                'longest_repeat': long_rep,
                'longest_repeat_dist_5prime': fiveprime_offset,
                'longest_repeat_dist_3prime': threeprime_offset,
                'longest_repeat_beginning_del_after_del': longest_repeat_beginning_del_after_del,
                'dist_beg_del-longest_repeat_beginning_del_after_del': repeat_dist_beg_del,
                'dist_after_del-longest_repeat_beginning_del_after_del': repeat_dist_after_del,
                'longest_repeat_before_del_end_del': longest_repeat_before_del_end_del,
                'dist_before_del-longest_repeat_before_del_end_del': repeat_dist_before_del,
                'dist_end_del-longest_repeat_before_del_end_del': repeat_dist_end_del,
                'longest_repeat_score': exact_del['longest_repeat_score'],
                'longest_a-d_ligation_repeat': exact_del['longest_a-d_ligation_repeat'],
                'within_one_fragment': exact_del['within_one_fragment'],
                'distance_to_fragment_overhang_5prime': exact_del['distance_to_fragment_overhang_5prime'],
                'distance_to_fragment_overhang_3prime': exact_del['distance_to_fragment_overhang_3prime'],

                'perfect_sequence_across_deletion': str(exact_del['perfect_sequence_across_deletion']),
                'energy_longest_repeat_beginning_del_after_del': exact_del['energy_longest_repeat_beginning_del_after_del'],
                'energy_longest_repeat_before_del_end_del': exact_del['energy_longest_repeat_before_del_end_del'],
                'Log2_exact_deletion_frequency': exact_del['Log2_exact_deletion'],
                'Log2_group_deletion_frequency': log2_group_percent#exact_del['Log2_deletion_group_frequency'],
            }


            # Exact deletion secondary structures
            while len(exact_del['5_prime_hairpins']) < 5:
                exact_del['5_prime_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            for b,c in enumerate(exact_del['5_prime_hairpins']):
                tmp[f'5_prime_hairpins_start_{b}'] = c['start'] 
                tmp[f'5_prime_hairpins_end_{b}'] = c['end'] 
                tmp[f'5_prime_hairpins_stem_{b}'] = c['stem']
                tmp[f'5_prime_hairpins_loop_{b}'] = c['loop']
                tmp[f'5_prime_hairpins_energy_{b}'] = c['energy']
            while len(exact_del['5_prime_near_hairpins']) < 5:
                exact_del['5_prime_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            for b,c in enumerate(exact_del['5_prime_near_hairpins']):
                tmp[f'5_prime_near_hairpins_start_{b}'] = c['start']
                tmp[f'5_prime_near_hairpins_end_{b}']=  c['end']
                tmp[f'5_prime_near_hairpins_stem_{b}'] = c['stem']
                tmp[f'5_prime_near_hairpins_loop_{b}'] =  c['loop']
                tmp[f'5_prime_near_hairpins_energy_{b}'] = c['energy']
                tmp[f'5_prime_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']
            

            while len(exact_del['3_prime_hairpins']) < 5:
                exact_del['3_prime_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            for b,c in enumerate(exact_del['3_prime_hairpins']):
                tmp[f'3_prime_hairpins_start_{b}'] = c['start'] 
                tmp[f'3_prime_hairpins_end_{b}'] = c['end'] 
                tmp[f'3_prime_hairpins_stem_{b}'] = c['stem']
                tmp[f'3_prime_hairpins_loop_{b}'] = c['loop']
                tmp[f'3_prime_hairpins_energy_{b}'] = c['energy']
            while len(exact_del['3_prime_near_hairpins']) < 5:
                exact_del['3_prime_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            for b,c in enumerate(exact_del['3_prime_near_hairpins']):
                tmp[f'3_prime_near_hairpins_start_{b}'] = c['start']
                tmp[f'3_prime_near_hairpins_end_{b}']=  c['end']
                tmp[f'3_prime_near_hairpins_stem_{b}'] = c['stem']
                tmp[f'3_prime_near_hairpins_loop_{b}'] =  c['loop']
                tmp[f'3_prime_near_hairpins_energy_{b}'] = c['energy']
                tmp[f'3_prime_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']

            #     # Exact deletion secondary structures
            # while len(exact_del['IHFL_hairpins']) < 2:
            #     exact_del['IHFL_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            # for b,c in enumerate(exact_del['IHFL_hairpins']):
            #     tmp[f'IHFL_hairpins_start_{b}'] = c['start'] 
            #     tmp[f'IHFL_hairpins_end_{b}'] = c['end'] 
            #     tmp[f'IHFL_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'IHFL_hairpins_loop_{b}'] = c['loop']
            #     tmp[f'IHFL_hairpins_energy_{b}'] = c['energy']
            # while len(exact_del['IHFL_near_hairpins']) < 2:
            #     exact_del['IHFL_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            # for b,c in enumerate(exact_del['IHFL_near_hairpins']):
            #     tmp[f'IHFL_near_hairpins_start_{b}'] = c['start']
            #     tmp[f'IHFL_near_hairpins_end_{b}']=  c['end']
            #     tmp[f'IHFL_near_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'IHFL_near_hairpins_loop_{b}'] =  c['loop']
            #     tmp[f'IHFL_near_hairpins_energy_{b}'] = c['energy']
            #     tmp[f'IHFL_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']

            # while len(exact_del['PDFL_hairpins']) < 2:
            #     exact_del['PDFL_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            # for b,c in enumerate(exact_del['PDFL_hairpins']):
            #     tmp[f'PDFL_hairpins_start_{b}'] = c['start'] 
            #     tmp[f'PDFL_hairpins_end_{b}'] = c['end'] 
            #     tmp[f'PDFL_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'PDFL_hairpins_loop_{b}'] = c['loop']
            #     tmp[f'PDFL_hairpins_energy_{b}'] = c['energy']
            # while len(exact_del['PDFL_near_hairpins']) < 2:
            #     exact_del['PDFL_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            # for b,c in enumerate(exact_del['PDFL_near_hairpins']):
            #     tmp[f'PDFL_near_hairpins_start_{b}'] = c['start']
            #     tmp[f'PDFL_near_hairpins_end_{b}']=  c['end']
            #     tmp[f'PDFL_near_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'PDFL_near_hairpins_loop_{b}'] =  c['loop']
            #     tmp[f'PDFL_near_hairpins_energy_{b}'] = c['energy']
            #     tmp[f'PDFL_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']

            # while len(exact_del['IHFR_hairpins']) < 2:
            #     exact_del['IHFR_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            # for b,c in enumerate(exact_del['IHFR_hairpins']):
            #     tmp[f'IHFR_hairpins_start_{b}'] = c['start'] 
            #     tmp[f'IHFR_hairpins_end_{b}'] = c['end'] 
            #     tmp[f'IHFR_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'IHFR_hairpins_loop_{b}'] = c['loop']
            #     tmp[f'IHFR_hairpins_energy_{b}'] = c['energy']
            # while len(exact_del['IHFR_near_hairpins']) < 2:
            #     exact_del['IHFR_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            # for b,c in enumerate(exact_del['IHFR_near_hairpins']):
            #     tmp[f'IHFR_near_hairpins_start_{b}'] = c['start']
            #     tmp[f'IHFR_near_hairpins_end_{b}'] = c['end']
            #     tmp[f'IHFR_near_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'IHFR_near_hairpins_loop_{b}'] =  c['loop']
            #     tmp[f'IHFR_near_hairpins_energy_{b}'] = c['energy']
            #     tmp[f'IHFR_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']

            # while len(exact_del['PDFR_hairpins']) < 2:
            #     exact_del['PDFR_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            # for b,c in enumerate(exact_del['PDFR_hairpins']):
            #     tmp[f'PDFR_hairpins_start_{b}'] = c['start'] 
            #     tmp[f'PDFR_hairpins_end_{b}'] = c['end'] 
            #     tmp[f'PDFR_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'PDFR_hairpins_loop_{b}'] = c['loop']
            #     tmp[f'PDFR_hairpins_energy_{b}'] = c['energy']
            # while len(exact_del['PDFR_near_hairpins']) < 2:
            #     exact_del['PDFR_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            # for b,c in enumerate(exact_del['PDFR_near_hairpins']):
            #     tmp[f'PDFR_near_hairpins_start_{b}'] = c['start']
            #     tmp[f'PDFR_near_hairpins_end_{b}'] = c['end']
            #     tmp[f'PDFR_near_hairpins_stem_{b}'] = c['stem']
            #     tmp[f'PDFR_near_hairpins_loop_{b}'] =  c['loop']
            #     tmp[f'PDFR_near_hairpins_energy_{b}'] = c['energy']
            #     tmp[f'PDFR_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']

            while len(exact_del['Overall_hairpins']) < 5:
                exact_del['Overall_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': ''})
            for b,c in enumerate(exact_del['Overall_hairpins']):
                tmp[f'Overall_hairpins_start_{b}'] = c['start'] 
                tmp[f'Overall_hairpins_end_{b}'] = c['end'] 
                tmp[f'Overall_hairpins_stem_{b}'] = c['stem']
                tmp[f'Overall_hairpins_loop_{b}'] = c['loop']
                tmp[f'Overall_hairpins_energy_{b}'] = c['energy']
            while len(exact_del['Overall_near_hairpins']) < 5:
                exact_del['Overall_near_hairpins'].append({'start': '', 'end': '', 'stem': '', 'loop': '', 'energy': '', 'mismatch_count': ''})
            for b,c in enumerate(exact_del['Overall_near_hairpins']):
                tmp[f'Overall_near_hairpins_start_{b}'] = c['start']
                tmp[f'Overall_near_hairpins_end_{b}'] = c['end']
                tmp[f'Overall_near_hairpins_stem_{b}'] = c['stem']
                tmp[f'Overall_near_hairpins_loop_{b}'] =  c['loop']
                tmp[f'Overall_near_hairpins_energy_{b}'] = c['energy']
                tmp[f'Overall_near_hairpins_mismatch_count_{b}'] = c['mismatch_count']
            
            while len(exact_del['Overall_repeats']) < 5:
                exact_del['Overall_repeats'].append({'start': '', 'end': '', 'length': '', 'energy': ''})
            for b,c in enumerate(exact_del['Overall_repeats']):
                tmp[f'Overall_repeats_start_{b}'] = c['start'] 
                tmp[f'Overall_repeats_end_{b}'] = c['end'] 
                tmp[f'Overall_repeats_length_{b}'] = c['length']
                tmp[f'Overall_repeats_energy_{b}'] = c['energy']
            while len(exact_del['Overall_near_repeats']) < 5:
                exact_del['Overall_near_repeats'].append({'start': '', 'end': '', 'length': '', 'energy': '', 'mismatch_count': ''})
            for b,c in enumerate(exact_del['Overall_near_repeats']):
                tmp[f'Overall_near_repeats_start_{b}'] = c['start']
                tmp[f'Overall_near_repeats_end_{b}'] = c['end']
                tmp[f'Overall_near_repeats_length_{b}'] = c['length']
                tmp[f'Overall_near_repeats_energy_{b}'] = c['energy']
                tmp[f'Overall_near_repeats_mismatch_count_{b}'] = c['mismatch_count']
            
            while len(exact_del['Overall_reverse_repeats']) < 5:
                exact_del['Overall_reverse_repeats'].append({'start': '', 'end': '', 'length': '', 'energy': ''})
            for b,c in enumerate(exact_del['Overall_reverse_repeats']):
                tmp[f'Overall_reverse_repeats_start_{b}'] = c['start'] 
                tmp[f'Overall_reverse_repeats_end_{b}'] = c['end'] 
                tmp[f'Overall_reverse_repeats_length_{b}'] = c['length']
                tmp[f'Overall_reverse_repeats_energy_{b}'] = c['energy']

                # # Deletion Group Data
                # 'group_strt': exact_del['group_start'],
                # 'group_end': exact_del['group_end'],
                # 'Log2_deletion_group_frequency': exact_del['Log2_deletion_group_percent_of_reads'],
                # 'Log2_not_in_spanning_del_group_frequency': exact_del['Log2_total_reads_not_in_spanning_del_group'],
                # # 'furthest_right': str(exact_del['del_start'] == exact_del['furthest_right_coordinates']),
                # # 'furthest_left': str(exact_del['del_end'] == exact_del['furthest_left_coordinates']),
                # 'furthest_right': str(exact_del['del_start'] == rightmost_strt),
                # 'furthest_left': str(exact_del['del_end'] == leftmost_end),
                # 'group_repeats': str(group_sec_struct[i]['repeats'][:5]),
                # 'group_near_repeats': str(group_sec_struct[i]['near_repeats'][:5]),
                # 'group_hairpins': str(group_sec_struct[i]['hairpins'][:5]),
                # 'group_near_hairpins': str(group_sec_struct[i]['near_hairpins'][:5]),

                # # Well-averaged data
                # 'well_avg_group_strt': sum(strts) / len(strts),
                # 'well_avg_group_end': sum(ends) / len(ends),
                # 'well_avg_group_repeats': str(avg_group_sec_struct[i]['repeats'][:5]),
                # 'well_avg_group_near_repeats': str(avg_group_sec_struct[i]['near_repeats'][:5]),
                # 'well_avg_group_hairpins': str(avg_group_sec_struct[i]['hairpins'][:5]),
                # 'well_avg_group_near_hairpins': str(avg_group_sec_struct[i]['near_hairpins'][:5]),

            # TOOK OUT FRAGMENT CONCENTRATIONS BECAUSE TWO "FIELDS" HAVE VARYING AMOUNT OF DATA AND WILL THROW OFF COLUMN NAMES
            # tmp['fragment_concentrations'] = exact_del['fragment_concentrations']
            # frag_conc = ast.literal_eval(str(exact_del['fragment_concentrations']))
            # log2_frag_conc = []
            # for fc in frag_conc:
            #     if not fc:
            #         log2_frag_conc.append('')
            #     else:
            #         log2_frag_conc.append(log2(float(fc)))
            # tmp['Log2_fragment_concentrations'] = log2_frag_conc
            
            
            # Left side for dynamic programming
            tmp['dp_Left_Sequence'] = dpLeftSequence
            tmp['dp_Left_Energy'] = dpLeftEnergy
            lim5 = 0
            for k,v in dpLeftPathDict.items():
                tmp[f'{k}_LEFT'] = v
                if k[-1] == '2':
                    lim5 += 1
                    if lim5 == 5: break
            if lim5 < 5:
                for _ in range(lim5, 5):
                    tmp[f'length_{_}'] = 0
                    tmp[f'energy_{_}'] = 0
                    tmp[f'sequence_{_}'] = ''
                    tmp[f'gap{_}_d1'] = 0
                    tmp[f'gap{_}_d2'] = 0
            
            # Right side for dynamic programming
            tmp['dp_Right_Sequence'] = dpRightSequence
            tmp['dp_Right_Energy'] = dpRightEnergy
            lim5 = 0
            for k,v in dpRightPathDict.items():
                tmp[f'{k}_RIGHT'] = v
                if k[-1] == '2':
                    lim5 -= 1
                    if lim5 == 5: break
            if lim5 < 5:
                for _ in range(lim5, 5):
                    tmp[f'length_{_}'] = 0
                    tmp[f'energy_{_}'] = 0
                    tmp[f'sequence_{_}'] = ''
                    tmp[f'gap{_}_d1'] = 0
                    tmp[f'gap{_}_d2'] = 0
            
            
            dumpable_data[gene_id][str(exact_del['del_start']) + '-' + str(exact_del['del_end'])] = tmp
            with open('data/consolidated_top_deletions.json', 'a') as td:
                json.dump(dumpable_data, td)#save whole dict to save key sequence names
                td.write(',\n')



def main():
    frag_data = Path('data/customer_items_block_quant_alphaprod.csv')
    frag_id_2_frag_len = get_frag_lengths2(frag_data)

    # instead of getting the st1's from prelive, list them off from alpha prod using
    # boto3 commands and iterate through those
    response = s3.list_objects_v2(
        Bucket=buck,
        Prefix='static/elegen_csv/st1_csv/'
    )
    full_seq_id_2_frag_pass = {}
    for content in response.get('Contents', []):
        if content['Key'][-1] == '/' or content['Key'] == 'static/elegen_csv/st1_csv/SVT1-final-qc_18-26-51_08-17-22.csv':#this st1 was eliminated due to misnaming of columns
            continue
        print(content['Key'])
        obj = s3.get_object(Bucket=buck, Key=content['Key'])
        st1_content = pd.read_csv(obj['Body'])
        full_seq_id_2_frag_pass.update(get_number_failed_frags(st1_content, frag_id_2_frag_len))
    print(len(full_seq_id_2_frag_pass))


    st_final_file_locations = Path('data/batch_st_final_filename_alphaprod.csv')

    geneid_2_conf = {}
    conf_batches = set()
    missing_st_final = 0
    conf_batch_cnt = 0
    with st_final_file_locations.open('r') as rd:
        csv.field_size_limit(sys.maxsize)
        reader = csv.DictReader(rd)
        for row in reader:
            try:
                obj = s3.get_object(Bucket=buck, Key='static/elegen_csv/st_final_csv/' + row['st_final_filename'])
                batchname = row['batch_number'].partition('_')[0]
                date = row['date'].partition(' ')[0]
                conf_batch_cnt += 1
            except:
                missing_st_final += 1
                continue
            st_final_content = pd.read_csv(obj['Body'])
            geneid_2_conf.update(get_conf_data(st_final_content, batchname, date))
            conf_batches.add(batchname)
            print(batchname)

    print(len(geneid_2_conf))
    print("missing: ", missing_st_final)
    print('batch count: ', conf_batch_cnt)

    # Process the original json
    error_lines = 0
    # with open('data/top_deletions_2.3.25_finished.json') as f:
    # with open('data/top_deletions_5perc_group_2.17_finished.json') as f:
    with open('data/top_deletions_5perc_group.json') as f:
        combined_wells = []
        total_wells = 0
        gene_in_prev_well = '2078_5_B274-1'
        broken_datapoints = 0
        for line in f:
            try:
                tmp = ast.literal_eval(line)
            except:
                error_lines += 1
                continue
            for seq,seq_dels in tmp[0].items():
                gene_id = seq.rsplit('_', 1)[0]
                # process deletion groups across all wells for a sequence
                if gene_id != gene_in_prev_well and gene_in_prev_well:
                    # process_all_wells(
                    #     combined_wells,
                    #     total_wells,
                    #     gene_in_prev_well,
                    #     full_seq_id_2_frag_pass[gene_in_prev_well]['sequence'],
                    # )
                    try:
                        process_all_wells(
                            combined_wells,
                            total_wells,
                            gene_in_prev_well,
                            full_seq_id_2_frag_pass[gene_in_prev_well]['sequence'],
                        )
                    except KeyError:
                        print(seq)

                    combined_wells = []
                    total_wells = 0
                    gene_in_prev_well = gene_id

                # tmp_dels = reprocess_well(seq_dels, gene_id, full_seq_id_2_frag_pass[gene_id]['concentrations'])
                # total_wells += 1
                # if total_wells == 1:
                #     combined_wells = [tmp_dels]
                # else:
                #     combined_wells.append(tmp_dels)
                try:
                    tmp_dels = reprocess_well(seq_dels, gene_id, full_seq_id_2_frag_pass[gene_id]['concentrations'])
                    total_wells += 1
                    if total_wells == 1:
                        combined_wells = [tmp_dels]
                    else:
                        combined_wells.append(tmp_dels)
                except:
                    broken_datapoints += 1
                    print('broke: ', broken_datapoints, f' {seq}')
                    continue
    


if __name__ == '__main__':
    main()