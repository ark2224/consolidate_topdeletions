import logging
import numpy as np
import pysam
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

from sequence_manipulation import *
from data_retrieval import *

import speedy


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
            overall_reps = find_only_repeats(seq[ds-35:de+35])
            overall_near_reps = find_only_near_repeats(seq[ds-35:de+35])
            if len(overall_reps):
                sec_struct_2_keep = []
                top5limit = 5
                for strukt in overall_reps:
                    if (
                        (strukt['positions'][0] >= 5+lR)
                        or (strukt['positions'][1]+strukt['length'] <= de-ds+65-rL)
                        # and (strukt['ATGC_score'] > 0)
                        ):
                        sec_struct_2_keep.append({
                            'start': strukt['positions'][0] - 35 - lR,
                            'end': 35 + de - ds - rL - strukt['positions'][1]+strukt['length'],
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
    # seq2019_4 = 'TAATACGACTCACTATAGGGATAATACTGTGCGCTGATGAGTCCGTGAGGACGAAACCCGGCGTACCGGGTCCGCACAGTGGATCCTAGGCTATTGGATTGCGCTTTGCTTTTATCATTTTAGCAGATAGCCTCAGCTTTTTGTTGCGCGCATCCAGCACAACAATCTGGCGATGAGTGCTTCAAAGGAAATAAAGTCCTTTCTGTGGACACAATCATTAAGGAGGGAATTATCTAGTTACTGTTCCAACATCAAACTACAGGTGGTGAAAGATGCCCAGGCTCTTTTGCATGGGCTTGACTTCTCCGAAGTCAGCAATGTTCAACGGTTGATGCGCAAGGAGAGAAGGGATGACAATGATTTAAAGCGGTTGAGGGACCTTAATCAAGCGGTCAACAACCTTGTTGAGTTAAAATCAACTCAACAAAAGAGTATATTGAGAGTTGGGACCTTAACCTCAGATGACTTACTAATCTTAGCTGCTGATCTAGAGAAGCTGAAATCAAAGGTGACTAGGACAGAAAGGCCATTGAGTGCAGGTGTCTACATGGGTAACCTGAGCTCACAGCAACTTGACCAGAGACGAGCTCTTCTGAGCATGATAGGGATGAGTGGTAGTAATCAAGGGGCTCGGGCTGGGAGAGATGGAGTGGTGAGAGTTTGGGATGTGAAAAATGCAGAGTTGCTCAATAATCAGTTCGGGACCATGCCAAGTCTTACACTGGCATGTTTGACAAAACAGGGGCAGGTTGACTTGAATGATGCAGTGCAAGCATTGACAGATTTGGGTCTGATCTACACGGCAAAGTATCCCAATACTTCAGACCTGGATAGGCTGACTCAAAGTCACCCCATCCTAAACATGATTGACACCAAGAAGAGCTCTTTAAACATCTCAGGTTATAATTTTAGCTTGGGTGCAGCTGTGAAGGCAGGAGCTTGCATGCTGGATGGTGGCAACATGTTGGAGACAATCAAGGTGTCACCTCAGACAATGGATGGTATCCTCAAATCCATCTTGAAGGTCAAGAAGGCTCTTGGGATGTTCATCTCAGACACCCCTGGTGAAAGGAATCCTTATGAGAACATACTCTATAAGATCTGTCTGTCAGGAGATGGGTGGCCATACATTGCATCAAGGACCTCAATAACAGGAAGGGCTTGGGAAAACACTGTCGTTGATCTGGAATCAGATGGGAAGCCACAGAAAGCTGGCAGCAACAATTCTAACAAATCCCTGCAGTCGGCAGGGTTTACCGCTGGGCTTACCTATTCTCAGCTAATGACTCTCAAGGATGCAATGCTGCAACTTGATCCAAATGCTAAGACCTGGATGGATATTGAAGGAAGACCTGAAGATCCAGTGGAAGTTGCCCTCTATCAACCAAGTTCAGGCTGTTACATTCACTTCTTCCGTGAACCTACTGATTTAAAGCAGTTCAAACAGGATGCTAAGTACTCACATGGGATTGATGTCACAGACCTCTTCGCTGCACAACCGGGCCTGACCAGTGCTGTCATAGATGCACTCCCTCGGAATATGGTCATTACCTGTCAGGGGTCCGATGACATAAGGAAACTCCTTGAATCACAAGGTAGGAAAGACATCAAACTAATTGATATTGCCCTCAGCAAAACTGATTCCAGGAAGTATGAAAATGCAGTCTGGGATCAGTATAAAGACTTATGCCACATGCACACAGGTGTCGTTGTTGAGAAGAAGAAAAGAGGTGGTAAAGAGGAAATAACCCCTCACTGTGCACTAATGGACTGCATCATGTTTGACGCAGCAGTGTCAGGAGGACTGAACACATCGGTTTTGAGAGCAGTGCTACCCAGGGATATGGTGTTCAGAACATCAACACCTAGAGTCGTTCTGTAAATGGACGCCCCCGTGACCCACCGCCAACAGGCGGTGGGTCACGGGGGCCCTGACAAGGGTCTCATCTTTTCCATCTCACAGGCACACCGGGCTGTTTGTAGAGTCCACAGGAACATATGCCCATGTGATTCAATCTGTGAGGTTTAGGACACGACTTGCCTACAATATGCCTATGAGTTGGTATTTTGACTAGGTGAAGGAAGATGCTGATGAGATAGAAACTTGTGCTGAACACAAAGAGGTCAACTAGGCCCAATGGTGTCTTCCCCTGCCTCTCCATATACTCCTTCTGTAACATCTCAGTGATCATGTTGTCAGCTTGTTGTTCAATATCATCAGAAAAGTGGGTCTCGTTCAAGTATGACCCATTTGAGACAAGCCAACATTTAGGCAGTGATGTTTTCCCAGTAGTTGTGTGGTTGAGGTACCAATATTTGCTGTAATTACAGTATGGGATTCCCATGATGTCCCGTAGATGGTTCTTCATTATAAGTTGATCATTTATCAAGGCATTTACTGCCTTGTTGATCAACTGAATGCTCATTTGTGCTTCAGCTTTCAACCTTTGAATGGCTTGTTTGTTGAAGTCAAACAGCCTCAACATGTCACAAAATTCCTCATCATGCTTCTCATTACATTTTGCTACAGCTGTGTTCCCGAAGCATTTTAATTCAGCTTCAATCAGCATCCATCTGGTCAGACAATATCCCCCTGGTGTGTCCTTGCCTTCAGAATCTGACAGTGTCCATGTGAATGTGCCTAGTAGTCTTCTACTGATATAAATATCTCTGGTCCTTTGTGAGAGAAGCCCAAGATAACCGATAGGAGATGGTCTTGAGAATTGGCAGTGATCTTCCCAGGTTGTATTTTGGATTATCAGATATTGATAACTAGTCATAATACAGTCCCACTTGCCATGGCCTGAGTCAAGAGCAATGTAGCTCCCACCCCAAGCCATCCTCATGAAAGTCTGTAAGACACCATTTGCAACAGTGCCACAATGGTTGGCAGCATCCCCAGCATAGCTGTGACTCAGGTTGTACTGCACACTAATCTTTCCCCCATTGAAATCACAGCTCATTGCTTCATACTGATTGAAGTTGGGGATGGACAAGTGGAAAGTTGAAATTATGCTCATAAGGGCATGGTCATAGAGGTTCTTTTTGTGGGCATCAGACAGATTGCAAAATTTATGATTAATAATGCTTGTGTTGGTCAAGGTCAGTTCTAGTCCTGTCTCATTACCTACCATTATATAATGATGACTGTTATTCTTTGTGCAGGAGAGAGGCATGGTCATATTGAGTGTCTCCATGTTTAACTCCAGAGTCTGAAGCTCATAGACCCCTTTATAAAGACTGGTTGTGCAAGACCTGCCACACAACAGGAGGAAAGTAACCAAACCAACGAGGCCACATGTTGCAAAGTTGTACAGACCTTTCAGCACTGCTAATATAGACAGCGCAATGAGGACAATGTTCATCACCTCTTCTATTACATGAGGCACTTCCTGGAAGAATGTCACTATCTGTCCCATTTTGAATAGGACACTTGAATTGCGCAACCAAAAATGCCTAGGATCCCCGGTGCGGGGTCGGCATGGCATCTCCACCTCCTCGCGGTCCGACCTGGGCATCCGAAGGAGGACGTCGTCCACTCGGATGGCTAAGGGAGAGCTCGGATCCGGCTGCTAACAAAGCCCGAAAGGAAGCTGAGTTGGCTGCTGCCACCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTAC'
    # hairpins = find_hairpins(seq2019_4)
    # print(hairpins)
    # print('Hairpins^')
    # nearhairpins = find_near_hairpins(seq2019_4)
    # nearhairpins = find_near_hairpins('')
    # print('Near Hairpins^')
    # print(nearhairpins)
    # repeats = find_repeats('GGATCTGGCGAGGGAGGATCTTCTA')
    # repeats = find_repeats(seq2019_4)
    # print(repeats)
    # return
    # print(repeats)
    # print('Repeats^')
    # nearrepeats = find_near_repeats(seq2019_4)
    # print(nearrepeats)
    # print('Near Repeats^')
    # return
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
    with open('data/top_deletions_2.3.25_finished.json') as f:
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
    
    # Convert into csv
    error_lines = 0
    with open('data/consolidated_top_deletions.json') as f:
        with open('data/consolidated_top_deletions_csv.csv', 'a') as outf:
            writer = csv.writer(outf)
            writer.writerow([
                'Sequence ID',
                'Coordinates',

                'Batch Number',
                'Total Wells',
                
                'Deletion Start',
                'Deletion End',
                'Deletion Appears in X Wells',
                'Similar to Past Deletion',
                'Read Count',
                'Percentage of Reads',
                'Matching Over Deletion',
                'Right Energy',
                'Right Offset',
                'Left Energy',
                'Left Offset',
                'Longest Repeat',
                'Longest Repeat - Distance from 5\' Deletion-End',
                'Longest Repeat - Distance from 3\' Deletion-End',
                'Longest Repeat Beginning Del-after Del',
                'Distance after Deletion Start',
                'Distance after Deletion End',
                'Longest Repeat Before Del-End Del',
                'Distance Before Deletion Start',
                'Distance Before Deletion End',
                'Longest Repeat AT-GC Score',
                'Longest A-B Ligation Repeat',
                'Within One Fragment',
                'Distance from Del Beginning to Fragment Overhang',
                'Distance from Del End to Fragment Overhang',
                'Perfect Sequence Across Deletion',
                'Energy Longest Repeat Beginning Del-after Del',
                'Energy Longest Repeat Before Del-End Del',
                'Log2 Exact Deletion Frequency',
                # Exact deletion secondary structures
                'IHFL_hairpin_start-Optimal',
                'IHFL_hairpin_end-Optimal',
                'IHFL_hairpin_stem-Optimal',
                'IHFL_hairpin_loop-Optimal',
                'IHFL_hairpin_energy-Optimal',

                'IHFL_hairpin_start-Sub_Optimal',
                'IHFL_hairpin_end-Sub_Optimal',
                'IHFL_hairpin_stem-Sub_Optimal',
                'IHFL_hairpin_loop-Sub_Optimal',
                'IHFL_hairpin_energy-Sub_Optimal',

                'IHFL_near_hairpin_start-Optimal',
                'IHFL_near_hairpin_end-Optimal',
                'IHFL_near_hairpin_stem-Optimal',
                'IHFL_near_hairpin_loop-Optimal',
                'IHFL_near_hairpin_energy-Optimal',
                'IHFL_near_hairpin_mismatch_count-Optimal',

                'IHFL_near_hairpin_start-Sub_Optimal',
                'IHFL_near_hairpin_end-Sub_Optimal',
                'IHFL_near_hairpin_stem-Sub_Optimal',
                'IHFL_near_hairpin_loop-Sub_Optimal',
                'IHFL_near_hairpin_energy-Sub_Optimal',
                'IHFL_near_hairpin_mismatch_count-Sub_Optimal',


                'PDFL_hairpin_start-Optimal',
                'PDFL_hairpin_end-Optimal',
                'PDFL_hairpin_stem-Optimal',
                'PDFL_hairpin_loop-Optimal',
                'PDFL_hairpin_energy-Optimal',

                'PDFL_hairpin_start-Sub_Optimal',
                'PDFL_hairpin_end-Sub_Optimal',
                'PDFL_hairpin_stem-Sub_Optimal',
                'PDFL_hairpin_loop-Sub_Optimal',
                'PDFL_hairpin_energy-Sub_Optimal',

                'PDFL_near_hairpin_start-Optimal',
                'PDFL_near_hairpin_end-Optimal',
                'PDFL_near_hairpin_stem-Optimal',
                'PDFL_near_hairpin_loop-Optimal',
                'PDFL_near_hairpin_energy-Optimal',
                'PDFL_near_hairpin_mismatch_count-Optimal',

                'PDFL_near_hairpin_start-Sub_Optimal',
                'PDFL_near_hairpin_end-Sub_Optimal',
                'PDFL_near_hairpin_stem-Sub_Optimal',
                'PDFL_near_hairpin_loop-Sub_Optimal',
                'PDFL_near_hairpin_energy-Sub_Optimal',
                'PDFL_near_hairpin_mismatch_count-Sub_Optimal',


                'IHFR_hairpin_start-Optimal',
                'IHFR_hairpin_end-Optimal',
                'IHFR_hairpin_stem-Optimal',
                'IHFR_hairpin_loop-Optimal',
                'IHFR_hairpin_energy-Optimal',

                'IHFR_hairpin_start-Sub_Optimal',
                'IHFR_hairpin_end-Sub_Optimal',
                'IHFR_hairpin_stem-Sub_Optimal',
                'IHFR_hairpin_loop-Sub_Optimal',
                'IHFR_hairpin_energy-Sub_Optimal',

                'IHFR_near_hairpin_start-Optimal',
                'IHFR_near_hairpin_end-Optimal',
                'IHFR_near_hairpin_stem-Optimal',
                'IHFR_near_hairpin_loop-Optimal',
                'IHFR_near_hairpin_energy-Optimal',
                'IHFR_near_hairpin_mismatch_count-Optimal',

                'IHFR_near_hairpin_start-Sub_Optimal',
                'IHFR_near_hairpin_end-Sub_Optimal',
                'IHFR_near_hairpin_stem-Sub_Optimal',
                'IHFR_near_hairpin_loop-Sub_Optimal',
                'IHFR_near_hairpin_energy-Sub_Optimal',
                'IHFR_near_hairpin_mismatch_count-Sub_Optimal',


                'PDFR_hairpin_start-Optimal',
                'PDFR_hairpin_end-Optimal',
                'PDFR_hairpin_stem-Optimal',
                'PDFR_hairpin_loop-Optimal',
                'PDFR_hairpin_energy-Optimal',

                'PDFR_hairpin_start-Sub_Optimal',
                'PDFR_hairpin_end-Sub_Optimal',
                'PDFR_hairpin_stem-Sub_Optimal',
                'PDFR_hairpin_loop-Sub_Optimal',
                'PDFR_hairpin_energy-Sub_Optimal',

                'PDFR_near_hairpin_start-Optimal',
                'PDFR_near_hairpin_end-Optimal',
                'PDFR_near_hairpin_stem-Optimal',
                'PDFR_near_hairpin_loop-Optimal',
                'PDFR_near_hairpin_energy-Optimal',
                'PDFR_near_hairpin_mismatch_count-Optimal',

                'PDFR_near_hairpin_start-Sub_Optimal',
                'PDFR_near_hairpin_end-Sub_Optimal',
                'PDFR_near_hairpin_stem-Sub_Optimal',
                'PDFR_near_hairpin_loop-Sub_Optimal',
                'PDFR_near_hairpin_energy-Sub_Optimal',
                'PDFR_near_hairpin_mismatch_count-Sub_Optimal',


                'Overall_hairpin_start-1',
                'Overall_hairpin_end-1',
                'Overall_hairpin_stem-1',
                'Overall_hairpin_loop-1',
                'Overall_hairpin_energy-1',

                'Overall_hairpin_start-2',
                'Overall_hairpin_end-2',
                'Overall_hairpin_stem-2',
                'Overall_hairpin_loop-2',
                'Overall_hairpin_energy-2',


                'Overall_hairpin_start-3',
                'Overall_hairpin_end-3',
                'Overall_hairpin_stem-3',
                'Overall_hairpin_loop-3',
                'Overall_hairpin_energy-3',

                'Overall_hairpin_start-4',
                'Overall_hairpin_end-4',
                'Overall_hairpin_stem-4',
                'Overall_hairpin_loop-4',
                'Overall_hairpin_energy-4',

                'Overall_hairpin_start-5',
                'Overall_hairpin_end-5',
                'Overall_hairpin_stem-5',
                'Overall_hairpin_loop-5',
                'Overall_hairpin_energy-5',


                'Overall_near_hairpin_start-1',
                'Overall_near_hairpin_end-1',
                'Overall_near_hairpin_stem-1',
                'Overall_near_hairpin_loop-1',
                'Overall_near_hairpin_energy-1',
                'Overall_near_hairpin_mismatch_count-1',

                'Overall_near_hairpin_start-2',
                'Overall_near_hairpin_end-2',
                'Overall_near_hairpin_stem-2',
                'Overall_near_hairpin_loop-2',
                'Overall_near_hairpin_energy-2',
                'Overall_near_hairpin_mismatch_count-2',

                'Overall_near_hairpin_start-3',
                'Overall_near_hairpin_end-3',
                'Overall_near_hairpin_stem-3',
                'Overall_near_hairpin_loop-3',
                'Overall_near_hairpin_energy-3',
                'Overall_near_hairpin_mismatch_count-3',

                'Overall_near_hairpin_start-4',
                'Overall_near_hairpin_end-4',
                'Overall_near_hairpin_stem-4',
                'Overall_near_hairpin_loop-4',
                'Overall_near_hairpin_energy-4',
                'Overall_near_hairpin_mismatch_count-4',

                'Overall_near_hairpin_start-5',
                'Overall_near_hairpin_end-5',
                'Overall_near_hairpin_stem-5',
                'Overall_near_hairpin_loop-5',
                'Overall_near_hairpin_energy-5',
                'Overall_near_hairpin_mismatch_count-5',



                'Overall_repeat_start-1',
                'Overall_repeat_end-1',
                'Overall_repeat_length-1',
                'Overall_repeat_energy-1',

                'Overall_repeat_start-2',
                'Overall_repeat_end-2',
                'Overall_repeat_length-2',
                'Overall_repeat_energy-2',

                'Overall_repeat_start-3',
                'Overall_repeat_end-3',
                'Overall_repeat_length-3',
                'Overall_repeat_energy-3',

                'Overall_repeat_start-4',
                'Overall_repeat_end-4',
                'Overall_repeat_length-4',
                'Overall_repeat_energy-4',

                'Overall_repeat_start-5',
                'Overall_repeat_end-5',
                'Overall_repeat_length-5',
                'Overall_repeat_energy-5',

                'Overall_near_repeat_start-1',
                'Overall_near_repeat_end-1',
                'Overall_near_repeat_length-1',
                'Overall_near_repeat_energy-1',
                'Overall_near_repeat_mismatch_count-1',

                'Overall_near_repeat_start-2',
                'Overall_near_repeat_end-2',
                'Overall_near_repeat_length-2',
                'Overall_near_repeat_energy-2',
                'Overall_near_repeat_mismatch_count-2',

                'Overall_near_repeat_start-3',
                'Overall_near_repeat_end-3',
                'Overall_near_repeat_length-3',
                'Overall_near_repeat_energy-3',
                'Overall_near_repeat_mismatch_count-3',

                'Overall_near_repeat_start-4',
                'Overall_near_repeat_end-4',
                'Overall_near_repeat_length-4',
                'Overall_near_repeat_energy-4',
                'Overall_near_repeat_mismatch_count-4',

                'Overall_near_repeat_start-5',
                'Overall_near_repeat_end-5',
                'Overall_near_repeat_length-5',
                'Overall_near_repeat_energy-5',
                'Overall_near_repeat_mismatch_count-5',
                



                # 'Deletion Group Start',
                # 'Deletion Group End',
                # 'Log2 Deletion Group Frequency',
                # 'Log2 Not in Spanning Del Group Frequency',
                # 'Furthest Right',
                # 'Furthest Left',
                # 'Group Repeats',
                # 'Group Near Repeats',
                # 'Group Hairpins',
                # 'Group Near Hairpins',
                
                # 'Group Start - Well Avg',
                # 'Group End - Well Avg',
                # 'Group Repeats - Well Avg',
                # 'Group Near Repeats - Well Avg',
                # 'Group Hairpins - Well Avg',
                # 'Group Near Hairpins - Well Avg',

                'Fragment Concentrations',
            ])
            for line in f:
                # print(line[1:9])
                # if line[1:9] == '"2078_5"':
                #     print(line)
                # print(line[1:20])
                # try:
                #     tmp = ast.literal_eval(line)
                # except:
                #     error_lines += 1
                #     continue
                tmp = ast.literal_eval(line)
                if not len(tmp[0]):
                    print('here')
                    continue
                for seq,seq_dels in tmp[0].items():
                    # conf_idx = seq[:-2]
                    # if conf_idx not in geneid_2_conf:
                    #     passing_dial_out = False
                    # else:
                    #     passing_dial_out = geneid_2_conf[conf_idx]['wells_passing'] > 0

                    for coordinates,data_fields in seq_dels.items():
                        tmprow = [seq, coordinates]
                        for k,v in data_fields.items():

                            tmprow.append(v)
                        # tmprow.append(passing_dial_out)
                        writer.writerow(tmprow)
    print(error_lines)


if __name__ == '__main__':
    main()