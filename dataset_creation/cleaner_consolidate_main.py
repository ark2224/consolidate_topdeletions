import numpy as np
import csv

from pathlib import Path
from math import log2
import os
import ast

import pandas as pd
import sys
import json

from cred import (
    s3,
    buck
)
from data_retrieval import (
    get_frag_lengths2,
    get_number_failed_frags,
)
from dp_path_finder import (
    find_graph,
    process_output_path
)
from sequence_manipulation import (
    longest_common_substring,
    ATGC_energy,
    find_mirrors,
    fix_energy_score,
    find_only_hairpins,
    find_only_repeats,
    longest_substring,
    find_complements
)


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


def process_all_wells(cw: dict, well_cnt: int, gene_batch_id: str, seq: str, element, overlaps):
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
                data_fields['appears_in_x_wells'] += 1  #Still VERY MUCH WRONG
                consolidated_deletions[coordinates]['Log2_exact_deletion'] += data_fields['Log2_exact_deletion']
                consolidated_deletions[coordinates]['perfect_sequence_across_deletion'] |= (data_fields['perfect_sequence_across_deletion'].lower() == 'true')
                frac = data_fields['percentage_of_reads']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(frac/(1-frac))'] += 10 + log2(max(2**-15, frac / (1 - frac))) if frac != 1 else 2**63 - 1
                group_frac = data_fields['group_percent']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(group_frac/(1-group_frac))'] += 10 + log2(max(2 **-15, group_frac / (1 - group_frac))) if group_frac != 1 else 2**63 - 1
                consolidated_deletions[coordinates]['Log2_deletion_group_percent_of_reads'] += max(log2(group_frac), -15)
            else:
                consolidated_deletions[coordinates] = data_fields
                consolidated_deletions[coordinates]['perfect_sequence_across_deletion'] = (data_fields['perfect_sequence_across_deletion'].lower() == 'true')
                frac = data_fields['percentage_of_reads']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(frac/(1-frac))'] = 10 + log2(max(2**-15, frac / (1 - frac))) if frac != 1 else 2**63 - 1
                group_frac = data_fields['group_percent']
                consolidated_deletions[coordinates]['Avg_wells_10+log2(group_frac/(1-group_frac))'] = 10 + log2(max(2**-15, group_frac / (1 - group_frac))) if group_frac != 1 else 2**63 - 1
                consolidated_deletions[coordinates]['Log2_deletion_group_percent_of_reads'] = max(log2(group_frac), -15)
            if group_frac < 0.05:
                print(' CHECK ERROR IN GROUP FRAC RECORDING')
    
    # Establish average group deletions across combined wells' deletions
    group_del_groups = []#[[group_strt, group_end, log2(frequency), well_ct_with_group ]]
    for coords,cd_data in consolidated_deletions.items():
        cd_data['Log2_exact_deletion'] /= cd_data['appears_in_x_wells']
        cd_data['Avg_wells_10+log2(frac/(1-frac))'] /= cd_data['appears_in_x_wells']
        cd_data['Avg_wells_10+log2(group_frac/(1-group_frac))'] /= cd_data['appears_in_x_wells']
        cd_data['Log2_deletion_group_percent_of_reads'] /= cd_data['appears_in_x_wells']

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


    # Now re-calculate top10 deletions + rightmost/leftmost exact-deletion for each group
    gene_id, batch_num = gene_batch_id.rsplit('_', 1)
    # dumpable_data = {gene_id: {}}
    # for i, (strts, ends, log2_group_percent, exact_del_data, rightmost_strt, leftmost_end, frequency) in enumerate(group_del_groups):
    for i, (strts, ends, log2_group_percent, exact_del_data, _) in enumerate(group_del_groups):
        exact_del_data = sorted(exact_del_data, key=lambda x: x['Log2_exact_deletion'], reverse=True)

        # TOP 10 DELETIONS
        for exact_del in exact_del_data:#[:1]:#[:10]:            
            dumpable_data = {gene_id: {}}
            # Exact-Del Secondary Structures
            fiveprime_end = (
                exact_del['bases_at_beginning_of_deletion'][:50] +
                exact_del['bases_at_beginning_of_deletion'][51:]
            )
            threeprime_end = (
                exact_del['bases_at_end_of_deletion'][:50] + 
                exact_del['bases_at_end_of_deletion'][51:]
            )
            
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

            # === Dynamic Programming ===
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

            # === Secondary Structures =================================================================================
            blank_hairpin = [{
                'start': '',
                'end': '',
                'stem': '',
                'loop': '',
                'energy': '',
                'sequence': '',
            }]
            blank_repeat = [{
                'start': '',
                'end': '',
                'length': '',
                'energy': '',
                'sequence': '',
            }]
            blank_near_hairpin = blank_hairpin.copy()
            blank_near_hairpin[0]['mismatch_count'] = ''
            blank_near_repeat = blank_repeat.copy()
            blank_near_repeat[0]['mismatch_count'] = ''
            blank_complement = blank_repeat.copy()

            # ===== +/-100bp around LL / 5-PRIME =====
            # 5 prime hairpins
            fiveprime_srch_window = seq[ds-lL-120:ds-lL+120]
            lL100_hairpins = find_only_hairpins(
                fiveprime_srch_window,
                120 - lL,
                120 - lL,
                1,
                1,
                5
            )
            exact_del['5_prime_hairpins'] = lL100_hairpins['hairpins']
            exact_del['5_prime_near_hairpins'] = lL100_hairpins['near_hairpins']
            # 5 prime repeats
            lL100_repeats = find_only_repeats(
                fiveprime_srch_window,
                120 - lL,
                120 - lL,
                1,
                1,
                5
            )
            exact_del['5_prime_repeats'] = lL100_repeats['repeats']
            exact_del['5_prime_near_repeats'] = lL100_repeats['near_repeats']
            # 5 prime mirrors
            exact_del['5_prime_mirrors'] = find_mirrors(
                fiveprime_srch_window,
                120 - lL,
                120 - lL,
                1,
                1,
                5
            )
            # 5 prime complements
            exact_del['5_prime_complements'] = find_complements(
                fiveprime_srch_window,
                120-lL,
                120-lL,
                1,
                1,
                5
            )

            # ===== +/-100bp around RL /  3-PRIME =====
            # 3 prime hairpins
            threeprime_srch_window = seq[de-rL-120:de-rL+120]
            rl100_hairpins = find_only_hairpins(
                threeprime_srch_window,
                120 - rL,
                120 - rL,
                -1,
                -1,
                5
            )
            exact_del['3_prime_hairpins'] = rl100_hairpins['hairpins']
            exact_del['3_prime_near_hairpins'] = rl100_hairpins['near_hairpins']
            # 3 prime repeats
            rl100_repeats = find_only_repeats(
                threeprime_srch_window,
                120 - rL,
                120 - rL,
                -1,
                -1,
                5
            )
            exact_del['3_prime_repeats'] = rl100_repeats['repeats']
            exact_del['3_prime_near_repeats'] = rl100_repeats['near_repeats']
            # 3 prime mirrors
            exact_del['3_prime_mirrors'] = find_mirrors(
                fiveprime_srch_window,
                120 - rL,
                120 - rL,
                -1,
                -1,
                5
            )
            # 3 prime complements
            exact_del['3_prime_complements'] = find_complements(
                threeprime_srch_window,
                120 - rL,
                120 - rL,
                -1,
                -1,
                5
            )

            # OVERALL SECONDARY STRUCTURES
            overall_srch_window = seq[fs:fe]
            overall_hairpins = find_only_hairpins(
                overall_srch_window,
                ds-fs-lL,
                de-fs-rL,
                1,
                -1
            )
            exact_del['Overall_hairpins'] = overall_hairpins['hairpins'][:5]
            exact_del['Overall_near_hairpins'] = overall_hairpins['near_hairpins'][:5]
            # special overall hairpins - long-range
            sec_struct_2_keep_left = {}
            sec_struct_2_keep_right = {}
            sec_struct_2_keep_both = {}
            sec_struct_2_keep_before5p = {}
            sec_struct_2_keep_after5p = {}
            for k,v in overall_hairpins.items():
                sec_struct_2_keep_left[k] = []
                sec_struct_2_keep_right[k] = []
                sec_struct_2_keep_both[k] = []
                sec_struct_2_keep_before5p[k] = []
                sec_struct_2_keep_after5p[k] = []
                for strukt in v:
                    if not strukt['start']:
                        continue
                    s5 = strukt['start']
                    e3 = strukt['end']
                    s3 = -(s5 + ds - de)
                    e5 = de - e3 - ds
                    l = strukt['stem']
                    
                    if len(sec_struct_2_keep_both[k]) < 2 and (s5 <= 0 <= s5+l) and (0 <= e3+l and e3 <= 40):
                        sec_struct_2_keep_both[k].append(strukt)
                    elif len(sec_struct_2_keep_left[k]) < 2 and (s5 <= 0 <= s5+l or e5-l <= 0 <= e5):
                        sec_struct_2_keep_left[k].append(strukt)
                    elif len(sec_struct_2_keep_right[k]) < 2 and ((0 <= e3+l <= 40+l) or (-l <= s3-l <= 40)):
                        sec_struct_2_keep_right[k].append(strukt)

                    if len(sec_struct_2_keep_before5p[k]) < 5 and s5 <= 0 and e3 >= s5:
                        sec_struct_2_keep_before5p[k].append(strukt)
                    elif len(sec_struct_2_keep_after5p[k]) < 5 and e3 < s5 and e3 <= 1:
                        sec_struct_2_keep_after5p[k].append(strukt)

            exact_del.update({
                'long_range_straddling_5_prime_hairpins': sec_struct_2_keep_left['hairpins'],
                'long_range_straddling_5_prime_near_hairpins': sec_struct_2_keep_left['near_hairpins'],
                'long_range_3_prime_hairpins': sec_struct_2_keep_right['hairpins'],
                'long_range_3_prime_near_hairpins': sec_struct_2_keep_right['near_hairpins'],
                'long_range_5_3_prime_hairpins': sec_struct_2_keep_both['hairpins'],
                'long_range_5_3_prime_near_hairpins': sec_struct_2_keep_both['near_hairpins'],
                'long_range_before_5prime_hairpins': sec_struct_2_keep_before5p['hairpins'],
                'long_range_before_5prime_near_hairpins': sec_struct_2_keep_before5p['near_hairpins'],
                'long_range_after_5prime_hairpins': sec_struct_2_keep_after5p['hairpins'],
                'long_range_after_5prime_near_hairpins': sec_struct_2_keep_after5p['near_hairpins'],
            })
            # else:
            #     exact_del['long_range_straddling_5_prime_hairpins'] = blank_hairpin.copy()
            #     exact_del['long_range_straddling_5_prime_near_hairpins'] = blank_near_hairpin.copy()
            #     exact_del['long_range_3_prime_hairpins'] = blank_hairpin.copy()
            #     exact_del['long_range_3_prime_near_hairpins'] = blank_near_hairpin.copy()
            #     exact_del['long_range_5_3_prime_hairpins'] = blank_hairpin.copy()
            #     exact_del['long_range_5_3_prime_near_hairpins'] = blank_near_hairpin.copy()
            
            # overall repeats
            overall_reps = find_only_repeats(
                overall_srch_window,
                ds-fs-lL,
                de-fs-rL,
                1,
                -1,
                5
            )
            exact_del['Overall_repeats'] = overall_reps['repeats']
            exact_del['Overall_near_repeats'] = overall_reps['near_repeats']
            # overall palindromes
            exact_del['Overall_mirrors'] = find_mirrors(
                overall_srch_window,
                ds-fs-lL,
                de-fs-rL,
                1,
                -1,
                5
            )
            # overall complements
            exact_del['Overall_complements'] = find_complements(
                overall_srch_window,
                ds-fs-lL,
                de-fs-rL,
                1,
                -1,
                5
            )
            
            # # checking for secondary structures in fragment ends if deletion is over 170bp away from ends
            exact_del['fragment_start_hairpins'] = blank_hairpin.copy()
            exact_del['fragment_start_near_hairpins'] = blank_near_hairpin.copy()
            exact_del['fragment_start_repeats'] = blank_repeat.copy()
            exact_del['fragment_start_near_repeats'] = blank_near_repeat.copy()
            exact_del['fragment_start_mirrors'] = blank_repeat.copy()
            exact_del['fragment_start_complements'] = blank_complement.copy()
            if ds - fs >= 80:
                fragment_start_window = seq[fs:fs+100]
                fragment_start_hairpins = find_only_hairpins(
                    fragment_start_window,
                    ds-lL-fs,
                    ds-lL-fs,
                    1,
                    1,
                    5
                )
                exact_del['fragment_start_hairpins'] = fragment_start_hairpins['hairpins']
                exact_del['fragment_start_near_hairpins'] = fragment_start_hairpins['near_hairpins']

                fragment_start_reps = find_only_repeats(
                    fragment_start_window,
                    ds-lL-fs,
                    ds-lL-fs,
                    1,
                    1,
                    5
                )
                exact_del['fragment_start_repeats'] = fragment_start_reps['repeats']
                exact_del['fragment_start_near_repeats'] = fragment_start_reps['near_repeats']

                exact_del['fragment_start_mirrors'] = find_mirrors(
                    fragment_start_window,
                    ds-lL-fs,
                    ds-lL-fs,
                    1,
                    1,
                    5
                )

                exact_del['fragment_start_complements'] = find_complements(
                    fragment_start_window,
                    ds-lL-fs,
                    ds-lL-fs,
                    1,
                    1,
                    5
                )
                    
            exact_del['fragment_end_hairpins'] = blank_hairpin.copy()
            exact_del['fragment_end_near_hairpins'] = blank_near_hairpin.copy()
            exact_del['fragment_end_repeats'] = blank_repeat.copy()
            exact_del['fragment_end_near_repeats'] = blank_near_repeat.copy()
            exact_del['fragment_end_mirrors'] = blank_repeat.copy()
            exact_del['fragment_end_complements'] = blank_complement.copy()
            if fe - de >= 80:
                fragment_end_window = seq[fe-100:fe]
                fragment_end_hairpins = find_only_hairpins(
                    fragment_end_window,
                    100+de-rL-fe,
                    100+de-rL-fe,
                    -1,
                    -1,
                    5
                )
                exact_del['fragment_end_hairpins'] = fragment_end_hairpins['hairpins']
                exact_del['fragment_end_near_hairpins'] = fragment_end_hairpins['near_hairpins']

                fragment_end_reps = find_only_repeats(
                    fragment_end_window,
                    100+de-rL-fe,
                    100+de-rL-fe,
                    -1,
                    -1,
                    5
                )
                exact_del['fragment_end_repeats'] = fragment_end_reps['repeats']
                exact_del['fragment_end_near_repeats'] = fragment_end_reps['near_repeats']

                exact_del['fragment_end_mirrors'] = find_mirrors(
                    fragment_end_window,
                    100+de-rL-fe,
                    100+de-rL-fe,
                    -1,
                    -1,
                    5
                )
                exact_del['fragment_end_complements'] = find_complements(
                    fragment_end_window,
                    100+de-rL-fe,
                    100+de-rL-fe,
                    -1,
                    -1,
                    5
                )


            # Fixing right/left energy+offset
            (
                right_energy,
                right_offset,
                left_energy,
                left_offset
            ) = fix_energy_score(exact_del['matching_over_deletion'])
            
            # fixing longest_repeat offset
            long_rep = longest_substring(fiveprime_end, threeprime_end)
            exact_del['longest_repeat'] = long_rep
            fiveprime_offset = (len(long_rep) + fiveprime_end.find(long_rep)) - (50 - lL)
            threeprime_offset = (49 - rL) - (threeprime_end.find(long_rep) + len(long_rep))

            after_del = threeprime_end[50-rL:]
            beginning_del = fiveprime_end[50+lR:]
            longest_repeat_beginning_del_after_del = longest_substring(beginning_del, after_del)
            repeat_dist_beg_del = (
                beginning_del.find(longest_repeat_beginning_del_after_del) +
                len(longest_repeat_beginning_del_after_del)
            )#can be zero-indexed
            repeat_dist_after_del = -(1 +
                                      after_del.find(longest_repeat_beginning_del_after_del) +
                                      len(longest_repeat_beginning_del_after_del)
                                      )

            before_del = fiveprime_end[:50+lR]
            end_del = threeprime_end[:50-rL]
            longest_repeat_before_del_end_del = longest_substring(before_del, end_del)
            repeat_dist_before_del = 50 + lR - (len(longest_repeat_before_del_end_del) + before_del.find(longest_repeat_before_del_end_del))
            repeat_dist_end_del = 50 - rL - (1 + len(longest_repeat_before_del_end_del) + end_del.find(longest_repeat_before_del_end_del))#can be zero-indexed

            # if rightmost_strt == exact_del['del_start']:
            #     exact_del['furthest_right'] = True
            # if leftmost_end == exact_del['del_end']:
            #     exact_del['furthest_left'] = True
            is_element = 0
            for i in range(1, len(overlaps)):
                if overlaps[i-1][0] <= ds and de <= overlaps[i][1]:
                    is_element = element[i-1]

            tmp = {
                # Sequence meta-data:
                'batch_num': batch_num,
                'total_wells': well_cnt,
                'Log2 Total Reads in Well': exact_del['log2_total_reads_well'],

                # Exact-Deletion Data
                'del_start': exact_del['del_start'],
                'del_end': exact_del['del_end'],
                'del_length': exact_del['del_end'] - exact_del['del_start'],
                'fragment_is_element': is_element,
                'longest_repeat': long_rep,
                'longest_repeat_dist_5prime': fiveprime_offset,
                'longest_repeat_dist_3prime': threeprime_offset,


                'Avg_wells_10+log2(frac/(1-frac))': exact_del['Avg_wells_10+log2(frac/(1-frac))'],
                'Avg_wells_10+log2(group_frac/(1-group_frac))': exact_del['Avg_wells_10+log2(group_frac/(1-group_frac))'],
                'Left_End_Score': exact_del['left_end_score'],
                'Right_End_Score': exact_del['right_end_score'],
                'LL': exact_del['LL'],
                'LR': exact_del['LR'],
                'RL': exact_del['RL'],
                'RR': exact_del['RR'],

                '10bp_before_LL_deletion': fiveprime_end[40-lL:50-lL],
                '10bp_after_LL_deletion_start': fiveprime_end[50-lL:60-lL],
                '10bp_before_RL_deletion_end': threeprime_end[40-rL:50-rL],
                '10bp_after_RL_deletion_end': threeprime_end[50-rL:60-rL],
                # 'del_in_x_wells': exact_del['appears_in_x_wells'],
                # 'similar_to_past_deletion': exact_del['similar_to_past_deletion'],
            }

            base_scoring = {
                'A': 0,
                'C': 1,
                'G': 2,
                'T': 3
            }
            for name, side in [('right_of_LL', fiveprime_end[50-lL:60-lL]), ('left_of_RL', threeprime_end[50-rL:60-rL])]:
                for b in range(10):
                    if b >= len(side) or side[b] not in base_scoring:
                        tmp[f'{b}_{name}'] = -1
                    else:
                        tmp[f'{b}_{name}'] = base_scoring[side[b]]
            tmp.update({
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
                'fragment_gc_content': exact_del['fragment_gc_content'],

                'perfect_sequence_across_deletion': str(exact_del['perfect_sequence_across_deletion']),
                'energy_longest_repeat_beginning_del_after_del': exact_del['energy_longest_repeat_beginning_del_after_del'],
                'energy_longest_repeat_before_del_end_del': exact_del['energy_longest_repeat_before_del_end_del'],
                'Log2_exact_deletion_frequency': exact_del['Log2_exact_deletion'],
                # 'Log2_group_deletion_frequency': log2_group_percent#exact_del['Log2_deletion_group_frequency'],
            })

            # Secondary Structures used for Heatmaps:
            for cat in [
                '5_prime',
                '3_prime',
                'Overall',
                'fragment_start',
                'fragment_end'
            ]:
                for s in [
                    'hairpins',
                    'near_hairpins',
                    'repeats',
                    'near_repeats',
                    'mirrors',
                    'complements'
                ]:
                    while len(exact_del[f'{cat}_{s}']) < 5:
                        if s == 'hairpins' or s == 'near_hairpins':
                            blank = {
                                'start': '',
                                'end': '',
                                'stem': '',
                                'loop': '',
                                'energy':'',
                                'sequence':'',
                            }
                        else:
                            blank = {
                                'start': '',
                                'end': '',
                                'length': '',
                                'energy':'',
                                'sequence':'',
                            }
                        if s[:4] == 'near':
                            blank['mismatch_count'] = '' 
                        exact_del[f'{cat}_{s}'].append(blank)
                    for b,c in enumerate(exact_del[f'{cat}_{s}']):
                        if s == 'hairpins' or s == 'near_hairpins':
                            tmp.update({
                                f'{cat}_{s}_start-{b}': c['start'],
                                f'{cat}_{s}_end-{b}': c['end'],
                                f'{cat}_{s}_stem-{b}': c['stem'],
                                f'{cat}_{s}_loop-{b}': c['loop'],
                                f'{cat}_{s}_energy-{b}':c['energy'],
                                f'{cat}_{s}_sequence-{b}':c['sequence'],
                            })
                        else:
                            tmp.update({
                                f'{cat}_{s}_start-{b}': c['start'],
                                f'{cat}_{s}_end-{b}': c['end'],
                                f'{cat}_{s}_length-{b}': c['length'],
                                f'{cat}_{s}_energy-{b}':c['energy'],
                                f'{cat}_{s}_sequence-{b}':c['sequence'],
                            })
                        if s[:4] == 'near':
                            tmp[f'{cat}_{s}_mismatch_count{b}'] = c['mismatch_count']
            # Misc Secondary Structures
            for cat in [
                'long_range_straddling_5_prime',
                'long_range_3_prime',
                'long_range_5_3_prime',
            ]:
                for s in ['hairpins', 'near_hairpins']:
                    while len(exact_del[f'{cat}_{s}']) < 2:
                        blank = {
                            'start': '',
                            'end': '',
                            'stem': '',
                            'loop': '',
                            'energy': '',
                            'sequence': '',
                        }
                        if s[:4] == 'near':
                            blank['mismatch_count'] = ''
                        exact_del[f'{cat}_{s}'].append(blank)
                    for b,c in enumerate(exact_del[f'{cat}_{s}']):
                        tmp[f'{cat}_{s}_start-{b}'] = c['start'] 
                        tmp[f'{cat}_{s}_end-{b}'] = c['end'] 
                        tmp[f'{cat}_{s}_stem-{b}'] = c['stem']
                        tmp[f'{cat}_{s}_loop-{b}'] = c['loop']
                        tmp[f'{cat}_{s}_energy-{b}'] = c['energy']
                        tmp[f'{cat}_{s}_sequence-{b}'] = c['sequence']
                        if s[:4] == 'near':
                            tmp[f'{cat}_{s}_mismatch_count{b}'] = c['mismatch_count']
            for cat in [
                'long_range_before_5prime',
                'long_range_after_5prime'
            ]:
                for s in ['hairpins', 'near_hairpins']:
                    while len(exact_del[f'{cat}_{s}']) < 5:
                        blank = {
                            'start': '',
                            'end': '',
                            'stem': '',
                            'loop': '',
                            'energy': '',
                            'sequence': '',
                        }
                        if s[:4] == 'near':
                            blank['mismatch_count'] = ''
                        exact_del[f'{cat}_{s}'].append(blank)
                    for b,c in enumerate(exact_del[f'{cat}_{s}']):
                        tmp[f'{cat}_{s}_start-{b}'] = c['start'] 
                        tmp[f'{cat}_{s}_end-{b}'] = c['end'] 
                        tmp[f'{cat}_{s}_stem-{b}'] = c['stem']
                        tmp[f'{cat}_{s}_loop-{b}'] = c['loop']
                        tmp[f'{cat}_{s}_energy-{b}'] = c['energy']
                        tmp[f'{cat}_{s}_sequence-{b}'] = c['sequence']
                        if s[:4] == 'near':
                            tmp[f'{cat}_{s}_mismatch_count{b}'] = c['mismatch_count']


            # Left side for dynamic programming
            tmp['dp_Left_Sequence'] = dpLeftSequence
            tmp['dp_Left_Energy'] = dpLeftEnergy
            lim5 = 0
            for k,v in dpLeftPathDict.items():
                tmp[f'{k}_LEFT'] = v
                if k[-2:] == 'd2':
                    lim5 += 1
                    if lim5 >= 5: break
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
                if k[-2:] == 'd2':
                    lim5 += 1
                    if lim5 >= 5: break
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

    # Process the original json
    error_lines = 0
    with open('data/top_deletions_5perc_group.json') as f:
        combined_wells = []
        total_wells = 0
        gene_in_prev_well = '2078_5_B274-1'
        broken_datapoints = 0

        # passed_stopping_pt = False
        
        for line in f:
            try:
                tmp = ast.literal_eval(line)
            except:
                error_lines += 1
                continue
            for seq,seq_dels in tmp[0].items():
                gene_id = seq.rsplit('_', 1)[0]
                
                # if gene_id == '2012_1_B265-1':
                #     passed_stopping_pt = True
                # elif not passed_stopping_pt: continue
                
                # process deletion groups across all wells for a sequence
                if gene_id != gene_in_prev_well and gene_in_prev_well:
                    frag_seq = full_seq_id_2_frag_pass[gene_in_prev_well]['frag_sequences'].copy()
                    frag_element = full_seq_id_2_frag_pass[gene_in_prev_well]['element_count'].copy()
                    ordered_frag_seq = []
                    for i in range(len(frag_seq)):
                        max_overlap = ''
                        max_j = 0
                        for j in range(len(frag_seq)):
                            if i == j: continue
                            current_overlap = longest_substring(frag_seq[i][-35:], frag_seq[j][:35])
                            if len(current_overlap) > len(max_overlap):
                                max_overlap = current_overlap
                                max_j = j
                        if len(max_overlap) < 18:
                            max_j = len(frag_seq)
                        ordered_frag_seq.append((max_overlap, max_j, frag_seq[i], frag_element[i]))
                    ordered_frag_seq = sorted(ordered_frag_seq, key=lambda x: x[1])

                    overlaps = [(0,0)]
                    frag_element_correct = []
                    rpos = 0
                    for tup in ordered_frag_seq:
                        overlaps.append([rpos + len(tup[2]) - len(tup[0]), rpos + len(tup[2])])
                        frag_element_correct.append(tup[3])
                        rpos += len(tup[2]) - len(tup[0])

                    process_all_wells(
                        combined_wells,
                        total_wells,
                        gene_in_prev_well,
                        full_seq_id_2_frag_pass[gene_in_prev_well]['sequence'],
                        frag_element_correct,
                        overlaps
                    )
                    # try:
                    #     process_all_wells(
                    #         combined_wells,
                    #         total_wells,
                    #         gene_in_prev_well,
                    #         full_seq_id_2_frag_pass[gene_in_prev_well]['sequence'],
                    #     )
                    # except KeyError:
                    #     print(seq)

                    combined_wells = []
                    total_wells = 0
                    gene_in_prev_well = gene_id

                tmp_dels = reprocess_well(seq_dels, gene_id, full_seq_id_2_frag_pass[gene_id]['concentrations'])
                total_wells += 1
                if total_wells == 1:
                    combined_wells = [tmp_dels]
                else:
                    combined_wells.append(tmp_dels)
                # try:
                #     tmp_dels = reprocess_well(seq_dels, gene_id, full_seq_id_2_frag_pass[gene_id]['concentrations'])
                #     total_wells += 1
                #     if total_wells == 1:
                #         combined_wells = [tmp_dels]
                #     else:
                #         combined_wells.append(tmp_dels)
                # except:
                #     broken_datapoints += 1
                #     print('broke: ', broken_datapoints, f' {seq}')
                #     continue
    

if __name__ == '__main__':
    main()
