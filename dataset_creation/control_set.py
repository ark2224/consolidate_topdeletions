import csv
from pathlib import Path
import gzip
import shutil
import os
import ast
import boto3
import numpy as np
import json
import random

import pandas as pd
import sys
import random

from cred import s3, buck

from data_retrieval import (
    get_frag_lengths2,
    get_number_failed_frags,
    get_conf_data
)

from sequence_manipulation import (
    longest_common_substring,
    max_ad,
    ATGC_energy,
    find_mirrors,
    fix_energy_score,
    find_only_hairpins,
    find_only_repeats,
    longest_substring,
    find_complements
)

from dp_path_finder import find_graph, process_output_path


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def analyzeFakeDeletion(
        strt_end_idx,
        gene_batch_id,
        seq,
        overlaps,
        frag_gc,
        well_ct,
        frag_element
    ):
    strt, end = strt_end_idx
    exact_deletion_data = {
        'del_start': strt,
        'del_end': end,
        'similar_to_past_deletion': 0,
        'read_ct': 0,
        'percentage_of_reads': 0.,
        'log2_total_reads_well': -15
    }

    bases_at_beginning_of_del = ''
    if strt - 50 < 0:
        for _ in range(abs(strt - 50)):
            bases_at_beginning_of_del += 'X'
    bases_at_beginning_of_del += seq[max(strt - 50, 0):strt] + '|' + seq[strt:strt + 50]
    bases_at_very_end = seq[end - 50:end] + '|' + seq[end:min(end + 50, len(seq))]
    if end + 50 > len(seq):
        for _ in range(end + 50 - len(seq)):
            bases_at_very_end += 'X'
    exact_deletion_data['bases_at_beginning_of_deletion'] = bases_at_beginning_of_del
    exact_deletion_data['bases_at_end_of_deletion'] = bases_at_very_end

    five_prime_end = bases_at_beginning_of_del
    three_prime_end = bases_at_very_end

    matching_over_deletion = ''
    for i, c in enumerate(five_prime_end):
        matching_over_deletion += c if i < len(three_prime_end) and c == three_prime_end[i] else '_'
    exact_deletion_data['matching_over_deletion'] = matching_over_deletion[20:81]

    # ===== Right/Left Energy and Offset =====
    # Energy is currently AT-GC Score
    right_energy, right_offset, current_score = 0, 0, 0
    not_found_matching = True
    for c in matching_over_deletion[51:]:
        if c == 'A' or c == 'T':
            current_score += 2
            not_found_matching = False
        elif c == 'G' or c == 'C':
            current_score += 3
            not_found_matching = False
        else:
            current_score -= 1
            if not_found_matching:
                right_offset += 1
        right_energy = max(right_energy, current_score)
    exact_deletion_data['right_energy'] = right_energy
    exact_deletion_data['right_offset'] = right_offset

    left_offset, left_energy, current_score = 0, 0, 0
    not_found_matching = True
    if right_offset == 0:
        j = 51
        while j < len(matching_over_deletion) and matching_over_deletion[j] != '_':
            left_offset -= 1
            if matching_over_deletion[j] == 'A' or matching_over_deletion[j] == 'T':
                current_score += 2
            elif matching_over_deletion[j] == 'G' or matching_over_deletion[j] == 'C':
                current_score += 3
            j += 1
        not_found_matching = False
    for c in matching_over_deletion[49:0:-1]:
        if c == 'A' or c == 'T':
            current_score += 2
            not_found_matching = False
        elif c == 'G' or c == 'C':
            current_score += 3
            not_found_matching = False
        else:
            current_score -= 1
            if not_found_matching:
                left_offset += 1
        left_energy = max(left_energy, current_score)
    exact_deletion_data['left_energy'] = left_energy
    exact_deletion_data['left_offset'] = left_offset

 # ========================================
    bases_at_beginning_of_del = ''
    if strt - 50 < 0:
        for _ in range(abs(strt - 50)):
            bases_at_beginning_of_del += 'X'
    bases_at_beginning_of_del = seq[max(strt - 50, 0):strt + 50]
    bases_at_very_end = seq[end - 50:min(end + 50, len(seq))]
    lR, rR = 0, 0
    idx = 50
    while idx < len(bases_at_very_end) and idx < len(bases_at_beginning_of_del) and bases_at_beginning_of_del[idx] == bases_at_very_end[idx]:
        idx += 1
        lR += 1
        rR += 1

    lL, rL = 0, 0
    idx = 49
    while idx >= 0  and idx < len(bases_at_very_end) and idx < len(bases_at_beginning_of_del) and bases_at_beginning_of_del[idx] == bases_at_very_end[idx]:
        idx -= 1
        lL += 1
        rL += 1

    # Add scoring rules here ========================================
    left = list((matching_over_deletion[:50] + matching_over_deletion[51:51+lR])[::-1])
    right = list((matching_over_deletion[50-rL:50] + matching_over_deletion[51:]))
    for i in range(1,len(left)-1):
        if left[i] != '_' and left[i-1] == '_' and left[i+1] == '_':
            left[i] = '_'
    for i in range(1,len(right)-1):
        if right[i] != '_' and right[i-1] == '_' and right[i+1] == '_':
            right[i] = '_'

    side_names = ['left', 'right']
    scores = [[], []]
    GC_scores = [[],[]]
    for i, side in enumerate([left, right]):
        if side[0].isalpha():
            if side[1].isalpha():
                scores[i].append(0)#single_match
                scores[i].append(0)#tail
            else:
                scores[i].append(1)#single_match
                side[0] = '_'
        else:
            scores[i].append(0)#single_match
        mismatch_score = 0
        match_score = 0
        gc = 0
        for j, c in enumerate(side):
            if c.isalpha():
                if mismatch_score:
                    scores[i].append(mismatch_score)
                    mismatch_score = 0
                match_score += 1
                if c == 'G' or c == 'C':
                    gc += 1
            else:
                if match_score:
                    scores[i].append(match_score)
                    match_score = 0
                    if j >= 30:
                        break
                    GC_scores[i].append(gc)
                    gc = 0
                if j >= 30:
                    break
                mismatch_score += 1
        while len(scores[i]) < 10:
            scores[i].append(0)
        while len(GC_scores[i]) < 5:
            GC_scores[i].append(0)

        for j in range(5):
            exact_deletion_data[f'{side_names[i]}_matching_score_mismatch_{j+1}'] = scores[i][2*j]
            exact_deletion_data[f'{side_names[i]}_matching_score_match_{j+1}'] = scores[i][2*j+1]
        for j in range(5):
            exact_deletion_data[f'{side_names[i]}_match_group_GC_score_content_{j+1}'] = GC_scores[i][j]

    # ===============================================================

    # ===== Max Repeat Anywhere Around Deletion Coordinates =====
    matching_repeat = longest_substring(bases_at_beginning_of_del, bases_at_very_end)
    dist_5prime = bases_at_beginning_of_del.find(matching_repeat) + len(matching_repeat) - 50 - lR
    dist_3prime = 50 - rL - (three_prime_end.find(matching_repeat) + len(matching_repeat))
    exact_deletion_data['longest_repeat'] = matching_repeat
    exact_deletion_data['longest_repeat_dist_5prime'] = dist_5prime
    exact_deletion_data['longest_repeat_dist_3prime'] = dist_3prime
    # ===========================================================

    # ===== Max Beginning-After Repeat =====
    after_del = bases_at_very_end[50-rL:]
    beginning_del = bases_at_beginning_of_del[50+lR:]
    longest_repeat_beginning_del_after_del = longest_substring(beginning_del, after_del)
    repeat_dist_beg_del = beginning_del.find(longest_repeat_beginning_del_after_del) + len(longest_repeat_beginning_del_after_del)#can be zero-indexed
    repeat_dist_after_del = -(1 + after_del.find(longest_repeat_beginning_del_after_del) + len(longest_repeat_beginning_del_after_del))
    exact_deletion_data['longest_repeat_beginning_del_after_del'] = longest_repeat_beginning_del_after_del
    exact_deletion_data['dist_beg_del-longest_repeat_beginning_del_after_del'] = repeat_dist_beg_del
    exact_deletion_data['dist_after_del-longest_repeat_beginning_del_after_del'] = repeat_dist_after_del
    # ======================================

    # ===== Max Before-End Repeat =====
    before_del = bases_at_beginning_of_del[:50+lR]
    end_del = bases_at_very_end[:50-rL]
    longest_repeat_before_del_end_del = longest_substring(before_del, end_del)
    repeat_dist_before_del = 50 + lR - (len(longest_repeat_before_del_end_del) + before_del.find(longest_repeat_before_del_end_del))
    repeat_dist_end_del = 50 - rL - (1 + len(longest_repeat_before_del_end_del) + end_del.find(longest_repeat_before_del_end_del))#can be zero-indexed
    exact_deletion_data['longest_repeat_before_del_end_del'] = longest_repeat_before_del_end_del
    exact_deletion_data['dist_before_del-longest_repeat_before_del_end_del'] = repeat_dist_before_del
    exact_deletion_data['dist_end_del-longest_repeat_before_del_end_del'] = repeat_dist_end_del
    # ==================================

    # ===== Longest Repeat Energy Score =====
    max_repeat = ''
    if len(matching_repeat) > len(longest_repeat_beginning_del_after_del):
        max_repeat = matching_repeat
    else:
        max_repeat = longest_repeat_beginning_del_after_del

    if len(max_repeat) < len(longest_repeat_before_del_end_del):
        max_repeat = longest_repeat_before_del_end_del

    ATGC_score = 0
    for c in max_repeat:
        if c == 'G' or c == 'C':
            ATGC_score += 3
        else:
            ATGC_score += 2
    exact_deletion_data['longest_repeat_score'] = ATGC_score
    # =======================================

    # ===== Max A-D Ligation Repeat =====
    exact_deletion_data['longest_a-d_ligation_repeat'] = max_ad(seq, strt, end)
    # ===================================

    # 3.
    dist_5prime = overlaps[0][0]
    dist_3prime = overlaps[-1][1]
    prev_begin = overlaps[0][0]
    # prev_end = overlaps[0][1]
    exact_deletion_data['within_one_fragment'] = ''
    exact_deletion_data['fragment_gc_content'] = ''
    exact_deletion_data['fragment_is_element'] = ''
    idx = 0
    for (overlap_begin, overlap_end) in overlaps[1:]:
        if prev_begin <= strt <= overlap_end:
            if overlap_begin <= strt:
                dist_5prime = overlap_begin - strt#negative
            else:
                dist_5prime = strt - prev_begin
        if prev_begin <= end <= overlap_end:
            exact_deletion_data['fragment_gc_content'] = frag_gc[idx]
            exact_deletion_data['fragment_is_element'] = frag_element[idx]
            if overlap_begin <= end:
                dist_3prime = end - overlap_end
            else:
                dist_3prime = overlap_end - end#negative
            if strt - abs(dist_5prime) == prev_begin:
                exact_deletion_data['within_one_fragment'] = 'True'
            else:
                exact_deletion_data['within_one_fragment'] = 'False'
        prev_begin = overlap_begin
        idx += 1
    exact_deletion_data['distance_to_fragment_overhang_5prime'] = dist_5prime
    exact_deletion_data['distance_to_fragment_overhang_3prime'] = dist_3prime
    exact_deletion_data['gc_content_frag_start_LR'] = ''
    if dist_5prime <= 30:
        tmp_seq = seq[strt-dist_5prime:strt]
        exact_deletion_data['gc_content_frag_start_LR'] = tmp_seq.count('G') + tmp_seq.count('C')


    offset = (20,23)
    cutoff = (50,0)
    period = (45,47)
    spikeend = (62,26)
    spikestart = (57,18)
    log_spike_ratio = 2.321928095
    decay = 0.6666666667
    cutoff_adjustment = (-0.6, -1.3)
    names = ('left', 'right')
    for i,side in enumerate([dist_5prime,dist_3prime]):
        if side >= cutoff[i]:
            period_position = ((side + offset[i]) % period[i]) - 0.5*period[i]
            Score = 4*(abs(period_position) / period[i])
            Score -= side*(decay / period[i])
        else:
            Score = cutoff_adjustment[i] - 2*(cutoff[i] - side) / period[i]
        
        if spikestart[i] <= side <= spikeend[i]:
            Score += log_spike_ratio
        exact_deletion_data[f'{names[i]}_end_score'] = Score


    exact_deletion_data['latest_start'] = 'False'
    exact_deletion_data['earliest_end'] = 'False'

    exact_deletion_data['perfect_sequence_across_deletion'] = 'True'


    # ==================================================================================================================
    # Now adding in things from main.py that would go into CONSOLIDATE_TOPDELETIONS
    fiveprime_end = (
        exact_deletion_data['bases_at_beginning_of_deletion'][:50] +
        exact_deletion_data['bases_at_beginning_of_deletion'][51:]
    )
    threeprime_end = (
        exact_deletion_data['bases_at_end_of_deletion'][:50] + 
        exact_deletion_data['bases_at_end_of_deletion'][51:]
    )

    # increase boundary to put repeat on right side of deletion; doing this for re-defining longest repeat on left-hand side of deletion
    lR = 50
    while lR < min(len(fiveprime_end), len(threeprime_end)) and fiveprime_end[lR] == threeprime_end[lR]:
        lR += 1
    rR = lR
    rep, dist5, dist3 = longest_common_substring(fiveprime_end[:lR], threeprime_end[:rR])
    exact_deletion_data['longest_repeat_before_del_end_del'] = rep
    exact_deletion_data['dist_before_del-longest_repeat_before_del_end_del'] = dist5 - lR
    exact_deletion_data['dist_end_del-longest_repeat_before_del_end_del'] = lR - dist3

    lL = 49
    while lR >= 0 and lR < min(len(fiveprime_end), len(threeprime_end)) and fiveprime_end[lL] == threeprime_end[lL]:
        lL -= 1
    rL = lL
    rep, dist5, dist3 = longest_common_substring(fiveprime_end[lL:], threeprime_end[rL:])
    exact_deletion_data['longest_repeat_beginning_del_after_del'] = rep
    exact_deletion_data['dist_beg_del-longest_repeat_beginning_del_after_del'] = dist5
    exact_deletion_data['dist_after_del-longest_repeat_beginning_del_after_del'] = -dist3

    # New fields for repeat energy scores
    if exact_deletion_data['longest_repeat_beginning_del_after_del'] == 5:
        exact_deletion_data['longest_repeat_beginning_del_after_del'] = ''
        # ^for a single line error. Figure out what happened there!
    exact_deletion_data['energy_longest_repeat_beginning_del_after_del'] = ATGC_energy(exact_deletion_data['longest_repeat_beginning_del_after_del'])
    exact_deletion_data['energy_longest_repeat_before_del_end_del'] = ATGC_energy(exact_deletion_data['longest_repeat_before_del_end_del'])




    exact_deletion_data['Log2_exact_deletion'] = 0
    exact_deletion_data['Avg_wells_10+log2(frac/(1-frac))'] = -5
    exact_deletion_data['Avg_wells_10+log2(group_frac/(1-group_frac))'] = -5

    # Secondary Structs:
    ds = exact_deletion_data['del_start']
    de = exact_deletion_data['del_end']
    if exact_deletion_data['distance_to_fragment_overhang_5prime'] > 0:
        fs = ds - exact_deletion_data['distance_to_fragment_overhang_5prime']
    else:
        fs = ds + exact_deletion_data['distance_to_fragment_overhang_5prime']
    if exact_deletion_data['distance_to_fragment_overhang_3prime'] > 0:
        fe = de + exact_deletion_data['distance_to_fragment_overhang_3prime']
    else:
        fe = de - exact_deletion_data['distance_to_fragment_overhang_5prime']




    gene_id, batch_num = gene_batch_id.rsplit('_', 1)
    dumpable_data = {gene_id: {}}
    # Exact-Del Secondary Structures

    
    # lL, lR, rL, rR are all offsets from the deletion boundary marked by minimap
    #   all 4 values are absolute (positive) and represent the range of possible "real" deletion boundaries
    #   between lL to lR and rL to rR should be a repeat that makes the deletion boundary ambiguous
    lR, rR = 0, 0
    idx = 50
    while idx < min(len(threeprime_end), len(fiveprime_end)) and fiveprime_end[idx] == threeprime_end[idx]:
        idx += 1
        lR += 1
        rR += 1

    lL, rL = 0, 0
    idx = 49
    while idx >= 0 and fiveprime_end[idx] == threeprime_end[idx]:
        idx -= 1
        lL += 1
        rL += 1
    
    exact_deletion_data['LL'] = -lL
    exact_deletion_data['LR'] = lR
    exact_deletion_data['RL'] = rL
    exact_deletion_data['RR'] = -rR

    # === Dynamic Programming ===
    seq1_dp_left = seq[ds+lR-50:ds+lR]#inboard
    seq2_dp_left = seq[de+rR-50:de+rR]#inside deletion
    dpLeftPath, dpLeftEnergy, dpLeftSequence = find_graph(seq1_dp_left,seq2_dp_left)
    dpLeftPathDict = process_output_path(dpLeftPath)

    seq2_dp_right = seq[ds-lL:ds-lL+50]#inside deletion
    seq1_dp_right = seq[de-rL:de-rL+50]#inboard
    dpRightPath, dpRightEnergy, dpRightSequence = find_graph(seq1_dp_right,seq2_dp_right)
    dpRightPathDict = process_output_path(dpRightPath)


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
    exact_deletion_data['5_prime_hairpins'] = lL100_hairpins['hairpins']
    exact_deletion_data['5_prime_near_hairpins'] = lL100_hairpins['near_hairpins']
    # 5 prime repeats
    lL100_repeats = find_only_repeats(
        fiveprime_srch_window,
        120 - lL,
        120 - lL,
        1,
        1,
        5
    )
    exact_deletion_data['5_prime_repeats'] = lL100_repeats['repeats']
    exact_deletion_data['5_prime_near_repeats'] = lL100_repeats['near_repeats']
    # 5 prime mirrors
    exact_deletion_data['5_prime_mirrors'] = find_mirrors(
        fiveprime_srch_window,
        120 - lL,
        120 - lL,
        1,
        1,
        5
    )
    # 5 prime complements
    exact_deletion_data['5_prime_complements'] = find_complements(
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
    exact_deletion_data['3_prime_hairpins'] = rl100_hairpins['hairpins']
    exact_deletion_data['3_prime_near_hairpins'] = rl100_hairpins['near_hairpins']
    # 3 prime repeats
    rl100_repeats = find_only_repeats(
        threeprime_srch_window,
        120 - rL,
        120 - rL,
        -1,
        -1,
        5
    )
    exact_deletion_data['3_prime_repeats'] = rl100_repeats['repeats']
    exact_deletion_data['3_prime_near_repeats'] = rl100_repeats['near_repeats']
    # 3 prime mirrors
    exact_deletion_data['3_prime_mirrors'] = find_mirrors(
        fiveprime_srch_window,
        120 - rL,
        120 - rL,
        -1,
        -1,
        5
    )
    # 3 prime complements
    exact_deletion_data['3_prime_complements'] = find_complements(
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
    exact_deletion_data['Overall_hairpins'] = overall_hairpins['hairpins'][:5]
    exact_deletion_data['Overall_near_hairpins'] = overall_hairpins['near_hairpins'][:5]
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

    exact_deletion_data.update({
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
    exact_deletion_data['Overall_repeats'] = overall_reps['repeats']
    exact_deletion_data['Overall_near_repeats'] = overall_reps['near_repeats']
    # overall palindromes
    exact_deletion_data['Overall_mirrors'] = find_mirrors(
        overall_srch_window,
        ds-fs-lL,
        de-fs-rL,
        1,
        -1,
        5
    )
    # overall complements
    exact_deletion_data['Overall_complements'] = find_complements(
        overall_srch_window,
        ds-fs-lL,
        de-fs-rL,
        1,
        -1,
        5
    )
    
    # # checking for secondary structures in fragment ends if deletion is over 170bp away from ends
    exact_deletion_data['fragment_start_hairpins'] = blank_hairpin.copy()
    exact_deletion_data['fragment_start_near_hairpins'] = blank_near_hairpin.copy()
    exact_deletion_data['fragment_start_repeats'] = blank_repeat.copy()
    exact_deletion_data['fragment_start_near_repeats'] = blank_near_repeat.copy()
    exact_deletion_data['fragment_start_mirrors'] = blank_repeat.copy()
    exact_deletion_data['fragment_start_complements'] = blank_complement.copy()
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
        exact_deletion_data['fragment_start_hairpins'] = fragment_start_hairpins['hairpins']
        exact_deletion_data['fragment_start_near_hairpins'] = fragment_start_hairpins['near_hairpins']

        fragment_start_reps = find_only_repeats(
            fragment_start_window,
            ds-lL-fs,
            ds-lL-fs,
            1,
            1,
            5
        )
        exact_deletion_data['fragment_start_repeats'] = fragment_start_reps['repeats']
        exact_deletion_data['fragment_start_near_repeats'] = fragment_start_reps['near_repeats']

        exact_deletion_data['fragment_start_mirrors'] = find_mirrors(
            fragment_start_window,
            ds-lL-fs,
            ds-lL-fs,
            1,
            1,
            5
        )

        exact_deletion_data['fragment_start_complements'] = find_complements(
            fragment_start_window,
            ds-lL-fs,
            ds-lL-fs,
            1,
            1,
            5
        )
            
    exact_deletion_data['fragment_end_hairpins'] = blank_hairpin.copy()
    exact_deletion_data['fragment_end_near_hairpins'] = blank_near_hairpin.copy()
    exact_deletion_data['fragment_end_repeats'] = blank_repeat.copy()
    exact_deletion_data['fragment_end_near_repeats'] = blank_near_repeat.copy()
    exact_deletion_data['fragment_end_mirrors'] = blank_repeat.copy()
    exact_deletion_data['fragment_end_complements'] = blank_complement.copy()
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
        exact_deletion_data['fragment_end_hairpins'] = fragment_end_hairpins['hairpins']
        exact_deletion_data['fragment_end_near_hairpins'] = fragment_end_hairpins['near_hairpins']

        fragment_end_reps = find_only_repeats(
            fragment_end_window,
            100+de-rL-fe,
            100+de-rL-fe,
            -1,
            -1,
            5
        )
        exact_deletion_data['fragment_end_repeats'] = fragment_end_reps['repeats']
        exact_deletion_data['fragment_end_near_repeats'] = fragment_end_reps['near_repeats']

        exact_deletion_data['fragment_end_mirrors'] = find_mirrors(
            fragment_end_window,
            100+de-rL-fe,
            100+de-rL-fe,
            -1,
            -1,
            5
        )
        exact_deletion_data['fragment_end_complements'] = find_complements(
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
    ) = fix_energy_score(exact_deletion_data['matching_over_deletion'])
    
    # fixing longest_repeat offset

    long_rep = longest_substring(fiveprime_end, threeprime_end)
    exact_deletion_data['longest_repeat'] = long_rep
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

    tmp = {
        # Sequence meta-data:
        'batch_num': batch_num,
        'total_wells': well_ct,
        'Log2 Total Reads in Well': exact_deletion_data['log2_total_reads_well'],

        # Exact-Deletion Data
        'del_start': exact_deletion_data['del_start'],
        'del_end': exact_deletion_data['del_end'],
        'del_length': exact_deletion_data['del_end'] - exact_deletion_data['del_start'],
        'fragment_is_element': exact_deletion_data['fragment_is_element'],
        'longest_repeat': long_rep,
        'longest_repeat_dist_5prime': fiveprime_offset,
        'longest_repeat_dist_3prime': threeprime_offset,


        'Avg_wells_10+log2(frac/(1-frac))': exact_deletion_data['Avg_wells_10+log2(frac/(1-frac))'],
        'Avg_wells_10+log2(group_frac/(1-group_frac))': exact_deletion_data['Avg_wells_10+log2(group_frac/(1-group_frac))'],
        'Left_End_Score': exact_deletion_data['left_end_score'],
        'Right_End_Score': exact_deletion_data['right_end_score'],
        'LL': exact_deletion_data['LL'],
        'LR': exact_deletion_data['LR'],
        'RL': exact_deletion_data['RL'],
        'RR': exact_deletion_data['RR'],

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
        'read_ct': exact_deletion_data['read_ct'],
        'percentage_of_reads': exact_deletion_data['percentage_of_reads'],
        'matching_over_deletion': exact_deletion_data['matching_over_deletion'],
        'left_Single_Match': exact_deletion_data['left_matching_score_mismatch_1'],
        'left_Tail': exact_deletion_data['left_matching_score_match_1'],
        'left-Match_1': exact_deletion_data['left_matching_score_mismatch_2'],
        'left-GC_content_1': exact_deletion_data['left_match_group_GC_score_content_1'],
        'left-Mismatch_1': exact_deletion_data['left_matching_score_match_2'],
        'left-Match_2': exact_deletion_data['left_matching_score_mismatch_3'],
        'left-GC_content_2': exact_deletion_data['left_match_group_GC_score_content_2'],
        'left-Mismatch_2': exact_deletion_data['left_matching_score_match_3'],
        'left-Match_3': exact_deletion_data['left_matching_score_mismatch_4'],
        'left-GC_content_3': exact_deletion_data['left_match_group_GC_score_content_3'],
        'left-Mismatch_3': exact_deletion_data['left_matching_score_match_4'],
        'left-Match_4': exact_deletion_data['left_matching_score_mismatch_5'],
        'left-GC_content_4': exact_deletion_data['left_match_group_GC_score_content_4'],
        'left-Mismatch_4': exact_deletion_data['left_matching_score_match_5'],

        'right_Single_Match': exact_deletion_data['right_matching_score_mismatch_1'],
        'right_Tail': exact_deletion_data['right_matching_score_match_1'],
        'right-Match_1': exact_deletion_data['right_matching_score_mismatch_2'],
        'right-GC_content_1': exact_deletion_data['right_match_group_GC_score_content_1'],
        'right-Mismatch_1': exact_deletion_data['right_matching_score_match_2'],
        'right-Match_2': exact_deletion_data['right_matching_score_mismatch_3'],
        'right-GC_content_2': exact_deletion_data['right_match_group_GC_score_content_2'],
        'right-Mismatch_2': exact_deletion_data['right_matching_score_match_3'],
        'right-Match_3': exact_deletion_data['right_matching_score_mismatch_4'],
        'right-GC_content_3': exact_deletion_data['right_match_group_GC_score_content_3'],
        'right-Mismatch_3': exact_deletion_data['right_matching_score_match_4'],
        'right-Match_4': exact_deletion_data['right_matching_score_mismatch_5'],
        'right-GC_content_4': exact_deletion_data['right_match_group_GC_score_content_4'],
        'right-Mismatch_4': exact_deletion_data['right_matching_score_match_5'],

        'right_energy': right_energy,
        'right_offset': right_offset,
        'left_energy': left_energy,
        'left_offset': left_offset,
        'gc_content_frag_start_LR': exact_deletion_data['gc_content_frag_start_LR'],

        'longest_repeat_beginning_del_after_del': longest_repeat_beginning_del_after_del,
        'dist_beg_del-longest_repeat_beginning_del_after_del': repeat_dist_beg_del,
        'dist_after_del-longest_repeat_beginning_del_after_del': repeat_dist_after_del,

        'longest_repeat_before_del_end_del': longest_repeat_before_del_end_del,
        'dist_before_del-longest_repeat_before_del_end_del': repeat_dist_before_del,
        'dist_end_del-longest_repeat_before_del_end_del': repeat_dist_end_del,
        
        'longest_repeat_score': exact_deletion_data['longest_repeat_score'],
        'longest_a-d_ligation_repeat': exact_deletion_data['longest_a-d_ligation_repeat'],
        'within_one_fragment': exact_deletion_data['within_one_fragment'],
        'distance_to_fragment_overhang_5prime': exact_deletion_data['distance_to_fragment_overhang_5prime'],
        'distance_to_fragment_overhang_3prime': exact_deletion_data['distance_to_fragment_overhang_3prime'],
        'fragment_gc_content': exact_deletion_data['fragment_gc_content'],

        'perfect_sequence_across_deletion': str(exact_deletion_data['perfect_sequence_across_deletion']),
        'energy_longest_repeat_beginning_del_after_del': exact_deletion_data['energy_longest_repeat_beginning_del_after_del'],
        'energy_longest_repeat_before_del_end_del': exact_deletion_data['energy_longest_repeat_before_del_end_del'],
        'Log2_exact_deletion_frequency': exact_deletion_data['Log2_exact_deletion'],
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
            'palindromes',
            'complements'
        ]:
            while len(exact_deletion_data[f'{cat}_{s}']) < 5:
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
                exact_deletion_data[f'{cat}_{s}'].append(blank)
            for b,c in enumerate(exact_deletion_data[f'{cat}_{s}']):
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
        'long_range_5_3_prime'
    ]:
        for s in ['hairpins', 'near_hairpins']:
            while len(exact_deletion_data[f'{cat}_{s}']) < 2:
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
                exact_deletion_data[f'{cat}_{s}'].append(blank)
            for b,c in enumerate(exact_deletion_data[f'{cat}_{s}']):
                tmp[f'{cat}_{s}_start-{b}'] = c['start'] 
                tmp[f'{cat}_{s}_end-{b}'] = c['end'] 
                tmp[f'{cat}_{s}_stem-{b}'] = c['stem']
                tmp[f'{cat}_{s}_loop-{b}'] = c['loop']
                tmp[f'{cat}_{s}_energy-{b}'] = c['energy']
                tmp[f'{cat}_{s}_sequence-{i}'] = c['sequence']
                if s[:4] == 'near':
                    tmp[f'{cat}_{s}_mismatch_count{b}'] = c['mismatch_count']
    for cat in [
        'long_range_before_5prime',
        'long_range_after_5prime'
    ]:
        for s in ['hairpins', 'near_hairpins']:
            while len(exact_deletion_data[f'{cat}_{s}']) < 5:
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
                exact_deletion_data[f'{cat}_{s}'].append(blank)
            for b,c in enumerate(exact_deletion_data[f'{cat}_{s}']):
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
        if k[-2:] == 'D2':
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
        if k[-2:] == 'D2':
            lim5 += 1
            if lim5 == 5: break
    if lim5 < 5:
        for _ in range(lim5, 5):
            tmp[f'length_{_}'] = 0
            tmp[f'energy_{_}'] = 0
            tmp[f'sequence_{_}'] = ''
            tmp[f'gap{_}_d1'] = 0
            tmp[f'gap{_}_d2'] = 0

    dumpable_data[gene_id][str(exact_deletion_data['del_start']) + '-' + str(exact_deletion_data['del_end'])] = tmp
    with open('data/consolidated_FAKE_top_deletions.json', 'a') as td:
        json.dump(dumpable_data, td, cls=NpEncoder)#save whole dict to save key sequence names
        td.write(',\n')


def processFakeDeletion(geneSet: dict, frag_data: dict, population, frequency):
    for gene, well_ct in geneSet.items():
        seq = frag_data[gene]['sequence']

        frag_seq = frag_data[gene]['frag_sequences'].copy()
        frag_pass = frag_data[gene]['passing_frags'].copy()
        frag_element = frag_data[gene]['element_count'].copy()
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
            ordered_frag_seq.append((max_overlap, max_j, frag_seq[i], frag_pass[i], frag_element[i]))
        ordered_frag_seq = sorted(ordered_frag_seq, key=lambda x: x[1])

        overlaps = [(0,0)]
        frag_pass_correct = []
        frag_element_correct = []
        rpos = 0
        frag_gc_content = []
        for tup in ordered_frag_seq:
            overlaps.append([rpos + len(tup[2]) - len(tup[0]), rpos + len(tup[2])])
            frag_pass_correct.append(tup[3])
            frag_element_correct.append(tup[4])
            rpos += len(tup[2]) - len(tup[0])
            gc_content = tup[2].count('G') + tup[2].count('C')
            frag_gc_content.append(gc_content / len(tup[2]))

        # pick a random start + end point based on distribution of deletion
        # lengths (use ~6000/10000 noDelGene sequences for training and find 5 random
        # deletions per sequence to get equal amounts of real and fake deletions).
        
        # For validation, there is 15210 real deletions, so use the remaining
        # 4000 no-del sequences and pick 4 deletions per sequence
        # current_top_deletions = {
        #     gene: {}
        # }
        for i in range(5):
            del_len = len(seq)
            while len(seq) - 30 - del_len < 0:
                del_len = np.random.choice(population, p=frequency)
            strt = random.randint(0, len(seq) - 30 - del_len)
            strt_end = (strt, strt+del_len)
            analyzeFakeDeletion(
                strt_end,
                gene,
                seq,
                overlaps,
                frag_gc_content,
                well_ct,
                frag_element_correct
            )
        # with open('data/control_deletions.json', 'a') as td:
        #     json.dump(current_top_deletions, td)#save whole dict to save key sequence names
        #     td.write(',\n')

    return


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
    # Getting Confirmation Data
    st_final_file_locations = Path('data/batch_st_final_filename_alphaprod.csv')
    geneid_2_conf = {}
    conf_batches = set()
    missing_st_final = 0
    conf_batch_cnt = 0
    with st_final_file_locations.open('r') as rd:
        maxInt = sys.maxsize
        while True:
            # decrease the maxInt value by factor 10 
            # as long as the OverflowError occurs.
            try:
                csv.field_size_limit(maxInt)
                break
            except OverflowError:
                maxInt = int(maxInt/10)
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

    # pick random sequences and some not-so-random
    random.seed(101)
    # find all sequences that have deletions in them already
    # SAVE SOME FOR VALIDATION FUCK.
    # Get 80% control set from unseen sequences, 20% from sequences with deletions but picking different end points (6% entirely different bounds, 7% with real deletion starts but false ends, 7% with real deletion ends but false starts)
    # Once collected sequences and indices, run dataset analyses on them to make control data points.
    delGenes = {}# gene_id: coordinates
    validation_set = {}
    del_len_population = {}
    total_deletions = 0
    with open('data/consolidated_top_deletions.json') as f:
        for line in f:
            linedict = ast.literal_eval(line)
            if not len(linedict[0]):
                print('here')
                continue
            for seq,seq_dels in linedict[0].items():
                for coordinates, datafields in seq_dels.items():
                    if str(seq[:4]) in {'2013', '1627', '1835', '1833', '1723', '2073', '1759', '2117', '1668', '1625', '1804'}:
                        if seq in validation_set:
                            validation_set[seq].add(coordinates)
                        else:
                            validation_set[seq] = set(coordinates)
                    else:
                        total_deletions += 1
                        if seq not in delGenes:
                            delGenes[seq] = set(coordinates)
                        else:
                            delGenes[seq].add(coordinates)
                        dl = datafields['del_length']
                        if dl in del_len_population:
                            del_len_population[dl] += 1
                        else:
                            del_len_population[dl] = 1
    

    # found all seq with deletions...
    noDelGenes = {}
    for gene_id, conf_data in geneid_2_conf.items():
        if gene_id not in full_seq_id_2_frag_pass:
            continue#need both conf and frag data to make a control set
        if gene_id not in delGenes:
            noDelGenes[gene_id] = conf_data['total_wells']
    print('With Deletions: ', len(delGenes))#1190
    print('Without Deletions: ', len(noDelGenes))#9710
    print('Validation - With Deletions: ', len(validation_set))#815

    # total_deletions = 0
    # for k,v in delGenes.items():
    #     total_deletions += len(v)
    print(total_deletions)#33633

    frequency = []
    population = []
    percent_left = 1.0
    for k,v in del_len_population.items():
        population.append(k)
        if percent_left - float(v / total_deletions) < 0.:
            frequency.append(percent_left - float(v / total_deletions))
            break
        else:
            frequency.append(float(v / total_deletions))
        percent_left -= float(v / total_deletions)
    
    if len(frequency) != len(population):
        return -1

    totalValDel = 0
    for k,v in validation_set.items():
        totalValDel += len(v)
    print(totalValDel)

    processFakeDeletion(
        noDelGenes,
        full_seq_id_2_frag_pass,
        population,
        frequency
    )

    return


if __name__ == '__main__':
    main()