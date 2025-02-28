import speedy
import secrets
from difflib import SequenceMatcher
import pysam
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import json
from math import log2


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


# def find_secondary_structures(del_segment: str, frag_seq: str) -> dict:
#     try:
#         repeats = find_repeats(del_segment)
#         hairpins = find_hairpins(del_segment)
#         near_repeats = find_near_repeats(del_segment)
#         near_hairpins = find_near_hairpins(del_segment)
#     except:
#         print(del_segment)
#         repeats, hairpins, near_repeats, near_hairpins = [], [], [], []

#     # repeats, hairpins, near_repeats, near_hairpins = [], [], [], []
#     # ToDo: calculate energies and re-order them
#     return {
#         'repeats': repeats,
#         'hairpins': hairpins,
#         'near_repeats': near_repeats,
#         'near_hairpins': near_hairpins
#     }


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




# The rest from NewDataset_12182024_cleaner.ipynb
# New coverage function for recording deletions
def coverage(
        ref_seq,
        samfile,
    ):
    '''
    Inputs:
      - reference sequence string
      - samfile from QC
    Returns:
      - Dict({dd: (dict) of cigar results,
              mapped: (int) of mapped
              unmapped: (int) of unmapped,
              supp: (int) of supp})
              deletions: (list of tuples) all overlapping deletions in reads
    '''
    deletions = []

    # dd := list of dicts, one for each reference nucleotide. Dict records each operation at given position
    mapped, unmapped, supp = 0, 0, 0
    # iterate over querries in samfile:
    for aln in samfile:
        # if NOT an error
        if not aln.flag & (pysam.FUNMAP | pysam.FSECONDARY | pysam.FSUPPLEMENTARY):
            mapped += 1
            rpos = aln.reference_start
            qpos = 0
            # iterate over cigarstring/tuples
            for op, siz in aln.cigartuples:
                # separate by operation, then size of operation
                # SOFTCLIP operation; only advances query
                if op == pysam.CSOFT_CLIP:
                    qpos += siz
                # INSERTION operation; only advances query
                elif op == pysam.CINS:
                    qpos += siz
                # MATCH/MISMATCH/DELETION
                elif op in (pysam.CMATCH, pysam.CDIFF, pysam.CDEL):
                    # iterate over size
                    if op == pysam.CDEL and siz >= 10:# used to be 20
                        strt = rpos
                        end = rpos + siz
                        added = False
                        strt_end_key = str(strt) + '-' + str(end)
                        for i in range(len(deletions)):
                            if (abs(strt - deletions[i][0]) <= 100) and (abs(end - deletions[i][1]) <= 100):
                                if strt_end_key in deletions[i][2]:
                                    deletions[i][2][strt_end_key][0] += 1
                                else:
                                    deletions[i][2][strt_end_key] = [1, aln.query_sequence[max(qpos - 50, 0):qpos] + '|' + aln.query_sequence[qpos:qpos+50]]
                                deletions[i][0] = (strt + deletions[i][0]) // 2
                                deletions[i][1] = (end + deletions[i][1]) // 2
                                deletions[i][3] += 1
                                added = True
                                break
                        if not added:
                            deletions.append([
                                strt,#start of deletion GROUP
                                end,#end of deletion GROUP
                                {strt_end_key: [1,#exact deletion within the deletion group
                                                aln.query_sequence[max(qpos - 50, 0):qpos] + '|' + aln.query_sequence[qpos:qpos+50],
                                                ]
                                },
                                1#deletion GROUP read count
                            ])
                    # MATCH/MISMATCH; advances both query and ref
                    if op in (pysam.CMATCH, pysam.CDIFF):
                        # # MATCH
                        rpos += siz
                        qpos += siz
                        # DELETION; only advances reference
                    elif op == pysam.CDEL:
                        rpos += siz
            # if by right here, rpos < len(ref_seq): collect deletion
            if rpos + 10 < len(ref_seq):
                siz = len(ref_seq) - rpos
                strt = rpos
                end = rpos + siz
                added = False
                strt_end_key = str(strt) + '-' + str(end)
                for i in range(len(deletions)):
                    if (abs(strt - deletions[i][0]) <= 100) and (abs(end - deletions[i][1]) <= 100):
                        if strt_end_key in deletions[i][2]:
                            deletions[i][2][strt_end_key][0] += 1
                        else:
                            deletions[i][2][strt_end_key] = [1, aln.query_sequence[max(qpos - 50, 0):qpos] + '|' + aln.query_sequence[qpos:qpos+50]]
                        deletions[i][0] = (strt + deletions[i][0]) // 2
                        deletions[i][1] = (end + deletions[i][1]) // 2
                        deletions[i][3] += 1
                        added = True
                        break
                if not added:
                    deletions.append([
                        strt,#start of deletion GROUP
                        end,#end of deletion GROUP
                        {strt_end_key: [1,#exact deletion within the deletion group
                                        aln.query_sequence[max(qpos - 50, 0):qpos] + '|',
                                        ]
                        },
                        1#deletion GROUP read count
                    ])
        # if unmapped
        elif aln.flag & pysam.FUNMAP:
            unmapped += 1
        # if other mapped
        elif aln.flag & (pysam.FSECONDARY | pysam.FSUPPLEMENTARY):
            supp += 1
    return {
        "mapped": mapped,
        "unmapped": unmapped,
        "supp": supp,
        "deletions": deletions,
    }


# def plot_cigar(dd,
#                fname,
#                boundaries,
#                flen,
#                fpass):
#     plt.clf()
#     if not dd:
#         plt.text(0.5, 0.5, f"no reads for\n{fname}", ha="center")
#         plt.savefig(fname)
#         return

#     # add minor tick-marks to X-axis
#     fig, ax = plt.subplots()
#     ax.xaxis.set_minor_locator(AutoMinorLocator())
#     ax.tick_params(which="both", width=1)
#     ax.tick_params(which="major", length=7)
#     ax.tick_params(which="minor", length=4, color="grey")

#     ref = fname.split(".")[0]
#     X = np.array(range(min(dd), max(dd) + 1))
#     bottom = np.zeros(len(X))
#     for nom, k in (
#         ("match", "m"),
#         ("del", "d"),
#         ("mismatch", "x"),
#         ("ins", "i"),
#         ("clip", "c"),
#     ):
#         Y = np.array([dd[i][k] for i in X])
#         # bottom is used to stack bar charts, and width=1 eliminates gaps between bars
#         plt.bar(X, Y, bottom=bottom, label=f"{ref} - {nom}", width=1.0, alpha=0.8)
#         bottom += Y
#     for i in boundaries['b']:
#         plt.axvline(i, color="black", linestyle="dashed", alpha=0.2)


#     # plotting fragment bounds
#     idx = 0.0
#     for i, (s,e) in enumerate(flen):
#         for j in range(s,e):
#             plt.axvline(j, color='black', linestyle='dashed', alpha=0.2)
#         if not fpass[i]:
#             plt.axvspan(idx, e, color='red', alpha=0.2)
#         idx = e


#     plt.ylabel("depth")
#     plt.xlabel("position")
#     plt.legend()
#     # plt.savefig(f'/content/drive/MyDrive/new_alignment_with_cigar_overlay_A2_B2_O2,40_E2_endbonus50/{fname}')


def find_perfect_sequences_across_deletion(ref_seq, samfile, deletions):
    # deletions = [
    #       [strt,end,{strtend_key: sample_sequence, total_reads}, total_reads_group],
    #       [strt,end...],
    #       ...,
    #   ]
    for i in range(len(deletions)):
        deletions[i].append(False)
    for aln in samfile:
        early_return = True
        errors = np.zeros((len(ref_seq) + 1))
        # if NOT an error
        if not aln.flag & (pysam.FUNMAP | pysam.FSECONDARY | pysam.FSUPPLEMENTARY):
            rpos = aln.reference_start
            qpos = 0
            # iterate over cigarstring/tuples
            for op, siz in aln.cigartuples:
                # separate by operation, then size of operation
                # SOFTCLIP/INSERTION operations; only advances query
                if op in (pysam.CSOFT_CLIP, pysam.CINS):
                    qpos += siz
                    errors[rpos] = 1
                # MATCH/MISMATCH/DELETION
                elif op in (pysam.CMATCH, pysam.CDIFF):
                # advance both ref+query
                    for i in range(siz):
                        if aln.query_sequence[qpos+i] != ref_seq[rpos+i]:
                            errors[rpos+i] = 1
                    rpos += siz
                    qpos += siz
                # DELETION; only advances reference
                elif op == pysam.CDEL:
                    errors[rpos:rpos+siz] = 1
                    rpos += siz

            for i, dels in enumerate(deletions):
                if dels[4]:
                    continue
                else:
                    early_return = False# all deletions contain a perfect sequence
                if sum(errors[dels[0]:dels[1]]) < 2:
                    deletions[i][4] = True
            if early_return:
                print('============ got an early return ==============')
                break
    for dels in deletions:
        print(dels[4])
        if not dels[4]:
            print(dels[:2])
    return deletions


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


def max_ad(seq, strt_idx, end_idx) -> str:
    q = [(1,1)]
    attempted = set((1,1))
    max_seq = ''
    while len(q) != 0:
        a,b = q.pop(0)
        if strt_idx-a < 0 or end_idx+b >= len(seq):
            continue
        curr_seq = seq[strt_idx - a:strt_idx] + seq[end_idx:end_idx + b]
        if seq[:strt_idx-a].find(curr_seq) != -1 or seq[end_idx+b:].find(curr_seq) != -1:
            newa = (a+1,b)
            newb = (a,b+1)
            if newa not in attempted:
                q.append(newa)
                attempted.add(newa)
            if newb not in attempted:
                q.append(newb)
                attempted.add(newb)
            if len(max_seq) < len(curr_seq):
                max_seq = curr_seq
    return max_seq


def analyze_exact_deletion(seq, strt_end_idx, info, total_reads, perfect, deletion_coords, overlaps) -> dict:
    strt, end = strt_end_idx.split('-')
    # Take out TAF (22bp long) to get correct indices relative to just the reference sequence
    strt = int(strt) - 22
    end = int(end) - 22
    similar_2_past_del = []
    for (s,e) in deletion_coords:
        if abs(strt - s) + abs(end - e) <= 5:
            similar_2_past_del = [s, e]
            break
    deletion_coords.append((strt, end))

    exact_deletion_data = {
        'del_start': strt,
        'del_end': end,
        'similar_to_past_deletion': str(similar_2_past_del) if len(similar_2_past_del) else 0,
        'read_ct': info[0],
        'percentage_of_reads': float(info[0] / total_reads),
        'log2_total_reads_well': log2(total_reads)
    }

    # 1. Bases +/- 50bp of deletion coordinates
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
    # exact_deletion_data['sample_read_sequence'] = info[1]

    # 2. Finding Repeats within +/- 30bp window
    # five_prime_end = bases_at_beginning_of_del[20:81]
    # three_prime_end = bases_at_very_end[20:81]
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
    # if repeat_idx_5prime < 20 and len(matching_repeat) > 20 - repeat_idx_5prime:
    #     dist_5prime = 0
    # elif repeat_idx_5prime < 20:
    #     dist_5prime = 20 - (repeat_idx_5prime + len(matching_repeat))
    # else:
    #     dist_5prime = 20 - repeat_idx_5prime
    dist_3prime = 50 - rL - (three_prime_end.find(matching_repeat) + len(matching_repeat))
    # if repeat_idx_3prime < 20 and len(matching_repeat) > 20 - repeat_idx_3prime:
    #     dist_3prime = 0
    # elif repeat_idx_3prime < 20:
    #     dist_3prime = 20 - (repeat_idx_3prime + len(matching_repeat))
    # else:
    #     dist_3prime = 20 - repeat_idx_3prime
    exact_deletion_data['longest_repeat'] = matching_repeat
    exact_deletion_data['longest_repeat_dist_5prime'] = dist_5prime
    exact_deletion_data['longest_repeat_dist_3prime'] = dist_3prime
    # ===========================================================

    # ===== Max Beginning-After Repeat =====
    # after_del = three_prime_end[21:]
    # beginning_del = five_prime_end[21:]
    # longest_repeat_beginning_del_after_del = longest_substring(beginning_del, after_del)
    # repeat_dist_beg_del = beginning_del.find(longest_repeat_beginning_del_after_del)
    # repeat_dist_after_del = after_del.find(longest_repeat_beginning_del_after_del)
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
    # before_del = five_prime_end[:20]
    # end_del = three_prime_end[:20]
    # longest_repeat_before_del_end_del = longest_substring(before_del, end_del)
    # repeat_dist_before_del = 20 - (len(longest_repeat_before_del_end_del) + before_del.find(longest_repeat_before_del_end_del))
    # repeat_dist_end_del = 20 - (len(longest_repeat_before_del_end_del) + end_del.find(longest_repeat_before_del_end_del))
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
    for (overlap_begin, overlap_end) in overlaps[1:]:
        if prev_begin <= strt <= overlap_end:
            if overlap_begin <= strt:
                dist_5prime = overlap_begin - strt#negative
            else:
                dist_5prime = strt - prev_begin

        if prev_begin <= end <= overlap_end:
            if overlap_begin <= end:
                dist_3prime = end - overlap_end
            else:
                dist_3prime = overlap_end - end#negative

            if strt - abs(dist_5prime) == prev_begin:
                exact_deletion_data['within_one_fragment'] = 'True'
            else:
                exact_deletion_data['within_one_fragment'] = 'False'


        # if overlap_begin <= strt:
        #     if strt <= overlap_end:
        #         dist_5prime = overlap_begin - strt#negative
        #     else:
        #         dist_5prime = strt - overlap_begin
        # if end <= overlap_end:
        #     if overlap_begin <= end:
        #         dist_3prime = end - overlap_end
        #     else:
        #         dist_3prime = overlap_end - end#negative
        # else:
        #     exact_deletion_data['within_one_fragment'] = 'True' if dist_5prime == prev_prev_begin else 'False'
        #     break
        prev_begin = overlap_begin
        # prev_end = overlap_end
    exact_deletion_data['distance_to_fragment_overhang_5prime'] = dist_5prime
    # exact_deletion_data['distance_to_fragment_overhang_5prime_2'] = dist_5prime[1]
    # exact_deletion_data['distance_to_fragment_overhang_3prime_1'] = dist_3prime[0]
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

    exact_deletion_data['perfect_sequence_across_deletion'] = str(perfect)

    return exact_deletion_data, deletion_coords


def process_deletions(
        seq,
        deletions,
        total_reads,
        gene_id,
        well_ct,
        overlaps,
        min_ratio=0.05
    ) -> None:
    # for top 5 deletions ==========================================
    #   1. find +/- 30bp from end,
    #   2. (max) # of bases matching within a 40 window,
    #   3. distance of each deletion from endpoint to the eB overhang (which can be negative)
    current_top_deletions = {f'{gene_id}_{well_ct}': {}}
    idx = 0
    # iterate through all deletion groups
    while idx < len(deletions):
        group_sum = deletions[idx][3]
        print('Del group count: ',  group_sum)
        print('Precentage of reads: ', group_sum / total_reads)

        # Checking that the group of deletions affects a minimum % of total reads
        if group_sum / total_reads < min_ratio:
            break

        # Sorting the exact deletions within a group by how many reads have this exact deletion
        top_deletions = {k: v for k,v in sorted(deletions[idx][2].items(), key=lambda item: item[1][0], reverse=True)}
        perfect = deletions[idx][4]
        inner_idx = 0
        latest_strt = 0
        earliest_end = len(seq)
        past_deletion_coords = []
        # Process top 10 exact deletions for a deletion group
        for strt_end_idx, info in top_deletions.items():
            # only recording deletion coordinates that correspond to at least 7 reads; used to do 4 but was too much information
            ct = info[0]
            # if ct < 7:
            #     continue
            if float(ct / total_reads) < 0.002:
                print('CUT!!!')
                break
            if inner_idx > 10:
                break
            current_top_deletions[f'{gene_id}_{well_ct}'][strt_end_idx], past_deletion_coords = analyze_exact_deletion(
                seq,
                strt_end_idx,
                info,
                total_reads,
                perfect,
                past_deletion_coords,
                overlaps
            )
            current_top_deletions[f'{gene_id}_{well_ct}'][strt_end_idx]['group_percent'] = float(group_sum / total_reads)

            inner_idx += 1

        # # Taking the below out because no longer interested in latest-start earliest-end for a deletion group
        # for strt_end_idx, info in top_deletions.items():
        #     strt, end = strt_end_idx.split('-')
        #     strt = int(strt) - 22
        #     end = int(end) - 22
        #     latest_strt = max(latest_strt, strt)
        #     earliest_end = min(earliest_end, end)
        # for strt_end_idx, info in top_deletions.items():
        #     strt, end = strt_end_idx.split('-')
        #     # Take out TAF (22bp long) to get correct indices relative to just the reference sequence
        #     strt = int(strt) - 22
        #     end = int(end) - 22
        #     if strt == latest_strt:
        #         if strt_end_idx not in current_top_deletions[f'{gene_id}_{well_ct}']:
        #             current_top_deletions[f'{gene_id}_{well_ct}'][strt_end_idx], _ = analyze_exact_deletion(seq, strt_end_idx, info, total_reads, perfect, [], overlaps)
        #         current_top_deletions[f'{gene_id}_{well_ct}'][strt_end_idx]['latest_start'] = 'True'
        #     if end == earliest_end:
        #         if strt_end_idx not in current_top_deletions[f'{gene_id}_{well_ct}']:
        #             current_top_deletions[f'{gene_id}_{well_ct}'][strt_end_idx], _ = analyze_exact_deletion(seq, strt_end_idx, info, total_reads, perfect, [], overlaps)
        #         current_top_deletions[f'{gene_id}_{well_ct}'][strt_end_idx]['earliest_end'] = 'True'
        idx += 1

    # save deletions
    with open('data/top_deletions_5perc_group.json', 'a') as td:
        json.dump(current_top_deletions, td)#save whole dict to save key sequence names
        td.write(',\n')