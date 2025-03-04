import csv
import json
import ast
from pathlib import Path
import sys

from data_retrieval import *


def main():
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
        # csv.field_size_limit(sys.maxsize)
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


    # Convert into csv
    error_lines = 0
    with open('data/consolidated_top_deletions_finished2.27.25.json') as f:
        with open('data/consolidated_top_deletions_csv.csv', 'a') as outf:

            # Convert into csv
            error_lines = 0
            completed_genes = set()
            repeats = 0
            writer = csv.writer(outf)
            tmprow = [
                'Sequence ID',
                'Coordinates',
                'Passing Dial Out',

                'Batch Number',
                'Total Wells',
                'Log2 Total Reads in Well',
                
                'Deletion Start',
                'Deletion End',
                'Deletion Length',
                'Longest Repeat',
                'Longest Repeat - Distance from 5\' Deletion-End',
                'Longest Repeat - Distance from 3\' Deletion-End',

                'Average Well 10+log2(frac / (1-frac))',
                'Average Well 10+log2(group_frac / (1-group_frac))',
                'Left End Score',
                'Right End Score',
                'LL',
                'LR',
                'RL',
                'RR',
                '10bp Before Deletion Start',
                '10bp After Deletion Start',
                '10bp Before Deletion End',
                '10bp After Deletion End',


                # 'Deletion Appears in X Wells',
                # 'Similar to Past Deletion',
                'Read Count for Exact Deletion',
                'Percentage of Reads',
                'Matching Over Deletion',
                'Left - Single Match',
                'Left - Tail',
                'Left - First Match',
                'Left - First GC',
                'Left - First Gap',
                'Left - Second Match',
                'Left - Second GC',
                'Left - Second Gap',
                'Left - Third Match',
                'Left - Third GC',
                'Left - Third Gap',
                'Left - Fourth Match',
                'Left - Fourth GC',
                'Left - Fourth Gap',
                # 'Left - Fifth Match',
                # 'Left - Fifth GC',
                # 'Left - Fifth Gap',
                'Right - Single Match',
                'Right - Tail',
                'Right - First Match',
                'Right - First GC',
                'Right - First Gap',
                'Right - Second Match',
                'Right - Second GC',
                'Right - Second Gap',
                'Right - Third Match',
                'Right - Third GC',
                'Right - Third Gap',
                'Right - Fourth Match',
                'Right - Fourth GC',
                'Right - Fourth Gap',
                # 'Right - Fifth Match',
                # 'Right - Fifth GC',
                # 'Right - Fifth Gap',

                'Right Energy',
                'Right Offset',
                'Left Energy',
                'Left Offset',
                'GC Content between fragment overhang - LR',



                # 'Longest Repeat',
                # 'Longest Repeat - Distance from 5\' Deletion-End',
                # 'Longest Repeat - Distance from 3\' Deletion-End',
                'Longest Repeat Beginning Del-after Del',
                'Distance after Deletion Start',
                'Distance after Deletion End',
                'Longest Repeat Before Del-End Del',
                'Distance Before Deletion Start',
                'Distance Before Deletion End',
                'Longest Repeat AT-GC Score',
                'Longest A-B Ligation Repeat',
                'Within One Fragment',
                'Distance to Fragment Beginning',
                'Distance to Fragment End',
                'Fragment GC Content',

                'Perfect Sequence Across Deletion',
                'Energy Longest Repeat Beginning Del-after Del',
                'Energy Longest Repeat Before Del-End Del',
                'Log2 Exact Deletion Frequency',
                # 'Log2 Deletion Group Frequency',

                # # Exact deletion secondary structures
                '5_prime_hairpin_start-1',
                '5_prime_hairpin_end-1',
                '5_prime_hairpin_stem-1',
                '5_prime_hairpin_loop-1',
                '5_prime_hairpin_energy-1',

                '5_prime_hairpin_start-2',
                '5_prime_hairpin_end-2',
                '5_prime_hairpin_stem-2',
                '5_prime_hairpin_loop-2',
                '5_prime_hairpin_energy-2',


                '5_prime_hairpin_start-3',
                '5_prime_hairpin_end-3',
                '5_prime_hairpin_stem-3',
                '5_prime_hairpin_loop-3',
                '5_prime_hairpin_energy-3',

                '5_prime_hairpin_start-4',
                '5_prime_hairpin_end-4',
                '5_prime_hairpin_stem-4',
                '5_prime_hairpin_loop-4',
                '5_prime_hairpin_energy-4',

                '5_prime_hairpin_start-5',
                '5_prime_hairpin_end-5',
                '5_prime_hairpin_stem-5',
                '5_prime_hairpin_loop-5',
                '5_prime_hairpin_energy-5',


                '5_prime_near_hairpin_start-1',
                '5_prime_near_hairpin_end-1',
                '5_prime_near_hairpin_stem-1',
                '5_prime_near_hairpin_loop-1',
                '5_prime_near_hairpin_energy-1',
                '5_prime_near_hairpin_mismatch_count-1',

                '5_prime_near_hairpin_start-2',
                '5_prime_near_hairpin_end-2',
                '5_prime_near_hairpin_stem-2',
                '5_prime_near_hairpin_loop-2',
                '5_prime_near_hairpin_energy-2',
                '5_prime_near_hairpin_mismatch_count-2',

                '5_prime_near_hairpin_start-3',
                '5_prime_near_hairpin_end-3',
                '5_prime_near_hairpin_stem-3',
                '5_prime_near_hairpin_loop-3',
                '5_prime_near_hairpin_energy-3',
                '5_prime_near_hairpin_mismatch_count-3',

                '5_prime_near_hairpin_start-4',
                '5_prime_near_hairpin_end-4',
                '5_prime_near_hairpin_stem-4',
                '5_prime_near_hairpin_loop-4',
                '5_prime_near_hairpin_energy-4',
                '5_prime_near_hairpin_mismatch_count-4',

                '5_prime_near_hairpin_start-5',
                '5_primenear_hairpin_end-5',
                '5_prime_near_hairpin_stem-5',
                '5_prime_near_hairpin_loop-5',
                '5_prime_near_hairpin_energy-5',
                '5_prime_near_hairpin_mismatch_count-5',


                '3_prime_hairpin_start-1',
                '3_prime_hairpin_end-1',
                '3_prime_hairpin_stem-1',
                '3_prime_hairpin_loop-1',
                '3_prime_hairpin_energy-1',

                '3_prime_hairpin_start-2',
                '3_prime_hairpin_end-2',
                '3_prime_hairpin_stem-2',
                '3_prime_hairpin_loop-2',
                '3_prime_hairpin_energy-2',


                '3_prime_hairpin_start-3',
                '3_prime_hairpin_end-3',
                '3_prime_hairpin_stem-3',
                '3_prime_hairpin_loop-3',
                '3_prime_hairpin_energy-3',

                '3_prime_hairpin_start-4',
                '3_prime_hairpin_end-4',
                '3_prime_hairpin_stem-4',
                '3_prime_hairpin_loop-4',
                '3_prime_hairpin_energy-4',

                '3_prime_hairpin_start-5',
                '3_prime_hairpin_end-5',
                '3_prime_hairpin_stem-5',
                '3_prime_hairpin_loop-5',
                '3_prime_hairpin_energy-5',


                '3_prime_near_hairpin_start-1',
                '3_prime_near_hairpin_end-1',
                '3_prime_near_hairpin_stem-1',
                '3_prime_near_hairpin_loop-1',
                '3_prime_near_hairpin_energy-1',
                '3_prime_near_hairpin_mismatch_count-1',

                '3_prime_near_hairpin_start-2',
                '3_prime_near_hairpin_end-2',
                '3_prime_near_hairpin_stem-2',
                '3_prime_near_hairpin_loop-2',
                '3_prime_near_hairpin_energy-2',
                '3_prime_near_hairpin_mismatch_count-2',

                '3_prime_near_hairpin_start-3',
                '3_prime_near_hairpin_end-3',
                '3_prime_near_hairpin_stem-3',
                '3_prime_near_hairpin_loop-3',
                '3_prime_near_hairpin_energy-3',
                '3_prime_near_hairpin_mismatch_count-3',

                '3_prime_near_hairpin_start-4',
                '3_prime_near_hairpin_end-4',
                '3_prime_near_hairpin_stem-4',
                '3_prime_near_hairpin_loop-4',
                '3_prime_near_hairpin_energy-4',
                '3_prime_near_hairpin_mismatch_count-4',

                '3_prime_near_hairpin_start-5',
                '3_primenear_hairpin_end-5',
                '3_prime_near_hairpin_stem-5',
                '3_prime_near_hairpin_loop-5',
                '3_prime_near_hairpin_energy-5',
                '3_prime_near_hairpin_mismatch_count-5',

                # 'IHFL_hairpin_start-Optimal',
                # 'IHFL_hairpin_end-Optimal',
                # 'IHFL_hairpin_stem-Optimal',
                # 'IHFL_hairpin_loop-Optimal',
                # 'IHFL_hairpin_energy-Optimal',

                # 'IHFL_hairpin_start-Sub_Optimal',
                # 'IHFL_hairpin_end-Sub_Optimal',
                # 'IHFL_hairpin_stem-Sub_Optimal',
                # 'IHFL_hairpin_loop-Sub_Optimal',
                # 'IHFL_hairpin_energy-Sub_Optimal',

                # 'IHFL_near_hairpin_start-Optimal',
                # 'IHFL_near_hairpin_end-Optimal',
                # 'IHFL_near_hairpin_stem-Optimal',
                # 'IHFL_near_hairpin_loop-Optimal',
                # 'IHFL_near_hairpin_energy-Optimal',
                # 'IHFL_near_hairpin_mismatch_count-Optimal',

                # 'IHFL_near_hairpin_start-Sub_Optimal',
                # 'IHFL_near_hairpin_end-Sub_Optimal',
                # 'IHFL_near_hairpin_stem-Sub_Optimal',
                # 'IHFL_near_hairpin_loop-Sub_Optimal',
                # 'IHFL_near_hairpin_energy-Sub_Optimal',
                # 'IHFL_near_hairpin_mismatch_count-Sub_Optimal',


                # 'PDFL_hairpin_start-Optimal',
                # 'PDFL_hairpin_end-Optimal',
                # 'PDFL_hairpin_stem-Optimal',
                # 'PDFL_hairpin_loop-Optimal',
                # 'PDFL_hairpin_energy-Optimal',

                # 'PDFL_hairpin_start-Sub_Optimal',
                # 'PDFL_hairpin_end-Sub_Optimal',
                # 'PDFL_hairpin_stem-Sub_Optimal',
                # 'PDFL_hairpin_loop-Sub_Optimal',
                # 'PDFL_hairpin_energy-Sub_Optimal',

                # 'PDFL_near_hairpin_start-Optimal',
                # 'PDFL_near_hairpin_end-Optimal',
                # 'PDFL_near_hairpin_stem-Optimal',
                # 'PDFL_near_hairpin_loop-Optimal',
                # 'PDFL_near_hairpin_energy-Optimal',
                # 'PDFL_near_hairpin_mismatch_count-Optimal',

                # 'PDFL_near_hairpin_start-Sub_Optimal',
                # 'PDFL_near_hairpin_end-Sub_Optimal',
                # 'PDFL_near_hairpin_stem-Sub_Optimal',
                # 'PDFL_near_hairpin_loop-Sub_Optimal',
                # 'PDFL_near_hairpin_energy-Sub_Optimal',
                # 'PDFL_near_hairpin_mismatch_count-Sub_Optimal',


                # 'IHFR_hairpin_start-Optimal',
                # 'IHFR_hairpin_end-Optimal',
                # 'IHFR_hairpin_stem-Optimal',
                # 'IHFR_hairpin_loop-Optimal',
                # 'IHFR_hairpin_energy-Optimal',

                # 'IHFR_hairpin_start-Sub_Optimal',
                # 'IHFR_hairpin_end-Sub_Optimal',
                # 'IHFR_hairpin_stem-Sub_Optimal',
                # 'IHFR_hairpin_loop-Sub_Optimal',
                # 'IHFR_hairpin_energy-Sub_Optimal',

                # 'IHFR_near_hairpin_start-Optimal',
                # 'IHFR_near_hairpin_end-Optimal',
                # 'IHFR_near_hairpin_stem-Optimal',
                # 'IHFR_near_hairpin_loop-Optimal',
                # 'IHFR_near_hairpin_energy-Optimal',
                # 'IHFR_near_hairpin_mismatch_count-Optimal',

                # 'IHFR_near_hairpin_start-Sub_Optimal',
                # 'IHFR_near_hairpin_end-Sub_Optimal',
                # 'IHFR_near_hairpin_stem-Sub_Optimal',
                # 'IHFR_near_hairpin_loop-Sub_Optimal',
                # 'IHFR_near_hairpin_energy-Sub_Optimal',
                # 'IHFR_near_hairpin_mismatch_count-Sub_Optimal',


                # 'PDFR_hairpin_start-Optimal',
                # 'PDFR_hairpin_end-Optimal',
                # 'PDFR_hairpin_stem-Optimal',
                # 'PDFR_hairpin_loop-Optimal',
                # 'PDFR_hairpin_energy-Optimal',

                # 'PDFR_hairpin_start-Sub_Optimal',
                # 'PDFR_hairpin_end-Sub_Optimal',
                # 'PDFR_hairpin_stem-Sub_Optimal',
                # 'PDFR_hairpin_loop-Sub_Optimal',
                # 'PDFR_hairpin_energy-Sub_Optimal',

                # 'PDFR_near_hairpin_start-Optimal',
                # 'PDFR_near_hairpin_end-Optimal',
                # 'PDFR_near_hairpin_stem-Optimal',
                # 'PDFR_near_hairpin_loop-Optimal',
                # 'PDFR_near_hairpin_energy-Optimal',
                # 'PDFR_near_hairpin_mismatch_count-Optimal',

                # 'PDFR_near_hairpin_start-Sub_Optimal',
                # 'PDFR_near_hairpin_end-Sub_Optimal',
                # 'PDFR_near_hairpin_stem-Sub_Optimal',
                # 'PDFR_near_hairpin_loop-Sub_Optimal',
                # 'PDFR_near_hairpin_energy-Sub_Optimal',
                # 'PDFR_near_hairpin_mismatch_count-Sub_Optimal',


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

                'Overall_reverse_repeat_start-1',
                'Overall_reverse_repeat_end-1',
                'Overall_reverse_repeat_length-1',
                'Overall_reverse_repeat_energy-1',

                'Overall_reverse_repeat_start-2',
                'Overall_reverse_repeat_end-2',
                'Overall_reverse_repeat_length-2',
                'Overall_reverse_repeat_energy-2',

                'Overall_reverse_repeat_start-3',
                'Overall_reverse_repeat_end-3',
                'Overall_reverse_repeat_length-3',
                'Overall_reverse_repeat_energy-3',

                'Overall_reverse_repeat_start-4',
                'Overall_reverse_repeat_end-4',
                'Overall_reverse_repeat_length-4',
                'Overall_reverse_repeat_energy-4',

                'Overall_reverse_repeat_start-5',
                'Overall_reverse_repeat_end-5',
                'Overall_reverse_repeat_length-5',
                'Overall_reverse_repeat_energy-5',




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

                # REMOVED FRAGMENT CONCENTRATION TEMPORARILY
                # 'Total Fragment Count',
                # 'Log 2 Fragment Concentrations',
            ]
            tmprow.extend(['DP Left Sequence', 'DP Left Energy'])
            for _ in range(5):
                tmprow.extend([
                    f'Match{_} Len - Left',
                    f'Match{_} Energy - Left',
                    f'Match{_} Sequence - Left',
                    f'Gap{_} D1 - Left',
                    f'Gap{_} D2 - Left',
                ])
            tmprow.extend(['DP Right Sequence', 'DP Right Energy'])
            for _ in range(5):
                tmprow.extend([
                    f'Match{_} Len - Right',
                    f'Match{_} Energy - Right',
                    f'Match{_} Sequence - Right',
                    f'Gap{_} D1 - Right',
                    f'Gap{_} D2 - Right',
                ])
            writer.writerow(tmprow)
            for line in f:
                tmp = ast.literal_eval(line)
                if not len(tmp[0]):
                    print('here')
                    continue
                for seq,seq_dels in tmp[0].items():
                    for coordinates,data_fields in seq_dels.items():
                        seq_coor = str(seq) + str(coordinates)
                        if seq_coor in completed_genes or data_fields['batch_num'] not in conf_batches:
                            continue
                        # try:
                        completed_genes.add(seq_coor)
                        if data_fields['within_one_fragment'] == 'False':
                            continue

                        conf_idx = seq + '_' + data_fields['batch_num']
                        if conf_idx not in geneid_2_conf:
                            passing_dial_out = 0
                        else:
                            passing_dial_out = int(geneid_2_conf[conf_idx]['wells_passing'] > 0)


                        tmprow = [seq,
                                  coordinates,
                                  passing_dial_out,
                                  ]
                        for k,v in data_fields.items():
                            if k == 'within_one_fragment':
                                tmprow.append(1 if v else 0)
                            elif k == 'perfect_sequence_across_deletion':
                                tmprow.append(1 if v else 0)
                            elif k == 'Log2_fragment_concentrations':
                                frag_conc = ast.literal_eval(str(v))
                                tmprow.append(len(frag_conc))
                                for conc in frag_conc:
                                    tmprow.append(conc)
                            elif k == 'Log2_group_deletion_frequency':
                                # tmprow.append(max(v, -15))
                                continue
                            elif k == 'del_in_x_wells':
                                continue
                            else:
                                tmprow.append(v)
                        writer.writerow(tmprow)
    print(error_lines)
                        # tmprow = [
                        #     seq, 
                        #     coordinates,
                        #     data_fields['batch_num'],
                        #     data_fields['total_wells'],
                        #     data_fields['del_start'],
                        #     data_fields['del_end'],
                        #     data_fields['del_length'],
                        #     data_fields['Avg_wells_10+log2(frac/(1-frac))'],
                        #     data_fields['Avg_wells_10+log2(group_frac/(1-group_frac))'],
                        #     data_fields['Left_End_Score'],
                        #     data_fields['Right_End_Score'],
                        #     data_fields['LL'],
                        #     data_fields['LR'],
                        #     data_fields['RL'],
                        #     data_fields['RR'],
                        #     data_fields['10bp_before_deletion'],
                        #     data_fields['10bp_after_deletion'],
                        #     data_fields['read_ct'],
                        #     data_fields['percentage_of_reads'],
                        #     data_fields['matching_over_deletion'],
                        #     data_fields['left_matching_score_mismatch_1'],
                        #     data_fields['left_matching_score_match_1'],
                        #     data_fields['left_match_group_GC_score_content_1'],
                        #     data_fields['left_matching_score_mismatch_2'],
                        #     data_fields['left_matching_score_match_2'],
                        #     data_fields['left_match_group_GC_score_content_2'],
                        #     data_fields['left_matching_score_mismatch_3'],
                        #     data_fields['left_matching_score_match_3'],
                        #     data_fields['left_match_group_GC_score_content_3'],
                        #     data_fields['left_matching_score_mismatch_4'],
                        #     data_fields['left_matching_score_match_4'],
                        #     data_fields['left_match_group_GC_score_content_4'],
                        #     # data_fields['left_matching_score_mismatch_5'],
                        #     # data_fields['left_matching_score_match_5'],
                        #     # data_fields['left_match_group_GC_score_content_5'],

                        #     data_fields['right_matching_score_mismatch_1'],
                        #     data_fields['right_matching_score_match_1'],
                        #     data_fields['right_match_group_GC_score_content_1'],
                        #     data_fields['right_matching_score_mismatch_2'],
                        #     data_fields['right_matching_score_match_2'],
                        #     data_fields['right_match_group_GC_score_content_2'],
                        #     data_fields['right_matching_score_mismatch_3'],
                        #     data_fields['right_matching_score_match_3'],
                        #     data_fields['right_match_group_GC_score_content_3'],
                        #     data_fields['right_matching_score_mismatch_4'],
                        #     data_fields['right_matching_score_match_4'],
                        #     data_fields['right_match_group_GC_score_content_4'],
                        #     # data_fields['right_matching_score_mismatch_5'],
                        #     # data_fields['right_matching_score_match_5'],
                        #     # data_fields['right_match_group_GC_score_content_5'],

                        #     data_fields['right_energy'],
                        #     data_fields['right_offset'],
                        #     data_fields['left_energy'],
                        #     data_fields['left_offset'],

                        #     data_fields['longest_repeat'],
                        #     data_fields['longest_repeat_dist_5prime'],
                        #     data_fields['longest_repeat_dist_3prime'],
                        #     data_fields['longest_repeat_beginning_del_after_del'],
                        #     data_fields['dist_beg_del-longest_repeat_beginning_del_after_del'],
                        #     data_fields['dist_after_del-longest_repeat_beginning_del_after_del'],
                        #     data_fields['longest_repeat_before_del_end_del'],
                        #     data_fields['dist_before_del-longest_repeat_before_del_end_del'],
                        #     data_fields['dist_end_del-longest_repeat_before_del_end_del'],
                        #     data_fields['longest_repeat_score'],
                        #     data_fields['longest_a-d_ligation_repeat'],
                        #     data_fields['within_one_fragment'],
                        #     data_fields['distance_to_fragment_overhang_5prime'],
                        #     data_fields['distance_to_fragment_overhang_3prime'],
                        #     # data_fields['gc_content_frag_start_LR'],
                        #     data_fields['left_end_score'],
                        #     data_fields['right_end_score'],
                        #     # data_fields['latest_start'],
                        #     # data_fields['earliest_end'],
                        #     data_fields['perfect_sequence_across_deletion'],
                        #     data_fields['energy_longest_repeat_beginning_del_after_del'],
                        #     data_fields['energy_longest_repeat_before_del_end_del'],
                        #     data_fields['Log2_exact_deletion_frequency'],
                        #     data_fields['Log2_deletion_group_frequency'],
                        #     data_fields['5_prime_hairpins_start_0'],
                        #     data_fields['5_prime_hairpins_end_0'],
                        #     data_fields['5_prime_hairpins_stem_0'],
                        #     data_fields['5_prime_hairpins_loop_0'],
                        #     data_fields['5_prime_hairpins_energy_0'],
                        #     data_fields['5_prime_hairpins_start_1'],
                        #     data_fields['5_prime_hairpins_end_1'],
                        #     data_fields['5_prime_hairpins_stem_1'],
                        #     data_fields['5_prime_hairpins_loop_1'],
                        #     data_fields['5_prime_hairpins_energy_1'],
                        #     data_fields['5_prime_hairpins_start_2'],
                        #     data_fields['5_prime_hairpins_end_2'],
                        #     data_fields['5_prime_hairpins_stem_2'],
                        #     data_fields['5_prime_hairpins_loop_2'],
                        #     data_fields['5_prime_hairpins_energy_2'],
                        #     data_fields['5_prime_hairpins_start_3'],
                        #     data_fields['5_prime_hairpins_end_3'],
                        #     data_fields['5_prime_hairpins_stem_3'],
                        #     data_fields['5_prime_hairpins_loop_3'],
                        #     data_fields['5_prime_hairpins_energy_3'],
                        #     data_fields['5_prime_hairpins_start_4'],
                        #     data_fields['5_prime_hairpins_end_4'],
                        #     data_fields['5_prime_hairpins_stem_4'],
                        #     data_fields['5_prime_hairpins_loop_4'],
                        #     data_fields['5_prime_hairpins_energy_4'],
                        #     data_fields['5_prime_near_hairpins_start_0'],
                        #     data_fields['5_prime_near_hairpins_end_0'],
                        #     data_fields['5_prime_near_hairpins_stem_0'],
                        #     data_fields['5_prime_near_hairpins_loop_0'],
                        #     data_fields['5_prime_near_hairpins_energy_0'],
                        #     data_fields['5_prime_near_hairpins_mismatch_count_0'],
                        #     data_fields['5_prime_near_hairpins_start_1'],
                        #     data_fields['5_prime_near_hairpins_end_1'],
                        #     data_fields['5_prime_near_hairpins_stem_1'],
                        #     data_fields['5_prime_near_hairpins_loop_1'],
                        #     data_fields['5_prime_near_hairpins_energy_1'],
                        #     data_fields['5_prime_near_hairpins_mismatch_count_1'],
                        #     data_fields['5_prime_near_hairpins_start_2'],
                        #     data_fields['5_prime_near_hairpins_end_2'],
                        #     data_fields['5_prime_near_hairpins_stem_2'],
                        #     data_fields['5_prime_near_hairpins_loop_2'],
                        #     data_fields['5_prime_near_hairpins_energy_2'],
                        #     data_fields['5_prime_near_hairpins_mismatch_count_2'],
                        #     data_fields['5_prime_near_hairpins_start_3'],
                        #     data_fields['5_prime_near_hairpins_end_3'],
                        #     data_fields['5_prime_near_hairpins_stem_3'],
                        #     data_fields['5_prime_near_hairpins_loop_3'],
                        #     data_fields['5_prime_near_hairpins_energy_3'],
                        #     data_fields['5_prime_near_hairpins_mismatch_count_3'],
                        #     data_fields['5_prime_near_hairpins_start_4'],
                        #     data_fields['5_prime_near_hairpins_end_4'],
                        #     data_fields['5_prime_near_hairpins_stem_4'],
                        #     data_fields['5_prime_near_hairpins_loop_4'],
                        #     data_fields['5_prime_near_hairpins_energy_4'],
                        #     data_fields['5_prime_near_hairpins_mismatch_count_4'],
                            
                        #     data_fields['3_prime_hairpins_start_0'],
                        #     data_fields['3_prime_hairpins_end_0'],
                        #     data_fields['3_prime_hairpins_stem_0'],
                        #     data_fields['3_prime_hairpins_loop_0'],
                        #     data_fields['3_prime_hairpins_energy_0'],
                        #     data_fields['3_prime_hairpins_start_1'],
                        #     data_fields['3_prime_hairpins_end_1'],
                        #     data_fields['3_prime_hairpins_stem_1'],
                        #     data_fields['3_prime_hairpins_loop_1'],
                        #     data_fields['3_prime_hairpins_energy_1'],
                        #     data_fields['3_prime_hairpins_start_2'],
                        #     data_fields['3_prime_hairpins_end_2'],
                        #     data_fields['3_prime_hairpins_stem_2'],
                        #     data_fields['3_prime_hairpins_loop_2'],
                        #     data_fields['3_prime_hairpins_energy_2'],
                        #     data_fields['3_prime_hairpins_start_3'],
                        #     data_fields['3_prime_hairpins_end_3'],
                        #     data_fields['3_prime_hairpins_stem_3'],
                        #     data_fields['3_prime_hairpins_loop_3'],
                        #     data_fields['3_prime_hairpins_energy_3'],
                        #     data_fields['3_prime_hairpins_start_4'],
                        #     data_fields['3_prime_hairpins_end_4'],
                        #     data_fields['3_prime_hairpins_stem_4'],
                        #     data_fields['3_prime_hairpins_loop_4'],
                        #     data_fields['3_prime_hairpins_energy_4'],
                        #     data_fields['3_prime_near_hairpins_start_0'],
                        #     data_fields['3_prime_near_hairpins_end_0'],
                        #     data_fields['3_prime_near_hairpins_stem_0'],
                        #     data_fields['3_prime_near_hairpins_loop_0'],
                        #     data_fields['3_prime_near_hairpins_energy_0'],
                        #     data_fields['3_prime_near_hairpins_mismatch_count_0'],
                        #     data_fields['3_prime_near_hairpins_start_1'],
                        #     data_fields['3_prime_near_hairpins_end_1'],
                        #     data_fields['3_prime_near_hairpins_stem_1'],
                        #     data_fields['3_prime_near_hairpins_loop_1'],
                        #     data_fields['3_prime_near_hairpins_energy_1'],
                        #     data_fields['3_prime_near_hairpins_mismatch_count_1'],
                        #     data_fields['3_prime_near_hairpins_start_2'],
                        #     data_fields['3_prime_near_hairpins_end_2'],
                        #     data_fields['3_prime_near_hairpins_stem_2'],
                        #     data_fields['3_prime_near_hairpins_loop_2'],
                        #     data_fields['3_prime_near_hairpins_energy_2'],
                        #     data_fields['3_prime_near_hairpins_mismatch_count_2'],
                        #     data_fields['3_prime_near_hairpins_start_3'],
                        #     data_fields['3_prime_near_hairpins_end_3'],
                        #     data_fields['3_prime_near_hairpins_stem_3'],
                        #     data_fields['3_prime_near_hairpins_loop_3'],
                        #     data_fields['3_prime_near_hairpins_energy_3'],
                        #     data_fields['3_prime_near_hairpins_mismatch_count_3'],
                        #     data_fields['3_prime_near_hairpins_start_4'],
                        #     data_fields['3_prime_near_hairpins_end_4'],
                        #     data_fields['3_prime_near_hairpins_stem_4'],
                        #     data_fields['3_prime_near_hairpins_loop_4'],
                        #     data_fields['3_prime_near_hairpins_energy_4'],
                        #     data_fields['3_prime_near_hairpins_mismatch_count_4'],

                        #     data_fields['Overall_hairpins_start_0'],
                        #     data_fields['Overall_hairpins_end_0'],
                        #     data_fields['Overall_hairpins_stem_0'],
                        #     data_fields['Overall_hairpins_loop_0'],
                        #     data_fields['Overall_hairpins_energy_0'],
                        #     data_fields['Overall_hairpins_start_1'],
                        #     data_fields['Overall_hairpins_end_1'],
                        #     data_fields['Overall_hairpins_stem_1'],
                        #     data_fields['Overall_hairpins_loop_1'],
                        #     data_fields['Overall_hairpins_energy_1'],
                        #     data_fields['Overall_hairpins_start_2'],
                        #     data_fields['Overall_hairpins_end_2'],
                        #     data_fields['Overall_hairpins_stem_2'],
                        #     data_fields['Overall_hairpins_loop_2'],
                        #     data_fields['Overall_hairpins_energy_2'],
                        #     data_fields['Overall_hairpins_start_3'],
                        #     data_fields['Overall_hairpins_end_3'],
                        #     data_fields['Overall_hairpins_stem_3'],
                        #     data_fields['Overall_hairpins_loop_3'],
                        #     data_fields['Overall_hairpins_energy_3'],
                        #     data_fields['Overall_hairpins_start_4'],
                        #     data_fields['Overall_hairpins_end_4'],
                        #     data_fields['Overall_hairpins_stem_4'],
                        #     data_fields['Overall_hairpins_loop_4'],
                        #     data_fields['Overall_hairpins_energy_4'],
                        #     data_fields['Overall_near_hairpins_start_0'],
                        #     data_fields['Overall_near_hairpins_end_0'],
                        #     data_fields['Overall_near_hairpins_stem_0'],
                        #     data_fields['Overall_near_hairpins_loop_0'],
                        #     data_fields['Overall_near_hairpins_energy_0'],
                        #     data_fields['Overall_near_hairpins_mismatch_count_0'],
                        #     data_fields['Overall_near_hairpins_start_1'],
                        #     data_fields['Overall_near_hairpins_end_1'],
                        #     data_fields['Overall_near_hairpins_stem_1'],
                        #     data_fields['Overall_near_hairpins_loop_1'],
                        #     data_fields['Overall_near_hairpins_energy_1'],
                        #     data_fields['Overall_near_hairpins_mismatch_count_1'],
                        #     data_fields['Overall_near_hairpins_start_2'],
                        #     data_fields['Overall_near_hairpins_end_2'],
                        #     data_fields['Overall_near_hairpins_stem_2'],
                        #     data_fields['Overall_near_hairpins_loop_2'],
                        #     data_fields['Overall_near_hairpins_energy_2'],
                        #     data_fields['Overall_near_hairpins_mismatch_count_2'],
                        #     data_fields['Overall_near_hairpins_start_3'],
                        #     data_fields['Overall_near_hairpins_end_3'],
                        #     data_fields['Overall_near_hairpins_stem_3'],
                        #     data_fields['Overall_near_hairpins_loop_3'],
                        #     data_fields['Overall_near_hairpins_energy_3'],
                        #     data_fields['Overall_near_hairpins_mismatch_count_3'],
                        #     data_fields['Overall_near_hairpins_start_4'],
                        #     data_fields['Overall_near_hairpins_end_4'],
                        #     data_fields['Overall_near_hairpins_stem_4'],
                        #     data_fields['Overall_near_hairpins_loop_4'],
                        #     data_fields['Overall_near_hairpins_energy_4'],
                        #     data_fields['Overall_near_hairpins_mismatch_count_4'],

                        #     data_fields['Overall_repeats_start_0'],
                        #     data_fields['Overall_repeats_end_0'],
                        #     data_fields['Overall_repeats_stem_0'],
                        #     data_fields['Overall_repeats_loop_0'],
                        #     data_fields['Overall_repeats_energy_0'],
                        #     data_fields['Overall_repeats_start_1'],
                        #     data_fields['Overall_repeats_end_1'],
                        #     data_fields['Overall_repeats_stem_1'],
                        #     data_fields['Overall_repeats_loop_1'],
                        #     data_fields['Overall_repeats_energy_1'],
                        #     data_fields['Overall_repeats_start_2'],
                        #     data_fields['Overall_repeats_end_2'],
                        #     data_fields['Overall_repeats_stem_2'],
                        #     data_fields['Overall_repeats_loop_2'],
                        #     data_fields['Overall_repeats_energy_2'],
                        #     data_fields['Overall_repeats_start_3'],
                        #     data_fields['Overall_repeats_end_3'],
                        #     data_fields['Overall_repeats_stem_3'],
                        #     data_fields['Overall_repeats_loop_3'],
                        #     data_fields['Overall_repeats_energy_3'],
                        #     data_fields['Overall_repeats_start_4'],
                        #     data_fields['Overall_repeats_end_4'],
                        #     data_fields['Overall_repeats_stem_4'],
                        #     data_fields['Overall_repeats_loop_4'],
                        #     data_fields['Overall_repeats_energy_4'],
                        #     data_fields['Overall_near_repeats_start_0'],
                        #     data_fields['Overall_near_repeats_end_0'],
                        #     data_fields['Overall_near_repeats_stem_0'],
                        #     data_fields['Overall_near_repeats_loop_0'],
                        #     data_fields['Overall_near_repeats_energy_0'],
                        #     data_fields['Overall_near_repeats_mismatch_count_0'],
                        #     data_fields['Overall_near_repeats_start_1'],
                        #     data_fields['Overall_near_repeats_end_1'],
                        #     data_fields['Overall_near_repeats_stem_1'],
                        #     data_fields['Overall_near_repeats_loop_1'],
                        #     data_fields['Overall_near_repeats_energy_1'],
                        #     data_fields['Overall_near_repeats_mismatch_count_1'],
                        #     data_fields['Overall_near_repeats_start_2'],
                        #     data_fields['Overall_near_repeats_end_2'],
                        #     data_fields['Overall_near_repeats_stem_2'],
                        #     data_fields['Overall_near_repeats_loop_2'],
                        #     data_fields['Overall_near_repeats_energy_2'],
                        #     data_fields['Overall_near_repeats_mismatch_count_2'],
                        #     data_fields['Overall_near_repeats_start_3'],
                        #     data_fields['Overall_near_repeats_end_3'],
                        #     data_fields['Overall_near_repeats_stem_3'],
                        #     data_fields['Overall_near_repeats_loop_3'],
                        #     data_fields['Overall_near_repeats_energy_3'],
                        #     data_fields['Overall_near_repeats_mismatch_count_3'],
                        #     data_fields['Overall_near_repeats_start_4'],
                        #     data_fields['Overall_near_repeats_end_4'],
                        #     data_fields['Overall_near_repeats_stem_4'],
                        #     data_fields['Overall_near_repeats_loop_4'],
                        #     data_fields['Overall_near_repeats_energy_4'],
                        #     data_fields['Overall_near_repeats_mismatch_count_4'],
                        #     data_fields['Log2_fragment_concentrations'],
                        # ]






if __name__ == '__main__':
    main()