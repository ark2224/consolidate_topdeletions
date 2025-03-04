import pysam
import csv
from pathlib import Path
import gzip
import shutil
import os
import ast
import boto3

import pandas as pd
import sys

from cred import s3, buck


from data_retrieval import (
    get_frag_lengths2,
    get_number_failed_frags,
    get_conf_data
)
from sequence_manipulation import (
    longest_substring,
    coverage,
    find_perfect_sequences_across_deletion,
    process_deletions,
    plot_cigar,
)

TAF = 'GCGAGTCTTAGCCTGCGACGCT'
TAR = 'AGACGCGCAGACGACGACACAC'

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

    # # Getting Confirmation Data
    # st_final_file_locations = Path('data/batch_st_final_filename_alphaprod.csv')
    # geneid_2_conf = {}
    # conf_batches = set()
    # missing_st_final = 0
    # conf_batch_cnt = 0
    # with st_final_file_locations.open('r') as rd:
    #     csv.field_size_limit(sys.maxsize)
    #     reader = csv.DictReader(rd)
    #     for row in reader:
    #         try:
    #             obj = s3.get_object(Bucket=buck, Key='static/elegen_csv/st_final_csv/' + row['st_final_filename'])
    #             batchname = row['batch_number'].partition('_')[0]
    #             date = row['date'].partition(' ')[0]
    #             conf_batch_cnt += 1
    #         except:
    #             missing_st_final += 1
    #             continue
    #         st_final_content = pd.read_csv(obj['Body'])
    #         geneid_2_conf.update(get_conf_data(st_final_content, batchname, date))
    #         conf_batches.add(batchname)
    #         print(batchname)

    # print(len(geneid_2_conf))
    # print("missing: ", missing_st_final)
    # print('batch count: ', conf_batch_cnt)

    # Collect the ReAmp Alignment Files location in AWS S3
    Sams_dict = {}
    csvfile = Path('data/conf_sams_250122.csv')
    with csvfile.open('r') as rd:
        csv.field_size_limit(sys.maxsize)
        reader = csv.DictReader(rd)
        for row in reader:
            batch = row['batch'].split('_')[0]
            sample_to_gene_dict = ast.literal_eval(row['sample_to_gene_id'])
            try:
                sams = ast.literal_eval(row['sams'])
            except:
                # some of the batches don't have sams
                continue
            for example in sams:
                plate_pos = example['left']
                sam = example['right']
                if plate_pos[-4:] != 'conf':
                    print(plate_pos, sam)
                    return
                gene_id = sample_to_gene_dict[plate_pos][0] + '_' + batch
                if gene_id not in Sams_dict:
                    Sams_dict[gene_id] = [sam]
                else:
                    Sams_dict[gene_id].append(sam)

    gene_count = {}
    Passed_last_run_sequence = False
    for k,v in Sams_dict.items():
        if k in full_seq_id_2_frag_pass:#will only be ~3300
            gene_id = k

            # for whenever you pick up, just look at aws to see how many files you've
            # already done and add name you want to pick up at:
            # if not Passed_last_run_sequence:
            #     if gene_id == '853_10_B20240303':
            #         Passed_last_run_sequence = True
            #     else:
            #         print(gene_id)
            #         continue

            # getting old sam
            for sam_name in v:
                if gene_id not in gene_count:
                    gene_count[gene_id] = 1
                else:
                    gene_count[gene_id] += 1
                sam = sam_name[14:]
                print(sam)

                obj = s3.get_object(Bucket='bfx-prod', Key=sam)
                with gzip.GzipFile(fileobj=obj['Body'], mode='rb') as gz_file:
                    with open('tmpfile.sam', 'wb') as outf:
                        shutil.copyfileobj(gz_file, outf)

                # creating new sam
                seq = full_seq_id_2_frag_pass[gene_id]['sequence']
                os.system('samtools fastq tmpfile.sam > tmp_file.fastq')
                with open('tmp_file.fasta', 'w') as fasta:
                    fasta.write(f'>{k}\n')
                    fasta.write(TAF + seq + TAR)
                os.system(f'minimap2 -ax map-ont -A 2 -B 2 -O 2,40 -E 2,0 --end-bonus 50 tmp_file.fasta tmp_file.fastq > {k}_{gene_count[gene_id]}.sam')

                # # save to aws:
                # try:
                #     s3.upload_file(
                #         Path(f'/content/{k}_{gene_count[gene_id]}.sam'),
                #         'bfx-prod',
                #         f'rnd/new_alignment_algorithm_results/A2_B2_O2,40_E2_endbonus50/{k}_{gene_count[gene_id]}.sam'
                #     )
                # except:
                #     print('failed to upload sam to aws')

                # create cigar plot
                frag_seq = full_seq_id_2_frag_pass[gene_id]['frag_sequences'].copy()
                frag_pass = full_seq_id_2_frag_pass[gene_id]['passing_frags'].copy()
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
                    ordered_frag_seq.append((max_overlap, max_j, frag_seq[i], frag_pass[i]))
                ordered_frag_seq = sorted(ordered_frag_seq, key=lambda x: x[1])

                overlaps = [(0,0)]
                frag_pass_correct = []
                rpos = 0
                frag_gc_content = []
                for tup in ordered_frag_seq:
                    overlaps.append([rpos + len(tup[2]) - len(tup[0]), rpos + len(tup[2])])
                    frag_pass_correct.append(tup[3])
                    rpos += len(tup[2]) - len(tup[0])
                    gc_content = tup[2].count('G') + tup[2].count('C')
                    frag_gc_content.append(gc_content / len(tup[2]))

                # COVERAGE =========================================================
                with pysam.Samfile(f'{k}_{gene_count[gene_id]}.sam') as samf:
                    print(f'Processing {gene_id}_{gene_count[gene_id]}...')
                    d = coverage(TAF + seq + TAR, samf)
                    total_reads = d['mapped'] + d['unmapped'] + d['supp']
                    print('TOTAL_READS', total_reads, '============================================')

                    # Plotting CIGAR Pile-Up
                    boundaries = {"b": [len(TAF + seq + TAR)]}
                    filename = gene_id + f'_cigarPlot{gene_count[gene_id]}.png'
                    plot_cigar(
                        d['depth'],
                        filename,
                        boundaries,
                        flen=overlaps,
                        # fpass=frag_pass
                        fpass=frag_pass_correct
                    )

                with pysam.Samfile(f'{k}_{gene_count[gene_id]}.sam') as samf:
                    ordered_deletions = find_perfect_sequences_across_deletion(
                        TAF + seq + TAR,
                        samf,
                        sorted(d['deletions'], key=lambda x: x[3], reverse=True)
                    )
                    ordered_deletions = sorted(ordered_deletions, key=lambda x: x[3], reverse=True)
                    process_deletions(
                        seq,
                        ordered_deletions,
                        total_reads,
                        gene_id,
                        gene_count[gene_id],
                        overlaps,
                        frag_gc_content
                    )
                    # Saving file
                    # try:
                    #     s3.upload_file(
                    #         Path(f'/content/drive/MyDrive/new_alignment_with_cigar_overlay_A2_B2_O2,40_E2_endbonus50/{filename}'),
                    #         'bfx-prod',
                    #         f'rnd/new_alignment_algorithm_results/A2_B2_O2,40_E2_endbonus50/new_alignment_cigar_plots/{filename}'
                    #     )
                    # except:
                    #     print('failed to upload file to aws')
                    os.remove('tmpfile.sam')
                    os.remove(f'{k}_{gene_count[gene_id]}.sam')
                    os.remove('tmp_file.fasta')
                    os.remove('tmp_file.fastq')

if __name__ == '__main__':
    main()