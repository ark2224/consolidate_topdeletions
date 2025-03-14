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
import os
import ast
import seaborn as sns


def find_energy_splits(all_energies: dict) -> dict:
    energy_lims = {}
    preset_lims = [0.06, 0.12, 0.25, 0.5]
    for k,v in all_energies.items():
        #highest energy first
        all_energies[k] = sorted(v, reverse=True)
        energy_indices = []
        for lim in preset_lims:
            energy_indices.append(
                all_energies[k][int(len(v) * lim)]
            )
        print(k, energy_lims[k])
    return energy_lims


def create_control_heatmaps(energy_lims: dict, file_name: str) -> dict:
    completed_genes = set()
    heatmap_matrices = {}
    stem5p_heatmaps = {}
    stem3p_heatmaps = {}
    unique_structures = {}
    structure_type = [
        'hairpins',
        'near_hairpins',
        'repeats',
        'near_repeats',
        'mirrors',
        'complements'
    ]
    regions = [
        '5_prime',
        '3_prime',
        'Overall',
        'fragment_start',
        'fragment_end'
    ]
    bound = 200
    for t in structure_type:
        for reg in regions:
            unique_structures[f'{reg}_{t}'] = set()
            heatmap_matrices[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            stem5p_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            stem3p_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))

    with open(file_name) as f:
        for line in f:
            tmp = ast.literal_eval(line)
            for seq,seq_dels in tmp[0].items():
                for coordinates,data_fields in seq_dels.items():
                    seq_coor = str(seq) + str(coordinates)
                    if seq_coor in completed_genes:
                        continue
                    completed_genes.add(seq_coor)

                    subcat = ['0', '1', '2', '3', '4']
                    for reg in regions:
                        for t in structure_type:
                            for sc in subcat:
                                if not data_fields[f'{reg}_{t}_energy-{sc}']:
                                    continue
                                hp = [
                                    data_fields[f'{reg}_{t}_energy-{sc}'],
                                    data_fields[f'{reg}_{t}_start-{sc}'],
                                    data_fields[f'{reg}_{t}_end-{sc}'],
                                ]
                                if reg == '3_prime':
                                    #            del
                                    # 5' + <------|------> - 3'
                                    hp[1] *= -1
                                    hp[2] *= -1
                                hp[1] += bound
                                hp[2] += bound

                                # Avoid repeats:
                                tup = tuple(hp)
                                if tup in unique_structures[f'{reg}_{t}']:
                                    continue
                                unique_structures[f'{reg}_{t}'].add(tup)

                                en = tup[0]
                                idx = 0
                                for e in energy_lims[f'{reg}_{t}']:
                                    if en < e:
                                        idx += 1
                                heatmap_matrices[f'{reg}_{t}'][idx, tup[1], tup[2]] += 1

                                if tup[1]+tup[3] >= 2*bound: continue
                                stem5p_heatmaps[f'{reg}_{t}'][idx, tup[1], tup[1]+tup[3]] += 1
                                # stem3p is 5'-end vs 3'-end relative to 3'-del (positive going INTO deletion)
                                if tup[2]+tup[3] >= 2*bound: continue
                                stem3p_heatmaps[f'{reg}_{t}'][idx, tup[2], tup[2]+tup[3]] += 1
    return heatmap_matrices, stem5p_heatmaps, stem3p_heatmaps


def fiveP_threeP_heatmaps(energy_lims: dict, controls: dict, file_name: str):
    completed_genes = set()
    heatmap_matrices = {}
    unique_structures = {}
    # Local Structures 3.4.25
    structure_type = [
        'hairpins',
        'near_hairpins',
        'repeats',
        'near_repeats',
        'mirrors',
        'complements'
    ]
    regions = [
        '5_prime',
        '3_prime',
    ]
    bound = 200

    for t in structure_type:
        for reg in regions:
            unique_structures[f'{reg}_{t}'] = set()
            heatmap_matrices[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            # for _ in range(0, 2*bound, 10):
            #     heatmap_matrices[f'{reg}_{t}'][:,bound,_] += 5
            #     heatmap_matrices[f'{reg}_{t}'][:,_,bound] += 5
            # heatmap_matrices[f'{reg}_{t}'][:,bound,bound] -= 5

    # with open('data/consolidated_top_deletions_finished_3.7.25.json') as f:
    # with open('data/consolidated_FAKE_top_deletions.json') as f:
    with open(file_name) as f:
        for line in f:
            tmp = ast.literal_eval(line)
            for seq,seq_dels in tmp[0].items():
                for coordinates,data_fields in seq_dels.items():
                    seq_coor = str(seq) + str(coordinates)
                    if seq_coor in completed_genes:
                        continue
                    completed_genes.add(seq_coor)

                    subcat = ['0', '1', '2', '3', '4']
                    for reg in regions:
                        for t in structure_type:
                            for sc in subcat:
                                if not data_fields[f'{reg}_{t}_energy-{sc}']:
                                    continue
                                hp = [
                                    data_fields[f'{reg}_{t}_energy-{sc}'],
                                    data_fields[f'{reg}_{t}_start-{sc}'],
                                    data_fields[f'{reg}_{t}_end-{sc}'],
                                ]
                                if reg == '3_prime':
                                    #            del
                                    # 5' + <------|------> - 3'
                                    hp[1] *= -1
                                    hp[2] *= -1
                                hp[1] += bound
                                hp[2] += bound

                                # Avoid repeats:
                                tup = tuple(hp)
                                if tup in unique_structures[f'{reg}_{t}']:
                                    continue
                                unique_structures[f'{reg}_{t}'].add(tup)

                                en = tup[0]
                                idx = 0
                                for e in energy_lims[f'{reg}_{t}']:
                                    if en < e:
                                        idx += 1
                                heatmap_matrices[f'{reg}_{t}'][idx, tup[1], tup[2]] += 1

    # Normalize 4x4 case window with a 6x6 control window
    normalized_heatmaps = {}
    for t in structure_type:
        for reg in regions:
            normalized_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
    for k,v in heatmap_matrices.items():
        for m in range(4):#channel
            for n in range(2, 2*bound-3):#row
                for o in range(2, 2*bound-3):#column
                    norm = 9*(1 + sum(heatmap_matrices[k][m, n:n+2, o:o+2]) / (1 + sum(controls[k][m,n:n+2,o:o+2])))
                    normalized_heatmaps[k][m, n:n+2, o:o+2] = heatmap_matrices[k][m, n:n+2, o:o+2] / norm
        for _ in range(0, 2*bound, 10):
            normalized_heatmaps[k][:,bound,_] += 5
            normalized_heatmaps[k][:,_,bound] += 5
        normalized_heatmaps[k][:,bound,bound] -= 5

    x = np.arange(-bound,bound)
    energies = {}
    for t in structure_type:
        for reg in regions:
            energies[f'{reg}_{t}'] = [f'{energy_lims[f'{reg}_{t}'][0]}_max']
            prev_e = energy_lims[f'{reg}_{t}'][0]
            for i,e in enumerate(energy_lims[f'{reg}_{t}'][1:]):
                energies[f'{reg}_{t}'].append(f'{prev_e}_{e}')
                prev_e = e

    for t in structure_type:
        for reg in regions:#excluding overall
            for i,e in enumerate(energies[f'{reg}_{t}']):
                # regular heatmap; only doing 5'-start and 3'-end for local regions...
                map_name = f'{reg}_{t}_{e}_energy'
                plt.figure(figsize=(35,30))
                im = sns.heatmap(
                    normalized_heatmaps[f'{reg}_{t}'][i],
                    xticklabels=x,
                    yticklabels=x,
                    # annot=True
                )
                im.set_xticks(im.get_xticks()[::20])
                im.set_xticklabels(x[::20])
                im.set_yticks(im.get_yticks()[::20])
                im.set_yticklabels(x[::20])
                fig = im.get_figure()

                fig.savefig(f'data/{t}_heatmaps_3.4.25/{map_name}.png')
                plt.close(fig)

                csv_name = f'{reg}_{t}'
                with open(f'data/{t}_heatmaps_3.4.25/{csv_name}_max_to_min.csv', 'ab') as f:
                    np.savetxt(f,
                            np.asarray(normalized_heatmaps[f'{reg}_{t}'][i]),
                            delimiter=",")
                    np.savetxt(f,
                            np.array([[],[],[]]),
                            delimiter=",")


def Overall_and_frag_end_heatmaps(
        energy_lims: dict,
        controls: dict,
        controls_stem5: dict,
        controls_stem3: dict,
        file_name: str
    ):
    # Overall Structures - 3.4.25
    structure_type = [
        'hairpins',
        'near_hairpins',
        'repeats',
        'near_repeats',
        'mirrors',
        'complements'
    ]
    regions = [
        'Overall',
        'fragment_start',
        'fragment_end',
    ]

    stem5p_heatmaps = {}
    stem3p_heatmaps = {}
    heatmap_matrices = {}
    unique_structures = {}
    completed_genes = set()
    
    bound = 200
    for t in structure_type:
        for reg in regions:
            unique_structures[f'{reg}_{t}'] = set()
            heatmap_matrices[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            stem5p_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            stem3p_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            # for _ in range(0, 2*bound, 10):
            #     heatmap_matrices[f'{reg}_{t}'][:,bound,_] += 5
            #     heatmap_matrices[f'{reg}_{t}'][:,_,bound] += 5

            #     stem5p_heatmaps[f'{reg}_{t}'][:,bound,_] += 5
            #     stem5p_heatmaps[f'{reg}_{t}'][:,_,bound] += 5

            #     stem3p_heatmaps[f'{reg}_{t}'][:,bound,_] += 5
            #     stem3p_heatmaps[f'{reg}_{t}'][:,_,bound] += 5

            # heatmap_matrices[f'{reg}_{t}'][:,bound,bound] -= 5
            # stem5p_heatmaps[f'{reg}_{t}'][:,bound,bound] -= 5
            # stem3p_heatmaps[f'{reg}_{t}'][:,bound,bound] -= 5

    with open(file_name) as f:
        for line in f:
            tmp = ast.literal_eval(line)
            for seq,seq_dels in tmp[0].items():
                for coordinates,data_fields in seq_dels.items():
                    seq_coor = str(seq) + str(coordinates)
                    if seq_coor in completed_genes:
                        continue
                    completed_genes.add(seq_coor)

                    subcat = ['0', '1', '2', '3', '4']
                    for t in structure_type:
                        for reg in regions:
                            for sc in subcat:
                                if not data_fields[f'{reg}_{t}_energy-{sc}']:
                                    continue
                                hp = [
                                    data_fields[f'{reg}_{t}_energy-{sc}'],
                                    data_fields[f'{reg}_{t}_start-{sc}'],
                                    data_fields[f'{reg}_{t}_end-{sc}'],
                                ]
                                hp[1] += bound
                                hp[2] += bound
                                if t[-4:] == 'pins':
                                    hp.append(data_fields[f'{reg}_{t}_stem-{sc}'])
                                else:
                                    hp.append(data_fields[f'{reg}_{t}_length-{sc}'])
                                hp[1] += bound
                                hp[2] += bound
                                if hp[1] >= 2*bound or hp[2] >= 2*bound or hp[1] < 0 or hp[2] < 0:
                                    continue
                                # Avoid repeats:
                                tup = tuple(hp)
                                if tup in unique_structures[f'{reg}_{t}']:
                                    continue
                                unique_structures[f'{reg}_{t}'].add(tup)

                                en = tup[0]
                                idx = 0
                                for e in energy_lims[f'{reg}_{t}']:
                                    if en < e:
                                        idx += 1
                                # heatmap is 5'-start relative to 5'-del VS 3'-end relative to 3'-del76y7uuj8io9lp0447
                                heatmap_matrices[f'{reg}_{t}'][idx, tup[1], tup[2]] += 1

                                if tup[1]+tup[3] >= 2*bound: continue
                                stem5p_heatmaps[f'{reg}_{t}'][idx, tup[1], tup[1]+tup[3]] += 1
                                # stem3p is 5'-end vs 3'-end relative to 3'-del (positive going INTO deletion)
                                if tup[2]+tup[3] >= 2*bound: continue
                                stem3p_heatmaps[f'{reg}_{t}'][idx, tup[2], tup[2]+tup[3]] += 1
    # Normalize 2x2 case window with a 6x6 control window
    normalized_heatmaps = {}
    normalized_stem5 = {}
    normalized_stem3 = {}
    for t in structure_type:
        for reg in regions:
            normalized_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            normalized_stem5[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            normalized_stem3[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))

    for k,v in heatmap_matrices.items():
        for m in range(4):#channel
            for n in range(2, 2*bound-3):#row
                for o in range(2, 2*bound-3):#column
                    norm = 9*(1 + sum(heatmap_matrices[k][m, n:n+2, o:o+2]) / (1 + sum(controls[k][m,n:n+2,o:o+2])))
                    normalized_heatmaps[k][m, n:n+2, o:o+2] = heatmap_matrices[k][m, n:n+2, o:o+2] / norm

                    norm_stem5 = 9*(1 + sum(stem5p_heatmaps[k][m, n:n+2, o:o+2]) / (1 + sum(controls_stem5[k][m,n:n+2,o:o+2])))
                    normalized_stem5[k][m, n:n+2, o:o+2] = stem5p_heatmaps[k][m, n:n+2, o:o+2] / norm_stem5
                    
                    norm_stem3 = 9*(1 + sum(stem3p_heatmaps[k][m, n:n+2, o:o+2]) / (1 + sum(controls_stem3[k][m,n:n+2,o:o+2])))
                    normalized_stem3[k][m, n:n+2, o:o+2] = stem3p_heatmaps[k][m, n:n+2, o:o+2] / norm_stem3
        for _ in range(0, 2*bound, 10):
            normalized_heatmaps[k][:,bound,_] += 5
            normalized_heatmaps[k][:,_,bound] += 5
            
            normalized_stem5[k][:,bound,_] += 5
            normalized_stem5[k][:,_,bound] += 5

            normalized_stem3[k][:,bound,_] += 5
            normalized_stem3[k][:,_,bound] += 5        
        normalized_heatmaps[k][:,bound,bound] -= 5
        normalized_stem5[k][:,bound,bound] -= 5
        normalized_stem3[k][:,bound,bound] -= 5

    
    # For OVERALL structures not including reverse_repeats
    energies = {}
    for t in structure_type:
        for reg in regions:
            energies[f'{reg}_{t}'] = [f'{energy_lims[f'{reg}_{t}'][0]}_max']
            prev_e = energy_lims[f'{reg}_{t}'][0]
            for i,e in enumerate(energy_lims[f'{reg}_{t}'][1:]):
                energies[f'{reg}_{t}'].append(f'{prev_e}_{e}')
                prev_e = e


    x = np.arange(-bound,bound)
    # x = np.arange(-50,550)
    # y = np.arange(-50,550)

    for t in structure_type:
        for reg in regions:
            for i,e in enumerate(energies[f'{reg}_{t}']):
                # regular heatmap; only doing 5'-start and 3'-end for local regions...
                if i < 1:
                    map_name = f'{reg}_{t}_{e}_energy'

                    plt.figure(figsize=(70,60))
                    im = sns.heatmap(
                        normalized_heatmaps[f'{reg}_{t}'][i],
                        xticklabels=x,
                        yticklabels=x,
                        # annot=True
                    )
                    im.set_xticks(im.get_xticks()[::20])
                    im.set_xticklabels(x[::20])
                    im.set_yticks(im.get_yticks()[::20])
                    im.set_yticklabels(x[::20])
                    fig = im.get_figure()
                    fig.savefig(f'data/{t}_heatmaps_3.4.25/{map_name}.png')
                    plt.close(fig)

                csv_name = f'{reg}_{t}'
                with open(f'data/{t}_heatmaps_3.4.25/{csv_name}_max_to_min.csv', 'ab') as f:
                    np.savetxt(f,
                                np.asarray(normalized_heatmaps[f'{reg}_{t}'][i]),
                                delimiter=",")
                    np.savetxt(f,
                                np.array([[],[],[]]),
                                delimiter=",")

                # Stem5p and Stem3p
                for stem_map, side in [(normalized_stem5[f'{reg}_{t}'][i], '5-end'), (normalized_stem3[f'{reg}_{t}'][i], '3-end')]:
                    # if i < 1:
                    map_name = f'STEM_LENGTH_{side}_{reg}_{t}_{e}_energy'
                    plt.figure(figsize=(35,30))
                    im = sns.heatmap(
                        stem_map,
                        xticklabels=x,
                        yticklabels=x,
                        # annot=True
                    )
                    im.set_xticks(im.get_xticks()[::20])
                    im.set_xticklabels(x[::20])
                    im.set_yticks(im.get_yticks()[::20])
                    im.set_yticklabels(x[::20])
                    fig = im.get_figure()
                    fig.savefig(f'data/{t}_heatmaps_3.4.25/{map_name}.png')
                    plt.close(fig)

                    csv_name = f'STEM_LENGTH_{side}_{reg}_{t}'
                    with open(f'data/{t}_heatmaps_3.4.25/{csv_name}_max_to_min.csv', 'ab') as f:
                        np.savetxt(f,
                                    np.asarray(stem_map),
                                    delimiter=",")
                        np.savetxt(f,
                                    np.array([[],[],[]]),
                                    delimiter=",")



def main():
    completed_genes = set()
    all_freq = []


    structure_type = [
        'hairpins',
        'near_hairpins',
        'repeats',
        'near_repeats',
        'mirrors',
        'complements'
    ]
    regions = [
        '5_prime',
        '3_prime',
        'Overall',
        'fragment_start',
        'fragment_end'
    ]
    all_energies = {}
    for t in structure_type:
        for reg in regions:
            all_energies[f'{reg}_{t}'] = []
    
    
    file_name = 'data/consolidated_top_deletions.json'

    # with open('data/consolidated_top_deletions_finished_3.7.25.json') as f:
    # with open('data/consolidated_top_deletions_740.json') as f:
    # with open('data/consolidated_top_deletions_finished_3.8.25.json') as f:
    # with open('data/consolidated_top_deletions_complete_3.9.25.json') as f:
    # with open('data/consolidated_top_deletions.json') as f:
    # with open('data/consolidated_FAKE_top_deletions.json') as f:
    # with open('data/consolidated_top_deletions_fixedlongrange.json') as f:
    with open(file_name) as f:
        for line in f:
            tmp = ast.literal_eval(line)
            for seq,seq_dels in tmp[0].items():
                for coordinates,data_fields in seq_dels.items():
                    seq_coor = str(seq) + str(coordinates)
                    if seq_coor in completed_genes:
                        continue
                    completed_genes.add(seq_coor)
                    all_freq.append(data_fields['percentage_of_reads'])
                    # now go through each hairpin:
                    subcat = ['0', '1', '2', '3', '4']
                    # hairpin
                    for t in structure_type:
                        for cat in regions:
                            for sc in subcat:
                                if not data_fields[f'{cat}_{t}_energy-{sc}']:
                                    break
                                all_energies[f'{cat}_{t}'].append(data_fields[f'{cat}_{t}_energy-{sc}'])

    energy_lims = find_energy_splits(all_energies)
    control_heatmaps, control_stem5, control_stem3 = create_control_heatmaps(energy_lims, file_name)
    fiveP_threeP_heatmaps(energy_lims, control_heatmaps, file_name)
    # overall_heatmaps(energy_lims, control_heatmaps, file_name)
    Overall_and_frag_end_heatmaps(energy_lims, control_heatmaps, control_stem5, control_stem3, file_name)



if __name__ == '__main__':
    main()



'''
# Scratchwork/old functions

def overall_heatmaps(energy_lims: dict, controls: dict, file_name: str):
    # Overall Structures - 3.4.25
    heatmap_matrices = {}
    unique_structures = {}
    completed_genes = set()

    structure_type = [
        'hairpins',
        'near_hairpins',
        'repeats',
        'near_repeats',
        'mirrors',
        'complements'
    ]
    regions = [
        'Overall'
    ]

    stem5p_heatmaps = {}
    stem3p_heatmaps = {}

    for t in structure_type:
        for reg in regions:
            bound = 200
            unique_structures[f'{reg}_{t}'] = set()
            heatmap_matrices[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            stem5p_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            stem3p_heatmaps[f'{reg}_{t}'] = np.zeros((4,2*bound,2*bound))
            for _ in range(0, 2*bound, 10):
                heatmap_matrices[f'{reg}_{t}'][:,bound,_] += 5
                heatmap_matrices[f'{reg}_{t}'][:,_,bound] += 5

                stem5p_heatmaps[f'{reg}_{t}'][:,bound,_] += 5
                stem5p_heatmaps[f'{reg}_{t}'][:,_,bound] += 5

                stem3p_heatmaps[f'{reg}_{t}'][:,bound,_] += 5
                stem3p_heatmaps[f'{reg}_{t}'][:,_,bound] += 5

            heatmap_matrices[f'{reg}_{t}'][:,bound,bound] -= 5
            stem5p_heatmaps[f'{reg}_{t}'][:,bound,bound] -= 5
            stem3p_heatmaps[f'{reg}_{t}'][:,bound,bound] -= 5

    # with open('data/consolidated_top_deletions_finished_3.7.25.json') as f:
    # with open('data/consolidated_top_deletions_finished_3.8.25.json') as f:
    # with open('data/consolidated_top_deletions.json') as f:
    # with open('data/consolidated_FAKE_top_deletions.json') as f:
    with open(file_name) as f:
        missing_lengths = 0
        for line in f:
            tmp = ast.literal_eval(line)
            for seq,seq_dels in tmp[0].items():
                for coordinates,data_fields in seq_dels.items():
                    seq_coor = str(seq) + str(coordinates)
                    if seq_coor in completed_genes:
                        continue
                    completed_genes.add(seq_coor)

                    subcat = ['0', '1', '2', '3', '4']
                    for t in structure_type:
                        for reg in regions:
                            for sc in subcat:
                                if not data_fields[f'{reg}_{t}_energy-{sc}']:
                                    continue
                                hp = [
                                    data_fields[f'{reg}_{t}_energy-{sc}'],
                                    data_fields[f'{reg}_{t}_start-{sc}'],
                                    data_fields[f'{reg}_{t}_end-{sc}'],
                                ]
                                if t[-4:] == 'pins':
                                    hp.append(data_fields[f'{reg}_{t}_stem-{sc}'])
                                else:
                                    hp.append(data_fields[f'{reg}_{t}_length-{sc}'])
                                hp[1] += bound
                                hp[2] += bound
                                if hp[1] >= 2*bound or hp[2] >= 2*bound or hp[1] < 0 or hp[2] < 0:
                                    continue
                                # Avoid repeats:
                                tup = tuple(hp)
                                if tup in unique_structures[f'{reg}_{t}']:
                                    continue
                                unique_structures[f'{reg}_{t}'].add(tup)

                                en = tup[0]
                                lower_bound = energy_lims[f'{reg}_{t}'][0]
                                upper_bound = energy_lims[f'{reg}_{t}'][1]
                                if en < lower_bound:
                                    idx = 0
                                elif lower_bound <= en < upper_bound:
                                    idx = 1
                                else:#>=upper_bound
                                    idx = 2
                                # heatmap is 5'-start relative to 5'-del VS 3'-end relative to 3'-del76y7uuj8io9lp0447
                                heatmap_matrices[f'{reg}_{t}'][idx, tup[1], tup[2]] += 1
                                # stem5p is 5'-start vs 3'-start relative to 5'-del (positive going INTO deletion)
                                # if type(tup[3]) != int:
                                # if not len(tup[3]):
                                # print(f'{reg}_{t}_{sc}')
                                #     missing_lengths += 1
                                #     continue
                                if tup[1]+tup[3] >= 2*bound: continue
                                stem5p_heatmaps[f'{reg}_{t}'][idx, tup[1], tup[1]+tup[3]] += 1
                                # stem3p is 5'-end vs 3'-end relative to 3'-del (positive going INTO deletion)
                                if tup[2]+tup[3] >= 2*bound: continue
                                stem3p_heatmaps[f'{reg}_{t}'][idx, tup[2], tup[2]+tup[3]] += 1
    # For OVERALL structures not including reverse_repeats
    print('MISSING LENGTHS: ', missing_lengths, '===================================================')
    energies = {}
    for t in structure_type:
        energies[f'Overall_{t}'] = [f'{energy_lims[f'Overall_{t}'][0]}_max']
        prev_e = energy_lims[f'Overall_{t}'][0]
        for i, e in enumerate(energy_lims[f'Overall_{t}'][1:]):
            energies[f'Overall_{t}'].append(f'{prev_e}_{e}')
            prev_e = e


    x = np.arange(-bound,bound)
    # y = np.arange(-50,550)

    for t in structure_type:
        for reg in regions:
            for i,e in enumerate(energies[f'{reg}_{t}']):
                # regular heatmap; only doing 5'-start and 3'-end for local regions...
                if i < 1:
                    map_name = f'{reg}_{t}_{e}_energy'

                    plt.figure(figsize=(70,60))
                    im = sns.heatmap(
                        heatmap_matrices[f'{reg}_{t}'][i],
                        xticklabels=x,
                        yticklabels=x,
                        # annot=True
                    )
                    im.set_xticks(im.get_xticks()[::20])
                    im.set_xticklabels(x[::20])
                    im.set_yticks(im.get_yticks()[::20])
                    im.set_yticklabels(x[::20])
                    fig = im.get_figure()
                    fig.savefig(f'data/{t}_heatmaps_3.4.25/{map_name}.png')
                    plt.close(fig)

                csv_name = f'Overall_{t}'
                with open(f'data/{t}_heatmaps_3.4.25/{csv_name}_max_to_min.csv', 'ab') as f:
                    np.savetxt(f,
                                np.asarray(heatmap_matrices[f'{reg}_{t}'][i]),
                                delimiter=",")
                    np.savetxt(f,
                                np.array([[],[],[]]),
                                delimiter=",")

                # Stem5p and Stem3p
                for stem_map, side in [(stem5p_heatmaps[f'{reg}_{t}'][i], '5-end'), (stem3p_heatmaps[f'{reg}_{t}'][i], '3-end')]:
                    # if i < 1:
                    map_name = f'STEM_LENGTH_{side}_{reg}_{t}_{e}_energy'
                    plt.figure(figsize=(35,30))
                    im = sns.heatmap(
                        stem_map,
                        xticklabels=x,
                        yticklabels=x,
                        # annot=True
                    )
                    im.set_xticks(im.get_xticks()[::20])
                    im.set_xticklabels(x[::20])
                    im.set_yticks(im.get_yticks()[::20])
                    im.set_yticklabels(x[::20])
                    fig = im.get_figure()

                    fig.savefig(f'data/{t}_heatmaps_3.4.25/{map_name}.png')
                    plt.close(fig)

                    csv_name = f'{reg}_{t}_{side}'
                    with open(f'data/{t}_heatmaps_3.4.25/{csv_name}_max_to_min.csv', 'ab') as f:
                        np.savetxt(f,
                                    np.asarray(stem_map),
                                    delimiter=",")
                        np.savetxt(f,
                                    np.array([[],[],[]]),
                                    delimiter=",")

'''