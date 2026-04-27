#!/usr/bin/env python
import argparse
import math
import numpy as np
import pandas as pd
import json
import warnings

import plots

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)

    parser = argparse.ArgumentParser(description = "Results of chromatin segmentation.")
    parser.add_argument("-c", metavar = "data", type = str, nargs = 1, help = "Data with histone modifcations and segmentation.")
    parser.add_argument("-j", metavar = "json", type = str, nargs = 1, help = "Json file with HMM parameters.")
    parser.add_argument("-t", metavar = "transtion_matrix", type = str, nargs = 1, help = "Plot transition matrix.")
    parser.add_argument("-e", metavar = "mean_emission", type = str, nargs = 1, help = "Plot mean emission.")
    parser.add_argument("-n", metavar = "norm_emission", type = str, nargs = 1, help = "Plot mean emission normalized per marker.")
    parser.add_argument("-m", metavar = "state_membership", type = str, nargs = 1, help = "Plot histogram of state memberships.")
    parser.add_argument("-l", metavar = "state_length", type = str, nargs = 1, help = "Plot state length distribution.")
    parser.add_argument("-s", metavar = "states", type = str, nargs = 1, help = "Column name with state sequence.")
    parser.add_argument("-d", metavar = "bed_data", type = str, nargs = 1, help = "Segmentation BED file for colors.")
    
    args = parser.parse_args()

    try:
        data = args.c[0]
        model = args.j[0]
        transitionMatrix = args.t[0]
        meanEmission = args.e[0]
        normEmission = args.n[0]
        stateMembership = args.m[0]
        stateLength = args.l[0]
        column = args.s[0]
        bed_data = pd.read_csv(args.d[0], sep='\t', header=None, skiprows=[0], converters={0:str})
        bed_states = list(set(bed_data[3]))

    except:
        parser.print_help()
        return

    HMM = dict()
    with open(model) as file:
        HMM = json.load(file)
    
    HMM['states'] = int(HMM['states'])
    # Safely handle the methylation flag
    HMM['methylation'] = True if HMM.get('methylation') == "true" or HMM.get('methylation') == True else False
    HMM['initial_state_distribution'] = np.array(HMM['initial_state_distribution'], dtype=float)
    HMM['transition_matrix'] = np.array(HMM['transition_matrix'], dtype=float)

    # Safely get marker length even if the key is missing in DNA-only runs
    m = len(HMM.get('marker', []))
    data = pd.read_csv(data, sep='\t', converters={0:str})
    data.columns = ['chr', 'start', 'end'] + list(data.columns)[3:]
    
    # Safe calculation for percentile (prevents IndexError on empty arrays)
    if m > 0:
        max_val = int(math.ceil(np.percentile(data.iloc[:,3:3+m].values, 99.5) / 100.0)) * 100
    else:
        max_val = 100 

    # DNA-Methylation logic
    if HMM['methylation']:
        data['DNA-Methylation'] = data['Meth'] / data['Cov']
        data = data.drop(columns=['Meth', 'Cov'])
        columns = list(data.columns)
        new_cols = columns[:m+3] + ['DNA-Methylation'] + columns[m+3:len(columns)-1]
        data=data[new_cols]

    start_states = range(0, HMM['states'])
    end_states = range(0, HMM['states'])

    if 'topology' in HMM:
        start_states = [int(HMM['topology'][str(i)][0]) for i in range(1, HMM['states']+1)]
        end_states = [int(HMM['topology'][str(i)][-1]) for i in range(1, HMM['states']+1)]

    transition_matrix = np.zeros((HMM['states'], HMM['states']))
    for i in range(HMM['states']):
        for j in range(HMM['states']):
            if i == j:
                transition_matrix[i, j] = HMM['transition_matrix'][end_states[i], end_states[i]]
            else:
                transition_matrix[i, j] = HMM['transition_matrix'][end_states[i], start_states[j]]

    # Palette generation
    palette = np.zeros((len(bed_states),3), dtype=float)
    for i in range(len(bed_states)):
        rgb = [int(x) for x in bed_data[bed_data[3]==bed_states[i]].iloc[0,8].split(',')]
        palette[i, :] = [x / 255.0 for x in rgb]

    plots.set_plot_style()
    plots.plot_mean_emission(data.iloc[:, 3:], column, m, HMM['states'], meanEmission, palette)
    plots.plot_state_memberships(data, column, HMM, stateMembership, palette)
    plots.plot_transition(transition_matrix, transitionMatrix, palette)
    plots.plot_state_length(data, column, stateLength, palette)
    plots.plot_mean_emission_norm(data.iloc[:, 3:], column, m, HMM['states'], normEmission, palette)

    if HMM['methylation']:
        file = meanEmission[:-16]
        plots.plot_average_methylation(data, column, file+'averageMethylation.png', palette)

if __name__ == "__main__":
    main()