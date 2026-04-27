#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import yaml
import warnings

from plots import plot_correlation, plot_histogram, plot_methylation, set_plot_style

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description = "Preprocess chromatin data.")
    parser.add_argument("-d", metavar = "data", type = str, nargs = 1, help = "YAML file with input file path.")
    parser.add_argument("-p", metavar = "histogram plot", type = str, nargs = 1, help = "Save histogram.")
    parser.add_argument("-c", metavar = "correlation plot", type = str, nargs = 1, help = "Save correlation plot.")
    parser.add_argument("-m", metavar = "methylation density", type = str, nargs = 1, help = "Save methylation density plot.")
    args = parser.parse_args()

    try:
        yamlFile = args.d[0]
        histogram= args.p[0]
        correlation = args.c[0]
        methDensity = args.m[0]
    except:
        parser.print_help()
        return

    marker = dict()
    with open(yamlFile) as file:
        marker = yaml.safe_load(file)
    files = marker['data']
    columns = [str(mark['name']) for mark in marker['marker_spec']]

    set_plot_style()

    data = pd.read_csv(files[0], sep='\t', converters={0:str})
    counts = data.iloc[:, 3:]
    counts = pd.DataFrame(counts, columns=columns)
    plot_correlation(counts, correlation)

    m = len(columns)
    _, axes = plt.subplots(m, len(files), figsize=(16 * len(files), 6*m), squeeze=False)
    for i in range(len(files)):
        data = pd.read_csv(files[i], sep='\t', converters={0:str})
        max = int(math.ceil(np.percentile(data.iloc[:, 3:].values, 99) / 100.0)) * 100
        counts = data.iloc[:, 3:]
        counts = pd.DataFrame(counts, columns=columns)

        for k in range(m):
            ax = axes[k][i]
            plot_histogram(ax, max, counts.iloc[:, k], 50)
    plt.savefig(histogram, bbox_inches = 'tight')

    if 'meth_data' in marker:
        files = marker['meth_data']
        _, axes = plt.subplots(1, len(files), figsize=(12 * len(files), 10), squeeze=False)
        for i in range(len(files)):
            data = pd.read_csv(files[i], sep='\t')
            data['prop'] = data['Meth'] / data['Cov']
            data['prop'] = data['prop'].fillna(-0.1)
            plot_methylation(axes[0][i], data)
        plt.savefig(methDensity, bbox_inches = 'tight')

if __name__ == "__main__":
    main()