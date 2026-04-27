#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import warnings

def set_plot_style():
    font = 'serif'
    tex_fonts = {
        # Use LaTeX to write all text
        'text.usetex': True,
        'font.family': font,
        'axes.labelsize': 26,
        'font.size': 26,
        'legend.fontsize': 26,
        'xtick.labelsize': 26,
        'ytick.labelsize': 26,
        'text.color':'black'
        }
    plt.rcParams.update(tex_fonts)

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    parser = argparse.ArgumentParser(description = "Plot state colors.")
    parser.add_argument("-d", metavar = "data", type = str, nargs = 1, help = "Segmentation.")
    parser.add_argument("-o", metavar = "output", type = str, nargs = 1, help = "Output filename.")
    args = parser.parse_args()

    data = pd.read_csv(args.d[0], sep='\t', header=None, skiprows=[0], converters={0:str})
    states = list(set(data[3]))

    palette = np.zeros((len(states),3), dtype=int)
    for i in range(len(states)):
        palette[i,:] = [int(x) for x in data[data[3]==states[i]].iloc[0,8].split(',')]

    # set_plot_style()

    I = np.array([[x for x in range(len(states))]], dtype=int)
    RGB = palette[I]
    plt.imshow(RGB)
    plt.yticks([])
    plt.xticks(list(range(0,len(states))), [str(x) for x in range(1, len(states)+1)])

    plt.savefig(args.o[0], bbox_inches = 'tight')

if __name__ == "__main__":
    main()