#!/usr/bin/env python
import argparse
import os
import yaml
import warnings

import pandas as pd
import seaborn as sb
import matplotlib.colors as mcolors

def combine_data(model, input, index, columns, state_column):
    prefix = os.path.basename(model['data'][index]) 

    combined = []
    data = pd.read_csv(model['data'][index], sep='\t', converters={0:str})

    data.columns = ['chr', 'start', 'end'] + list(data.columns)[3:]
    combined.append(pd.DataFrame(data, columns=columns))

    if 'dna_methylation' in model:
        combined.append(pd.read_csv(model['meth_data'][index], sep='\t'))

    if state_column == 'viterbi':
        combined.append(pd.read_csv(input+state_column+'_'+prefix, names=[state_column], sep='\t'))
    elif state_column == 'posterior':
        combined.append(pd.read_csv(input+state_column+'_'+prefix, names=[state_column, 'score'], sep='\t'))
    else:
        print("Unknown column")

    return pd.concat(combined, axis = 1), prefix

def to_bed(data, out_dir, prefix, column):
    data = data.set_index('chr')
    combined = []
    for chr in set(data.index):
        subset = data.loc[[chr]]
        label_groups = subset[column].ne(subset[column].shift()).cumsum()
        subset = (subset.groupby(label_groups).agg({'start':'min', 'end':'max', column:'first'}).reset_index(drop=True))
        subset.index = [chr] * subset.shape[0]
        combined.append(subset)

    combined = pd.concat(combined, ignore_index=False)

    combined['score'] = 0
    combined['strand'] = '.'
    combined['thickStart'] = combined['start']
    combined['thickEnd'] = combined['end']

    s = max(combined[column])
    #alternatively use "cubehelix" or "husl"
    cmap = sb.color_palette("hls", s)
    cmap = [mcolors.to_hex(c) for c in cmap]
    cmap = [tuple(int(hex[i+1:i+3], 16) for i in (0, 2, 4)) for hex in cmap]

    combined['rgb']=[','.join(map(str, cmap[x-1])) for x in combined[column]]
    combined.index = 'chr' + combined.index
    combined.to_csv(out_dir+column+'_'+os.path.splitext(prefix)[0]+'.bed.gz', sep='\t', header=False)

    if column == 'posterior':
        bedGraph = data.reset_index().loc[:,['chr', 'start', 'end', 'score']]
        bedGraph['chr'] = 'chr' + data.index
        bedGraph.to_csv(out_dir+os.path.splitext(prefix)[0]+'.bedGraph.gz', sep='\t', header=False, index=False)

def to_bed_methylation(model, output, sequence, column):
    data = pd.read_csv(model['data'], sep='\t', converters={0:str})

    if column == 'viterbi':
        data[column] = pd.read_csv(sequence, names=[column], sep='\t')[column]
    elif column == 'posterior':
        decoding = pd.read_csv(sequence, names=[column, 'score'], sep='\t')
        data[column] = decoding[column]
        data['score'] = decoding['score']
    else:
        print("Unknown column")
        return

    prefix = os.path.basename(model['data']) 
    data.to_csv(args.o[0]+prefix, index=False, sep='\t')
    to_bed(data, output, prefix, column)

warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser(description = "Create bed files from segmentation and data.")
parser.add_argument("-d", metavar = "model", type = str, nargs = 1, help = "Yaml file with path to data files.")
parser.add_argument("-o", metavar = "output", type = str, nargs = 1, help = "Output directory for bed files.")
parser.add_argument("-i", metavar = "input", type = str, nargs = 1, help = "Input directory with state sequences.")
parser.add_argument("-c", metavar = "columns", type = str, nargs = 1, help = "Column names for state sequence.")
args = parser.parse_args()

if not os.path.exists(args.o[0]):
    os.makedirs(args.o[0])

model = dict()
with open(args.d[0]) as file:
    model = yaml.safe_load(file)

if 'marker_spec' in model:
    columns = ['chr', 'start', 'end'] + [mark['name'] for mark in model['marker_spec']]
    state_columns = args.c[0]

    for i in range(len(model['data'])):
        df, prefix = combine_data(model, args.i[0], i, columns, state_columns)
        if i == 0:
            df.to_csv(args.o[0]+prefix, index=False, sep='\t')
        to_bed(df, args.o[0], prefix, state_columns)
else:
    to_bed_methylation(model, args.o[0], args.i[0], args.c[0])
