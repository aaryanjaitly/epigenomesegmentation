#!/usr/bin/env python
import pandas as pd
import yaml 
import argparse
import os
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser(description = "Create count matrices for several files.")
parser.add_argument("-d", metavar = "model", type = str, nargs = 1, help = "Yaml file with path to data files.")
parser.add_argument("-o", metavar = "output", type = str, nargs = 1, help = "Output directory.")
args = parser.parse_args()

try:
    data = args.d[0]
    out_dir = args.o[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
except:
    parser.print_help()
    exit()

marker = dict()
with open(data) as file:
    marker = yaml.safe_load(file)
files = marker['data']

methFiles = []
if 'meth_data' in marker:
    methFiles = marker['meth_data']
    if (len(files) != len(methFiles)):
        raise Exception('Number of histone counts and methylation must correspond!')

columns = [str(mark['name']) for mark in marker['marker_spec']]
for i in range(len(files)):
    data = pd.read_csv(files[i], sep='\t', converters={0:str})
    data.columns = ['chr', 'start', 'end'] + list(data.columns)[3:]

    series = data['chr'].ne(data['chr'].shift())
    index = series[series].index
    filename = os.path.basename(files[i])
    regions=open(out_dir+'/regions_'+filename,'w')
    for r in index:
        regions.write(str(r))
        regions.write('\n')
    regions.close()

    data = data.iloc[:, 3:]
    data = pd.DataFrame(data, columns=columns)
    data.to_csv(out_dir+'/counts_'+filename, sep='\t', header=False, index=False)

    if len(methFiles) > 0:
        methData = pd.read_csv(methFiles[i], sep='\t')
        methData.to_csv(out_dir+"/meth_"+filename, sep='\t', header=False, index=False)
    else:
        with open(out_dir+"/meth_"+filename, 'w') as _:
            pass