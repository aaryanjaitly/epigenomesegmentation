#!/usr/bin/env python
import pandas as pd
import yaml 
import argparse
import os
import warnings

def all_counts(data, counts, regions):
    series = data['chr'].ne(data['chr'].shift())
    index = series[series].index

    with open(regions, 'w') as file:
        for r in index:
            file.write(str(r))
            file.write('\n')

    data = data.iloc[:, 3:]
    data.to_csv(counts, sep='\t', header=False, index=False)

def train_counts(config, data, counts, regions):
    chr = ['pilot_hg38']
    if 'chr' in config:
        chr = config['chr']
        if type(chr) is list:
            chr = [str(x) for x in chr]
        else:
            chr = [str(chr)]

    train = []
    
    # FIX: Dynamically find the script directory to locate pilot beds
    script_dir = os.path.dirname(os.path.realpath(__file__))
    
    if chr[0].startswith('pilot'):
        if chr[0] == 'pilot_hg19':
            pilot_path = os.path.join(script_dir, "encode_pilot_regions", "hg19.bed")
            pilot = pd.read_csv(pilot_path, sep="\t", names=['chr', 'start', 'end', 'name'])
        else:
            pilot_path = os.path.join(script_dir, "encode_pilot_regions", "hg38.bed")
            pilot = pd.read_csv(pilot_path, sep="\t", names=['chr', 'start', 'end', 'name'])

        for row in pilot['chr'].index:
            chr_data = data[data["chr"] == pilot.loc[row]['chr']]
            chr_data = chr_data.loc[(chr_data['start'] >= pilot.loc[row]['start']) & (chr_data['end'] < pilot.loc[row]['end'])]
            chr_data = chr_data.iloc[:, 3:]
            train.append(chr_data)

    else:
        for c in chr:
            chr_data = data[data["chr"] == c]
            chr_data = chr_data.iloc[:, 3:]
            train.append(chr_data)
            
    train = pd.concat(train, axis=0)
    train.to_csv(counts, sep='\t', header=False, index=False)

    with open(regions, 'w') as file:
            file.write('0')

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description = "Create count matrices.")
    parser.add_argument("-d", metavar = "config", type = str, nargs = 1, help = "Yaml file with path to data file.")
    parser.add_argument("-c", metavar = "trainCounts", type = str, nargs = 1, help = "Output file with train counts.")
    parser.add_argument("-r", metavar = "trainRegions", type = str, nargs = 1, help = "Output file train regions.")
    parser.add_argument("-C", metavar = "counts", type = str, nargs = 1, help = "Output file with counts.")
    parser.add_argument("-R", metavar = "regions", type = str, nargs = 1, help = "Output file with regions.")
    args = parser.parse_args()

    try:
        config = args.d[0]
        trainCounts = args.c[0]
        trainRegions = args.r[0]
        counts = args.C[0]
        regions = args.R[0]
    except:
        parser.print_help()
        exit()

    conifg_dict = dict()
    with open(config) as file:
        conifg_dict = yaml.safe_load(file)
    data_file = conifg_dict['data']

    # Unpack from list if necessary (since YAML parses '["/path"]' as a list)
    if isinstance(data_file, list):
        data_file = data_file[0]

    data = pd.read_csv(data_file, sep='\t', converters={0:str, 1:int})
    data.columns = ['chr', 'start', 'end', 'Cov', 'Meth']

    all_counts(data, counts, regions)
    train_counts(conifg_dict, data, trainCounts, trainRegions)

if __name__ == "__main__":
    main()