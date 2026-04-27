#!/usr/bin/env python
import pandas as pd
import yaml 
import argparse
import warnings
import os

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description = "Create count matrices.")
    parser.add_argument("-d", metavar = "model", type = str, nargs = 1, help = "Yaml file with path to data files.")
    parser.add_argument("-c", metavar = "train", type = str, nargs = 1, help = "Output file.")
    parser.add_argument("-m", metavar = "trainMeth", type = str, nargs = 1, help = "Output file for methylation.")
    parser.add_argument("-r", metavar = "regions", type = str, nargs = 1, help = "Output for index to indicate start of separate regions.")
    args = parser.parse_args()

    try:
        data = args.d[0]
        trainFile = args.c[0]
        trainMethFile = args.m[0]
        regionFile = args.r[0]
    except:
        parser.print_help()
        return

    marker = dict()
    with open(data) as file:
        marker = yaml.safe_load(file)
    files = marker['data']

    chr = ['pilot_hg38']
    if 'chr' in marker:
        chr = marker['chr']
        if type(chr) is list:
            chr = [str(x) for x in chr]
        else:
            chr = [str(chr)]

    methFiles = []
    if 'meth_data' in marker:
        methFiles = marker['meth_data']
        if (len(files) != len(methFiles)):
            raise Exception('Number of histone counts and methylation must correspond!')

    columns = [str(mark['name']) for mark in marker['marker_spec']]

    train = []
    trainMeth = []
    index =  []
    start = 0
    
    # Dynamically find the directory where this script lives
    script_dir = os.path.dirname(os.path.realpath(__file__))

    for i in range(len(files)):
        data = pd.read_csv(files[i], sep='\t', converters={0:str})
        data.columns = ['chr', 'start', 'end'] + list(data.columns)[3:]

        if chr[0].startswith('pilot'):
            if chr[0] == 'pilot_hg19':
                pilot_path = os.path.join(script_dir, "encode_pilot_regions", "hg19.bed")
                pilot = pd.read_csv(pilot_path, sep="\t", names=['chr', 'start', 'end', 'name'])
            else:
                pilot_path = os.path.join(script_dir, "encode_pilot_regions", "hg38.bed")
                pilot = pd.read_csv(pilot_path, sep="\t", names=['chr', 'start', 'end', 'name'])

            rows = 0
            for row in pilot['chr'].index:
                chr_data = data[data["chr"] == pilot.loc[row]['chr']]
                chr_data = chr_data.loc[(chr_data['start'] >= pilot.loc[row]['start']) & (chr_data['end'] < pilot.loc[row]['end'])]
                chr_data = chr_data.iloc[:, 3:]
                chr_data = pd.DataFrame(chr_data, columns=columns)
                train.append(chr_data)

                if 'dna_methylation' in marker:
                    methData = pd.read_csv(methFiles[i], sep='\t')
                    methData = methData.iloc[chr_data.index[0]:chr_data.index[-1]+1,:]
                    trainMeth.append(methData)
                
                rows += chr_data.shape[0] 

            index.append(start)
            start += rows

        else:
            for c in chr:
                chr_data = data[data["chr"] == c]
                chr_data = chr_data.iloc[:, 3:]
                chr_data = pd.DataFrame(chr_data, columns=columns)
                train.append(chr_data)

                if 'dna_methylation' in marker:
                    methData = pd.read_csv(methFiles[i], sep='\t')
                    methData = methData.iloc[chr_data.index[0]:chr_data.index[-1]+1,:]
                    trainMeth.append(methData)

                index.append(start)
                start += chr_data.shape[0]
            
    train = pd.concat(train, axis=0)
    train.to_csv(trainFile, sep='\t', header=False, index=False)

    regions=open(regionFile,'w')
    for r in index:
        regions.write(str(r))
        regions.write('\n')
    regions.close()

    if 'dna_methylation' in marker:
        trainMeth = pd.concat(trainMeth, axis=0)
        trainMeth.to_csv(trainMethFile, sep='\t', header=False, index=False)
    else:
        with open(trainMethFile, 'w') as _:
            pass

if __name__ == "__main__":
    main()