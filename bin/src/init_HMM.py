#!/usr/bin/env python
import pandas as pd 
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import json
import yaml
import argparse
import warnings

from distribution import param, meth_param

def init_k_means(data, k, marker, distributions, file, methylation=None):
    json_data = {}

    json_data["states"] = k
    
    json_data["marker"] = marker

    json_data["methylation"] = False
    if 'dna_methylation' in marker:
        json_data["methylation"] = True
        json_data["marker"] = marker[:-1]

    json_data["emission"] = []
    
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(data)
    try:
        kmeans = KMeans(init="random", n_clusters=k, n_init=10, max_iter=300, random_state=42)
        kmeans.fit(scaled_features)

        mean = np.empty((k, data.shape[1]))
        std = np.empty((k, data.shape[1]))
        zero = np.empty((k, data.shape[1]))
        n = np.empty((k, data.shape[1]))
        for i in range(k):
            subset = data[kmeans.labels_==i]
            values = subset.values
            for j in range(data.shape[1]):
                mean[i, j] = values[:, j].mean()
                std[i, j] = values[:, j].std()
                zero[i, j] = np.count_nonzero(values[:, j]==0) / len(values[:, j])
                n[i, j] = max(values[:, j])

        if json_data["methylation"]:
            mean_DNA, std_DNA, zero_DNA = np.empty(k), np.empty(k), np.empty(k)
            for i in range(k):
                subset = methylation[kmeans.labels_==i]
                mean_DNA[i] = subset.mean(axis=0, skipna=True).to_list()[0]
                std_DNA[i] = subset.std(axis=0, skipna=True).to_list()[0]
                zero_DNA[i] = subset.isna().sum() / subset.shape[0]
        
        maxN = np.empty(data.shape[1])
        for i in range(data.shape[1]):
            maxN[i] = max(data.values[:, i])

        m = len(json_data["marker"])
        for i in range(k):
            state = []
            for j in range(m):
                state.append({"distribution": distributions[j], "parameters": param(distributions[j], mean[i, j], std[i, j], zero[i, j], n[i, j], maxN[j])})
            if json_data["methylation"]:
                state.append({"distribution": distributions[m], "parameters": meth_param(distributions[m], mean_DNA[i], std_DNA[i], zero_DNA[i])})
            json_data["emission"].append(state)
    except:
        m = len(json_data["marker"])
        for i in range(k):
            state = []
            for j in range(m):
                state.append({"distribution": 'PO', "parameters": {"lambda": 1}})
            if json_data["methylation"]:
                state.append({"distribution": 'BI', "parameters": {"p": 0.5}})
            json_data["emission"].append(state)
    
    with open(file, 'w') as jsonFile:
        json.dump(json_data, jsonFile, indent=4)

def init_k_means_meth(data, k, distribution, file):
    json_data = {}
    json_data["states"] = k
    json_data["marker"] = []
    json_data["emission"] = []
    json_data["methylation"] = True
    
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(data)
    try:
        kmeans = KMeans(init="random", n_clusters=k, n_init=10, max_iter=300, random_state=42)
        kmeans.fit(scaled_features)

        for i in range(k):
            subset = data[kmeans.labels_==i]
            mean = subset.mean(axis=0, skipna=True).to_list()[0]
            std = subset.std(axis=0, skipna=True).to_list()[0]
            zero = subset.isna().sum() / subset.shape[0]
            json_data["emission"].append([{"distribution": distribution, "parameters": meth_param(distribution, mean, std, zero)}])
    except:
        for i in range(k):
            state = []
            state.append([{"distribution": 'BI', "parameters": {"p": 0.5}}])
            json_data["emission"].append(state)
    
    with open(file, 'w') as jsonFile:
        json.dump(json_data, jsonFile, indent=4)

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)

    parser = argparse.ArgumentParser(description = "Initialize HMM.")
    parser.add_argument("-d", metavar = "data", type = str, nargs = 1, default=[''], help = "Input data.")
    parser.add_argument("-e", metavar = "methylation-data", type = str, nargs = 1, default='', help = "Optional methylation data.")
    parser.add_argument("-m", metavar = "marker", type = str, nargs = 1, help = "YAML file with marker information (name and distributional assumption).")
    parser.add_argument("-j", metavar = "json", type = str, nargs = 1, help = "Output json file with initial HMM parameters.")
    args = parser.parse_args()

    try:
        markerFile = args.m[0]
        dataFile = args.d[0]
        out = args.j[0]
    except:
        parser.print_help()
        return
        
    marker = dict()
    with open(markerFile) as file:
        marker = yaml.safe_load(file)
    
    if dataFile != '':
        k = marker['states']
        m = marker['marker']
        data = pd.read_csv(dataFile, sep='\t', header=None)
        names = []
        distribution = []
        for i in range(m):
            distribution.append(marker['marker_spec'][i]['distribution'])
            names.append(marker['marker_spec'][i]['name'])

        if 'dna_methylation' in marker:
            meth =  pd.read_csv(args.e[0], sep='\t', names = ['Cov', 'Meth'])
            meth['prop'] = meth['Meth'] / meth['Cov']
            meth = meth.drop(['Meth', 'Cov'], axis=1)
            distribution.append(marker['dna_methylation'])
            names.append('dna_methylation')
            init_k_means(data, k, names, distribution, out, meth)
        else:
            init_k_means(data, k, names, distribution, out)
    else:
        k = marker['states']
        meth =  pd.read_csv(args.e[0], sep='\t', names = ['Cov', 'Meth'])
        distribution = marker['distribution']
        meth['prop'] = meth['Meth'] / meth['Cov']
        meth.drop(['Meth', 'Cov'], axis=1, inplace=True)
        meth = meth.loc[~meth['prop'].isna()]
        init_k_means_meth(meth, k, distribution, out)

if __name__ == "__main__":
    main()

