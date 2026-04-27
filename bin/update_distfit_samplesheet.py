#!/usr/bin/env python
import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser(description="Update samplesheet with best distribution fits.")
    parser.add_argument("-l", "--log_likelihoods", required=True, help="Path to log_likelihoods.txt")
    parser.add_argument("-c", "--original_csv", required=True, help="Path to original samplesheet CSV")
    parser.add_argument("-s", "--sample_id", required=True, help="Sample ID to process")
    parser.add_argument("-o", "--output", required=True, help="Output CSV filename")
    args = parser.parse_args()

    # Parse likelihoods and map best distributions
    if os.path.exists(args.log_likelihoods) and os.path.getsize(args.log_likelihoods) > 0:
        scores = pd.read_csv(args.log_likelihoods, sep='\t', names=['Dist', 'Mark', 'Score'])
        scores['Score'] = pd.to_numeric(scores['Score'], errors='coerce')
        best_dists = scores.loc[scores.groupby('Mark')['Score'].idxmax()].set_index('Mark')['Dist'].to_dict()
    else:
        best_dists = {}

    # Update original samplesheet
    df = pd.read_csv(args.original_csv)
    df_sample = df[df['sample_id'].astype(str) == args.sample_id].copy()

    # Map the new distributions, falling back to original if no new one exists
    df_sample['distribution'] = df_sample['epigenetic_mark'].map(best_dists).fillna(df_sample['distribution'])
    
    # Save the updated samplesheet
    df_sample.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()