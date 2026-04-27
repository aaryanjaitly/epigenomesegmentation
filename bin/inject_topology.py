#!/usr/bin/env python
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description="Inject topology mapping into HMM pre-initialization JSON.")
    parser.add_argument("-i", "--input", required=True, help="Path to preinit JSON")
    parser.add_argument("-o", "--output", required=True, help="Path for output init JSON")
    args = parser.parse_args()

    # Load pre-initialized HMM
    with open(args.input, 'r') as f:
        hmm = json.load(f)
    
    # Inject topology
    N = int(hmm.get('states', 0))
    hmm['topology'] = {str(s+1): [s] for s in range(N)}
    
    # Save initialized HMM
    with open(args.output, 'w') as f:
        json.dump(hmm, f, indent=4)

if __name__ == "__main__":
    main()