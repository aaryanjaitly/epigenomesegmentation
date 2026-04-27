#!/usr/bin/env python
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Generate YAML config for STD_PREPARE.")
    parser.add_argument("-p", "--prefix", required=True, help="Output prefix for the YAML file")
    parser.add_argument("-i", "--histone", required=True, help="Path to histone data file")
    parser.add_argument("-m", "--meth", required=False, default="", help="Path to methylation data file (optional)")
    parser.add_argument("-s", "--state", required=True, help="Number of states")
    parser.add_argument("-c", "--chr_params", required=True, help="Chromosome parameters")
    parser.add_argument("--dist_hist", required=True, help="Default histone distribution")
    parser.add_argument("--dist_meth", required=True, help="Default methylation distribution")
    parser.add_argument("--overrides", required=False, default="", help="Comma-separated overrides (e.g. H3K4me3:NBI,WGBS:B)")
    args = parser.parse_args()

    # Parse distribution overrides
    overrides = {}
    if args.overrides:
        for pair in args.overrides.split(','):
            if ':' in pair:
                k, v = pair.split(':', 1)
                overrides[k] = v

    # Extract markers from the first row (assuming first 3 cols are chr, start, end)
    with open(args.histone, 'r') as f:
        header = f.readline().strip().split()
        markers = header[3:]

    # Build the YAML content dynamically
    yaml_lines = [
        f"states: {args.state}",
        f"marker: {len(markers)}",
        "marker_spec:"
    ]

    for mark in markers:
        dist = overrides.get(mark, args.dist_hist)
        yaml_lines.append(f"  - name: {mark}")
        yaml_lines.append(f"    distribution: {dist}")

    yaml_lines.append(f"data: [{os.path.abspath(args.histone)}]")
    yaml_lines.append(f"chr: {args.chr_params}")

    # Append Methylation data if the file exists and is not empty
    if args.meth and os.path.exists(args.meth) and os.path.getsize(args.meth) > 0:
        meth_dist = overrides.get("WGBS", args.dist_meth)
        yaml_lines.append(f"dna_methylation: {meth_dist}")
        yaml_lines.append(f"meth_data: [{os.path.abspath(args.meth)}]")

    # Write to file
    with open(f"{args.prefix}.yaml", 'w') as f:
        f.write("\n".join(yaml_lines) + "\n")

if __name__ == "__main__":
    main()