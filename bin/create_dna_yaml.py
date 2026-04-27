#!/usr/bin/env python
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Generate YAML config for DNA_PREPARE.")
    parser.add_argument("-p", "--prefix", required=True, help="Output prefix for the YAML file")
    parser.add_argument("-m", "--meth", required=True, nargs='+', help="Path(s) to methylation data file(s)")
    parser.add_argument("-s", "--state", required=True, help="Number of states")
    parser.add_argument("--dist_meth", required=True, help="Default methylation distribution")
    parser.add_argument("--overrides", required=False, default="", help="Comma-separated overrides (e.g. WGBS:B)")
    args = parser.parse_args()

    # Parse distribution overrides
    overrides = {}
    if args.overrides:
        for pair in args.overrides.split(','):
            if ':' in pair:
                k, v = pair.split(':', 1)
                overrides[k] = v

    # Fallback to default if no specific override exists for WGBS
    dist = overrides.get("WGBS", args.dist_meth)

    # Handle file path formatting (single string vs list)
    abs_meth_paths = [os.path.abspath(m) for m in args.meth]
    if len(abs_meth_paths) == 1:
        data_str = f'"{abs_meth_paths[0]}"'
    else:
        data_str = "[" + ", ".join([f'"{p}"' for p in abs_meth_paths]) + "]"

    # Build the YAML content
    yaml_lines = [
        f"distribution: {dist}",
        f"states: {args.state}",
        f"data: {data_str}",
        "marker: WGBS"
    ]

    # Write to file
    with open(f"{args.prefix}.yaml", 'w') as f:
        f.write("\n".join(yaml_lines) + "\n")

if __name__ == "__main__":
    main()