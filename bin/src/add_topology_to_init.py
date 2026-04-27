#!/usr/bin/env python
import json
from argparse import ArgumentParser


def main(args):
    HMM = dict()
    with open(args.input) as file:
        HMM = json.load(file)

    N = int(HMM["states"])
    topology = dict()
    for s in range(N):
        topology[str(s + 1)] = [s]
    HMM["topology"] = topology

    with open(args.output, 'w') as jsonFile:
        json.dump(HMM, jsonFile, indent=4)


if __name__ == '__main__':
    p = ArgumentParser(description = "Add topology to initial HMM.")
    p.add_argument("--input", help = "Json file with inital parameters.")
    p.add_argument("--output", help = "Output Json file with topology.")
    main(p.parse_args())