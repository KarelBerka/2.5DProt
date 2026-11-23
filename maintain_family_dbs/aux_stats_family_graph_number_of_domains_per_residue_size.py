#!/usr/bin/env python3
# Usage:
# ./aux_stats_family_graph_number_of_domains_per_residue_size.py --family 1.10.630.10 --filename generated-2024-03-13T10_01_20/cath-b-newest-latest-release

import argparse


def print_res(res):
    r = list(res.keys())
    for c in sorted(r):
        print(f"           {c*10}-{c*10+9} : {res[c]}")


# read arguments
parser = argparse.ArgumentParser(
    description='return statistics regarding the number of residues of domains in the given family')
parser.add_argument('--filename', action="store", help="cath database file name")
parser.add_argument('--family', action="store", help="family name")
args = parser.parse_args()
rozl = {}

with open(args.filename, "r") as fd:
    line = fd.readline()
    while line:
        sline = line.split(" ")
        if (sline[2] == args.family):
            rec = sline[3][:-1].split(",")
            dif = 0
            for r in rec:
                se = (r.split("-"))

                start = se[0]
                if start == "":
                    start = se[1]
                    see = se[2].split(":")
                else:
                    see = se[1].split(":")
                end = see[0]

                dif += float(end) - float(start)

            dif = int(dif / 10)

            if dif not in rozl:
                rozl[dif] = 1
            else:
                rozl[dif] += 1
        line = fd.readline()

print("number of residues : number of domains of such size")
print_res(rozl)
