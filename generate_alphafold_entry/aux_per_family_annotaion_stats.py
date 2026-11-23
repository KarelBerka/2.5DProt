#!/usr/bin/env python3
# EXAMPLE
# ./aux_per_family_annotaion_stats.py --dir generated-3.20.20.700-2024-04-21T18_29_17/annotation/
# 10000 annotation files processed
# 20000 annotation files processed
# 30000 annotation files processed
# average number of sses in domain is 25.100833902939165
# sses which have > 0.7 ration:
# E51 0.8139439507860561
# H77 0.933151059466849
# E84 0.9448257006151743
# H102 0.947423103212577
# E107 0.9739986329460014
# H117 0.9321394395078606
# E119 0.9837047163362953
# H129 0.9326315789473684
# E134 0.9639371155160629
# H147 0.9445249487354751
# E152 0.9351469583048531
# H172 0.8918386876281613
# E177 0.7685850991114149
# average number of SSEs in domain is                      - 25.100833902939165
# average ratio for SSEs (lower - worse annotation) is     - 0.07070657437447656

import argparse
import json
import os


# read arguments
parser = argparse.ArgumentParser(
    description='show statistics regarding annotation quality for fiven family annotation dir')
parser.add_argument('--dir', action="store", help="dir containing input annoation jsons")
args = parser.parse_args()


stats = {}
sses_number_per_domain = {}
counter = 0
sses_sum = 0
for f in os.listdir(args.dir):
    if f.endswith("-annotated.sses.json"):

        counter += 1
        sses_n = 0

        with open(f"{args.dir}/{f}", "r") as fd:
            data = json.load(fd)
            k = list(data.keys())
            cc = 0
            for sse in data[k[0]]["secondary_structure_elements"]:
                cc += 1
                if sse['label'] not in stats:
                    stats[sse['label']] = [f]
                else:
                    stats[sse['label']].append(f)

            if cc not in sses_number_per_domain:
                sses_number_per_domain[cc] = 1
            else:
                sses_number_per_domain[cc] += 1
            sses_sum += cc

        if counter % 10000 == 0:
            print(f"{counter} annotation files processed")

print(f"average number of sses in domain is {sses_sum / counter}")

print("sses which have > 0.7 ration:")
rsum = 0
rcounter = 0
for i in stats:
    ratio = float(len(stats[i]))/counter
    rsum += ratio
    rcounter += 1
    if ratio > 0.70:
        print(i, float(len(stats[i]))/counter)

print(f"average number of SSEs in domain is                      - {sses_sum / counter}")
print(f"average ratio for SSEs (lower - worse annotation) is     - {rsum / rcounter}")
