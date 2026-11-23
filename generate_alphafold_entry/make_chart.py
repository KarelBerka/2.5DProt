#!/usr/bin/env python3
# tool which generate for given alpha fold stats subdirectory
# summary stats record for all its subdirectory
import sys
import argparse
import time
import os
import json
import math

sys.path.append('..')
import utils
import utils_sses_types


def add_to_dict(input_dict, key, key_value):
    if key not in input_dict:
        input_dict[key] = key_value
    else:
        input_dict[key] += key_value


def add_stats(filename, fdir, stats):
    # add statistics item from detected-sses.json files
    name = filename[:-19]
    strands = 0
    helices = 0
    with open(f"{fdir}/{filename}", "r") as fd:
        data = json.load(fd)
        sses_list = data[name]["secondary_structure_elements"]
        for sses in sses_list:
            sses_len = int(math.dist(sses['start_vector'], sses['end_vector']))
            if utils_sses_types.is_helic(sses["type"]):
                helices += 1
                add_to_dict(stats["helices_length"], sses_len, 1)
            if utils_sses_types.is_sheet(sses["type"]):
                strands += 1
                add_to_dict(stats["strands_length"], sses_len, 1)

        add_to_dict(stats["helices"], helices, 1)
        add_to_dict(stats["strands"], strands, 1)
        add_to_dict(stats["sses"], strands + helices, 1)
        return True

    return False


def average(val_list):
    val_sum = 0
    val_number = 0
    for i in val_list:
        val_sum += i * val_list[i]
        val_number += val_list[i]
    return val_sum / val_number


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--stats_dir', required=True, action="store",
                    help="generate database stats for alfa fold structures in the given directory")
parser.add_argument('--tag', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
res_file = f"stats-{args.tag}-{timestamp}"
stats = {}
ssum = {}
f_list = []
phn = 0
pstr = 0
pcount = 0
for filename in os.listdir(args.stats_dir):
    if filename.endswith(".json"):
        f_list.append(filename)
        with open(f"{args.stats_dir}/{filename}", "r") as fd:
            data = json.load(fd)
            print(filename, data['avg_helices_number'])
            phn += data['avg_helices_number']
            pstr += data['avg_strands_number']
            pcount += 1
            data_list = data[args.tag]
            ssum[filename] = 0
            for i in data_list:
                ssum[filename] += data_list[i]
                if i not in stats:
                    stats[i] = {filename: data_list[i]}
                else:
                    stats[i][filename] = data_list[i]
print("avg helices number", phn / pcount)
print("avg strands number", pstr / pcount)

with open(f"{args.stats_dir}/{res_file}", "w") as fd:
    fd.write(" ")
    for j in f_list:
        if len(j)<6:
            end = len(j)
        else:
            end = 5
        fd.write(f"{j[0:end]} ")
    fd.write("\n")
    for i in stats:
        fd.write(f"{i} ")
        for j in f_list:
            if j in stats[i]:
                fd.write(f"{stats[i][j]/ssum[j]*100} ")
            else:
                fd.write("0 ")
        fd.write("\n")

print(f"{args.stats_dir}/{res_file}")
