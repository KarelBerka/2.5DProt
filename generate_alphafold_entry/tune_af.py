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
import utils_pdbe_api


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
parser.add_argument('--best', required=True, action="store",
                    help="generate database stats for alfa fold structures in the given directory")
parser.add_argument('--af', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--pdb', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--tag', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
res_file = f"stats-{args.tag}-{timestamp}"
stats = {}
ssum = {}
f_list = []
for filename in os.listdir(args.best):
    if not filename.startswith("dir"):
#        f_list.append(filename)
        print("----------------------------")
        with open(f"{args.best}/{filename}", "r") as fd:
            data = json.load(fd)
#            print(data.keys())
            keys = data.keys()
#            print(list(keys)[0])
            for i in data[list(keys)[0]]:
                protein = i['pdb_id']
                chain_id = i['chain_id']
                try:
                    struct_asym_id = utils_pdbe_api.transform_uniprot_chain_id_to_struct_asym_id(protein, list(keys)[0], chain_id)
                except:
                    print ("KKK -print ", protein)
                    continue
                
#                print(f"{protein}-{struct_asym_id}-{i['start']}-{i['end']}-detected.sses.json", args.pdb)
                name=f"{protein}-{struct_asym_id}-{i['start']}-{i['end']}-detected.sses.json"
                for j in os.listdir(args.pdb):
                    if j == name:
                       print("je to tam", j)
                       print("args.af", args.af)
                       for k in os.listdir(args.af):
                           if k == f"{list(keys)[0]}-detected.sses.json":
                                print("prislusne", k, i["unp_start"], i["unp_end"])
                                with open(f"{args.af}/{k}", "r") as fd:
                                     data = json.load(fd)
                                     pom = data[list(keys)[0]]["secondary_structure_elements"]
                                     sses= []
                                     for l in pom: 
                                         if (l["start"] > i["unp_start"] or l["end"]>i["unp_start"]+3)  and (l["end"] < i["unp_end"] or i["start"]+3<i["unp_end"]):
                                             sses.append(l)
                                             print("J", end="")
                                         else:
                                             print("N", end="")
                                     print("KKK")
                                if not os.path.isdir(f"{args.af[:-1]}_vyber"):
                                     os.mkdir(f"{args.af[:-1]}_vyber")
                                with open(f"{args.af[:-1]}_vyber/{list(keys)[0]}-{i['unp_start']}_{i['unp_end']}-detected.sses.json", "w") as fd:
                                     data_o = data
                                     data_o[list(keys)[0]]["secondary_structure_elements"]=sses
                                     json.dump(data, fd, indent=4, separators=(',', ': '))
        

#            data_list = data[args.tag[
#            ssum[filename] = 0
#            for i in data_list:
#                ssum[filename] += data_list[i]
#                if i not in stats:
#                    stats[i] = {filename: data_list[i]}
#                else:
#                    stats[i][filename] = data_list[i]

sys.exit()

with open(f"{args.stats_dir}/{res_file}", "w") as fd:
    fd.write(" ")
    for j in f_list:
        if len(j)<11:
            end = len(j)-5
        else:
            end = 6
        fd.write(f"{j[1:end]} ")
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
