#!/usr/bin/env python3
# TODO: fix this tool which generate for given protein its directory with
# * protein images
# and image per chain


import sys
import argparse
import time
import os
import json
import collections
import shutil

import utils_common
import utils_annotation
import utils_pdbe_api

sys.path.append('..')
import utils_sses_types


def create_af_list_from_dir(af_dir):
    # from the input directory containing .cif files, find out the list of proteins
    af_list = []
    for filename in os.listdir(af_dir):
        if filename.endswith(".cif"):
            af_list.append(filename[3:-13])
    return af_list


def create_af_list_from_det_dir(af_dir):
    # from the input directory containing -detect files, find out the list of proteins
    af_list = []
    annotation_dir = af_dir
    for filename in os.listdir(annotation_dir):
        if filename.endswith("-F2-detected.sses.json"):
            prot_name = filename[:-22]
            counter = 2
            with open(f"{annotation_dir}/{prot_name}-F1-detected.sses.json", "r") as fd:
                data = json.load(fd)
            # read the first sses
            while os.path.isfile(f"{annotation_dir}/{prot_name}-F{counter}-detected.sses.json"):
                # add another one
                with open(f"{annotation_dir}/{prot_name}-F{counter}-detected.sses.json", "r") as fd:
                    data_additional = json.load(fd)
                for i in data_additional[f"{prot_name}-F{counter}"]["secondary_structure_elements"]:
                    if i["end"] > 1200:
                        i["end"] = i["end"] + (counter - 1) * 200
                        i["start"] = i["start"] + (counter - 1) * 200
                        data[f"{prot_name}-F1"]["secondary_structure_elements"].append(i)
                os.remove(f"{annotation_dir}/{prot_name}-F{counter}-detected.sses.json")
                counter += 1

            # fflush it to F1 one
            with open(f"{annotation_dir}/{prot_name}-F1-detected.sses.json", "w") as fd:
                json.dump(data, fd, indent=4, separators=(',', ':'))

    for filename in os.listdir(af_dir):
        if filename.endswith("-detected.sses.json"):
            af_list.append(filename[:-19])
        if filename.endswith("-annotated.sses.json"):
            af_list.append(filename[:-20])
    return af_list


def create_pdb_list_from_dir(pdb_dir):
    # from the input pdb best matches files dir download cifs and find out the list of proteins
    pdb_list = []
    for filename in os.listdir(pdb_dir):
        with open(f"{pdb_dir}/{filename}", "r") as fd:
            data = json.load(fd)
            pcd = -1
            for j in data.keys():
                fkey = j
            for i in data[fkey]:
                if i["coverage"] > 0.8:
                    protein = i["pdb_id"]
                    chain_id = i["chain_id"]
                    domain = [i["start"], i["end"]]
                    try:
                        struct_asym_id = utils_pdbe_api.transform_uniprot_chain_id_to_struct_asym_id(protein, fkey, chain_id)
                    except:
                        struct_asym_id = chain_id
                    pcd = [protein, struct_asym_id, domain, filename, chain_id]
                break
        if pcd != -1:
            pdb_list.append(pcd)
        else:
            print(f"for {filename} only pcd with low coverage")
    return pdb_list


def add_to_dict(input_dict, key, key_value):
    if key not in input_dict:
        input_dict[key] = key_value
    else:
        input_dict[key] += key_value


def add_plain_stats(filename, fdir, stats, args):
    # add statistics item from detected-sses.json files
    strands = 0
    helices = 0
    with open(f"{fdir}/{filename}", "r") as fd:
        data = json.load(fd)
        dkey = list(data.keys())[0]
        sses_list = data[dkey]["secondary_structure_elements"]
        for sses in sses_list:
            sses_res_len = sses['end'] - sses['start']
            if type(sses['start_vector'][0]) == "NaN":
                # problem with A0A0K0ENM3-F1 annotation
                break
            if utils_sses_types.is_helic(sses["type"]) and utils_sses_types.is_nontrivial_sse(sses["label"], sses_res_len):
                helices += 1
                if sses_res_len > 200:
                    sf = filename.split("-")
                    utils_common.curl_exists_file(f"https://alphafold.ebi.ac.uk/files/AF-{sf[0]}-{sf[1]}-model_v4.cif")
                stats["helices_len_list"].append(sses_res_len)
                add_to_dict(stats["helices_len"], sses_res_len, 1)
            if utils_sses_types.is_sheet(sses["type"]) and utils_sses_types.is_nontrivial_sse(sses["label"], sses_res_len):
                strands += 1
                stats["strands_len_list"].append(sses_res_len)
                add_to_dict(stats["strands_len"], sses_res_len, 1)

        add_to_dict(stats["helices"], helices, 1)
        add_to_dict(stats["strands"], strands, 1)
        add_to_dict(stats["sses"], strands + helices, 1)

        stats["helices_list"].append(helices)
        stats["strands_list"].append(strands)
        sf = filename.split("-")
#        name = f"{sf[0]}-{sf[1]}.cif"
#        if (strands == 2) and (helices == 10):
#           print(f"dvadesitka {sf[0]}")
#        if strands+helices >90:
#            print( helices, strands, filename)
        stats["sses_list"].append(strands + helices)
        return True
    return False


def average(val_list):
    val_sum = 0
    val_number = 0
    for i in val_list:
        val_sum += i * val_list[i]
        val_number += val_list[i]
    if val_number == 0:
        return 0
    else:
        return val_sum / val_number


def median(val_list):
    sorted_val_list = collections.OrderedDict(sorted(val_list.items()))
    val_sum = 0
    res = 0
    for i in val_list:
        val_sum += val_list[i]
    val_number = 0
    for i in sorted_val_list:
        if val_number*2 > val_sum:
            res = i
            break
        val_number += sorted_val_list[i]
    return res


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--af_dir', action="store",
                    help="generate database stats for alfa fold .cif structures in the given directory")
parser.add_argument('--d_dir', action="store",
                    help="generate database stats for alfa fold structures in the given detected directory")
parser.add_argument('--input_cif', action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--name', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--pdb_dir', action="store",
                    help="generate database stats for pdb https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/ input files dir")

args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
name = args.name

new_dir = f"stats-{name}-without_{timestamp}"
os.makedirs(new_dir, exist_ok=True)

if args.d_dir:
    af_list = create_af_list_from_det_dir(args.d_dir)

    # Annotate protein structure
    annotation_dir = f"{new_dir}/annotation"
    os.mkdir(annotation_dir)
    for af in af_list:
        if os.path.isfile(f"{args.d_dir}/{af}-detected.sses.json"):
            os.symlink(f"../../{args.d_dir}/{af}-detected.sses.json", f"{annotation_dir}/{af}-detected.sses.json")
        if os.path.isfile(f"{args.d_dir}/{af}-annotated.sses.json"):
            os.symlink(f"../../{args.d_dir}/{af}-annotated.sses.json", f"{annotation_dir}/{af}-detected.sses.json")

elif args.af_dir:
    af_list = create_af_list_from_dir(args.af_dir)

    # we want to download read relevant cif files
    cif_dir = f"{new_dir}/cif"
    os.mkdir(cif_dir)
    # simlink is better - not to duplicate the data
    utils_common.symlink_alphafold_list(cif_dir, args.af_dir, af_list)

    annotation_dir = f"{new_dir}/annotation"
    os.mkdir(annotation_dir)
    ret = utils_annotation.detect_list(annotation_dir, cif_dir, af_list, new_dir, whole_protein=True)
    if (ret == 0):
        sys.exit()

    for filename in os.listdir(annotation_dir):
        if filename.endswith("-F2-detected.sses.json"):
            prot_name = filename[:-22]
            counter = 2
            with open(f"{annotation_dir}/{prot_name}-F1-detected.sses.json", "r") as fd:
                data = json.load(fd)
            # read the first sses
            while os.path.isfile(f"{annotation_dir}/{prot_name}-F{counter}-detected.sses.json"):
                # add another one
                with open(f"{annotation_dir}/{prot_name}-F{counter}-detected.sses.json", "r") as fd:
                    data_additional = json.load(fd)
                for i in data_additional[f"{prot_name}-F{counter}"]["secondary_structure_elements"]:
                    if i["end"] > 1200:
                        i["end"] = i["end"] + (counter - 1) * 200
                        i["start"] = i["start"] + (counter - 1) * 200
                        data[f"{prot_name}-F1"]["secondary_structure_elements"].append(i)
                os.remove(f"{annotation_dir}/{prot_name}-F{counter}-detected.sses.json")
                counter += 1

            # fflush it to F1 one
            with open(f"{annotation_dir}/{prot_name}-F1-detected.sses.json", "w") as fd:
                json.dump(data, fd, indent=4, separators=(',', ':'))

else:
    pdb_list = create_pdb_list_from_dir(args.pdb_dir)
    cif_dir = f"{new_dir}/cif"
    os.mkdir(cif_dir)
    for i in pdb_list:
        utils_pdbe_api.download_protein_cif(cif_dir, i[0])

    annotation_dir = f"{new_dir}/annotation"
    os.mkdir(annotation_dir)
    number = 0
    for pdb in pdb_list:
        number += 1
        ret = utils_annotation.detect_list(annotation_dir, cif_dir, [pdb], new_dir, whole_protein="pdbe_best_structures")
        if (ret == 0):
            print("bug {pdb}")
        else:
            shutil.move(f"{annotation_dir}/{pdb[0]}-detected.sses.json", f"{annotation_dir}/{number}-{pdb[0]}-detected.sses.json")

# now we have to go through all -detected- files and count stats
stats = {}
stats["sses"] = {}
stats["sses_list"] = []
stats["helices"] = {}
stats["helices_list"] = []
stats["strands"] = {}
stats["strands_list"] = []
stats["helices_len"] = {}
stats["helices_len_list"] = []
stats["strands_len"] = {}
stats["strands_len_list"] = []
counter = 0
for filename in os.listdir(annotation_dir):
    # -detected files are in annotation_dir
    if filename.endswith("-detected.sses.json"):
        counter += 1
        ret = add_plain_stats(filename, annotation_dir, stats, args)

        if ret is False:
            print(f"{filename} can't be parsed")
            sys.exit()

stats["avg_helices_number"] = average(stats["helices"])
stats["avg_strands_number"] = average(stats["strands"])
stats["avg_ses_number"] = average(stats["sses"])
stats["avg_helices_len"] = average(stats["helices_len"])
stats["avg_strands_len"] = average(stats["strands_len"])
# print("med_helices_number")
# stats["med_helices_number"] = median(stats["helices"])
# print("med_strands_number")
# stats["med_strands_number"] = median(stats["strands"])
# print("med_ses_number")
# stats["med_ses_number"] = median(stats["sses"])
# stats["med_helices_len"] = median(stats["helices_len"])
# stats["med_strands_len"] = median(stats["strands_len"])

stats["proteins"] = counter

with open(f"{new_dir}/{name}_without.json", "w") as fd:
    json.dump(stats, fd, indent=4, separators=(',', ':'))
