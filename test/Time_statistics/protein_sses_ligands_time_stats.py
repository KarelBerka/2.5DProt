#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse
import sys
import os
import os.path
import json


def process_protein_generated_dir(directory, min_time):

    prot = {"name": "", "sses": "", "ligands": "", "time": ""}
    for filename in os.listdir(directory):
        fn = f"{directory}/{filename}"
        if filename.endswith(".json") and not filename.startswith("chains"):
            with open(fn, "r") as fd:
                data = json.load(fd)
            if "sses" in data:
                prot["sses"] = len(data["sses"])
            else:
                prot["sses"] = 0
            if "ligands" in data:
                prot["ligands"] = len(data["ligands"])
            else:
                prot["ligands"] = 0
            prot["name"] = filename[0:4]

        if os.path.isfile(fn) and filename.startswith("main"):
            with open(fn, "r") as fd:
                line = fd.readline()
                while line:
                    if "Timestamp" in line:
                        val = line[:-5].rfind(" ")
                        c = int(line[15])
                        if c == 8:
                            prot["name"] = filename[0:4]
                            prot["time"] = float(line[val:-5])

                    line = fd.readline()
    if prot["time"] > min_time:
        print(prot)


# read arguments
parser = argparse.ArgumentParser(
description="print the list of per protein records containing \n"
            "protein name, number of sses, number of ligands, computational time\n"
            "like:\n"
            "{'name': '8si4', 'sses': 248, 'ligands': 65, 'time': 1.8543624043464662}\n"
            "{'name': '8si5', 'sses': 252, 'ligands': 65, 'time': 2.068746340274811}\n"
            "{'name': '8si6', 'sses': 252, 'ligands': 65, 'time': 2.1019378662109376}\n"
            "..")
parser.add_argument('--dir', action="store", required=True,
                    help="path where computed protein diagram dirs are stored (in multilevel hierarchy)")
parser.add_argument('--min', action="store",
                    help="the records with at least this time are returned (in seconds)")

args = parser.parse_args()


# Load params - input family generation dir
if len(sys.argv) < 2:
    print("ERROR: Argument missing")
    print("Usage:")
    print("    ./protein_time_stats protein enerated dbs dir")
    sys.exit()

input_dir = args.dir

stat = {}
label = {}
single_images = 0
protein_sizes = []

if args.min:
    min_time = int(args.min)
else:
    min_time = 0

for dirname1 in os.listdir(input_dir):
    for dirname2 in os.listdir(f"{input_dir}/{dirname1}"):
        for dirname3 in os.listdir(f"{input_dir}/{dirname1}/{dirname2}"):
            if dirname3.startswith("generated-"):
                # protein diagram generation directory
                process_protein_generated_dir(f"{input_dir}/{dirname1}/{dirname2}/{dirname3}", min_time)

