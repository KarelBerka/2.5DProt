#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import json
import os


def add_records(data_list, input_data, column, organism_name, vg_split_column=""):
    # add a records to data_list structure
    # from input data structure
    # column is column name and added data
    # organism_name is the name of organism which will be displayed on diagram
    # if the s[lit plot is wanted use vg_split and vg_split_column
    if vg_split_column == "pdf":
        vg_split_column = "pdbe"
    elif vg_split_column == "pdb":
        vg_split_column = "pdbe"
    elif vg_split_column == "af":
        vg_split_column = "AlphaFold"

    if len(organism_name[0]) > 5:
        if len(organism_name) == 4:
            org_name = f"{organism_name[0][:2]}_{organism_name[1][:2]}"
        else:
            org_name = f"{organism_name[0][:4]}"
    else:
        org_name = organism_name[0]
    for sses_num in data[f"{column}_list"]:
        add = {
            "type": "sses",
            f"{column} number": sses_num,
            "organism": org_name,
            "protein database": vg_split_column, }
        data_list.append(add)


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--dir', action="store", required="True", help="input json")
args = parser.parse_args()

counter = 0
for f in os.listdir(args.dir):
     if f.endswith("json"):
         counter += 1

print("counter", counter)