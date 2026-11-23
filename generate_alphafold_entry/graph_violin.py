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
parser.add_argument('--filen', action="store", help="input json")
parser.add_argument('--dirn', action="store", help="dir containing input jsons")
parser.add_argument('--column', action="store",
                    help="which column should be displayed (sses, strands, helices)")
parser.add_argument('--split', action="store_true",
                    help="if the wanted graph should contain split violin graph (the pairs are detected by the letter afret _ sign \n" \
                         "e.g. zea_af.json and zea_pdb.json")
args = parser.parse_args()

p = sns.load_dataset('titanic')

# set which column should be displayed
if not args.column:
    column = "sses"
else:
    column = args.column

# read input jsons
data_list = {}
if (args.filen):
    # there is only one input json file
    with open(args.filen, 'r') as infile:
        data = json.load(infile)
        data_list["sses"] = data["sses_list"]
        data_list["strands"] = data["strands_list"]
        data_list["helices"] = data["helices_list"]
        print(data)
elif (args.dirn):
    # read jsons from the whole dir
    data_list = []
    for filename in os.listdir(args.dirn):
        fnsplit = filename.split("_")
        with open(f"{args.dirn}/{filename}", 'r') as infile:
            data = json.load(infile)
            if args.split:
                vg_split = True
                vg_split_column = fnsplit[1].split(".")[0]
            else:
                vg_split = False
                vg_split_column = ""
            add_records(data_list, data, args.column, fnsplit, vg_split_column)


df = pd.DataFrame(data_list)
print(df)

# order_tag = ["Arabidopsis", "Caenorhabditis", "Candida", "Danio", "Dictyostelium", "Drosophila", "Escherichia", "Glycine",
#              "Homo", "Methanocaldococcus", "Mus", "Oryza", "Rattus", "Saccharomyces", "Schizosaccharomyces", "Zea"]
order_tag = ["Arab", "Caen", "Cand", "Danio", "Dict", "Dros", "Esch", "Glyc",
             "Homo", "Meth", "Mus", "Oryza", "Ratt", "Sacc", "Schi", "Zea"]

if args.filen:
    sns.violinplot(x=df[column])
elif args.dirn:
    if args.split:
        sns.violinplot(data=df, y=f"{column} number", x="organism", hue="protein database", split=True, order=order_tag, cut=0, bw=.05)
#        sns.violinplot(data=df, y=f"{column} number", x="organism", hue="protein database", split=True, cut=0, bw=.03)
#        sns.boxplot(data=df, y=f"{column} number", x="organism", hue="protein database", split=True, order=order_tag, cut=0, bw=.03)
#         inner="stick)
    else:
        sns.violinplot(data=df, y=f"{column} number", x="organism", hue="protein database", bw=0.02)

plt.show()
