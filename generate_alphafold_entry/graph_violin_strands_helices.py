#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import json
import os


def count_rmsd(data):
    #   TODO:
    aux = order_tag

    for org in aux:
        diffh = 0
        diffs = 0
        for i in range(10, 20):
            if str(i) not in data["helices"][org]:
                val_org1 = 0
            else:
                val_org1 = data["helices"][org][str(i)]
            diffh += val_org1
        for i in range(0, 5):
            if str(i) not in data["strands"][org]:
                val_org1 = 0
            else:
                val_org1 = data["strands"][org][str(i)]
            diffs += val_org1

    for org1 in aux:
        print(org1, end=", ")
        for org2 in aux:
            diff = 0
            for i in range(51):
                if str(i) not in data["helices"][org1]:
                    val_org1 = 0
                else:
                    val_org1 = data["helices"][org1][str(i)]
                if str(i) not in data["helices"][org2]:
                    val_org2 = 0
                else:
                    val_org2 = data["helices"][org2][str(i)]
                diff += (val_org1 - val_org2) * (val_org1 - val_org2)

            print(f"--- {org2}, {int(diff)}", end=", ")
        print()


def add_records(data_list, rmsd_data_list, name_list, input_data, organism_name, all_o):
    # add a records to data_list structure
    # from input data structure
    # column is column name and added data
    # organism_name is the name of organism which will be displayed on diagram
    # if the s[lit plot is wanted use vg_split and vg_split_column
    if len(organism_name) > 5:
        org_name = f"{organism_name[:5]}."
    else:
        if len(organism_name) == 4:
            org_name = f"{organism_name[0][:2]}_{organism_name[1][:2]}"
        else:
            org_name = organism_name[0][:4]
    if all_o:
        org_name = all_o
    spliton = organism_name[0].split(" ")
    if len(spliton) == 1:
        org_name = f"{spliton[0][:4]}. "
    elif len(spliton[1]) < 8:
        org_name = f"{spliton[0][:4]}. {spliton[1]}"
    else:
        org_name = f"{spliton[0][:4]}. {spliton[1][:7]}."
    name_list.add(org_name)

    for sses_num in data["strands_list"]:
        add = {"type": "sses", "sses number": sses_num, "organism": org_name, "sses type": "strands"}
        if sses_num < 51:
            data_list.append(add)
    for sses_num in data["helices_list"]:
        add = {"type": "sses", "sses number": sses_num, "organism": org_name, "sses type": "helices"}
        if sses_num < 51:
            data_list.append(add)
    for sses_num in data["helices"]:
        if int(sses_num) < 51:
            if org_name not in rmsd_data_list["helices"]:
                rmsd_data_list["helices"][org_name] = {}
            rmsd_data_list["helices"][org_name][sses_num] = data["helices"][sses_num]/(data["proteins"])*100
    for sses_num in data["strands"]:
        if int(sses_num) < 51:
            if org_name not in rmsd_data_list["strands"]:
                rmsd_data_list["strands"][org_name] = {}
            rmsd_data_list["strands"][org_name][sses_num] = data["strands"][sses_num]/(data["proteins"])*100
    for sses_num in data["sses"]:
        if int(sses_num) < 51:
            if org_name not in rmsd_data_list["sses"]:
                rmsd_data_list["sses"][org_name] = {}
            rmsd_data_list["sses"][org_name][sses_num] = data["sses"][sses_num]/(data["proteins"])*100
    for sses_num in data["helices_len"]:
        if int(sses_num) < 51:
            if org_name not in rmsd_data_list["len_helices"]:
                rmsd_data_list["len_helices"][org_name] = {}
            rmsd_data_list["len_helices"][org_name][sses_num] = data["helices_len"][sses_num]/(data["proteins"]*data["avg_helices_number"])*100
    for sses_num in data["strands_len"]:
        if int(sses_num) < 51:
            if org_name not in rmsd_data_list["len_strands"]:
                rmsd_data_list["len_strands"][org_name] = {}
            rmsd_data_list["len_strands"][org_name][sses_num] = data["strands_len"][sses_num]/(data["proteins"]*data["avg_strands_number"])*100


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--dirn', required=True, action="store", help="dir containing input jsons")
parser.add_argument('--all', action="store", help="dir containing input jsons")
args = parser.parse_args()

# read jsons from the whole dir
data_list = []
rmsd_data_list = {"helices": {}, "strands": {}, "sses": {}, "len_helices": {}, "len_strands": {}, "avg_helices_number": {}, "proteins": {}}
diff_data_list = {}
prot_number = 0
name_list = set()
for filename in os.listdir(args.dirn):
    fnsplit = filename.split("_")
    with open(f"{args.dirn}/{filename}", 'r') as infile:
        data = json.load(infile)
        prot_number += 1
        add_records(data_list, rmsd_data_list, name_list, data, fnsplit, args.all)

rmsd_data_list['helices']['all'] = {}
rmsd_data_list['strands']['all'] = {}
for org in name_list:
    for c in rmsd_data_list['helices'][org]:
        if c in rmsd_data_list['helices']['all']:
            rmsd_data_list['helices']['all'][c] += rmsd_data_list['helices'][org][c]
        else:
            rmsd_data_list['helices']['all'][c] = rmsd_data_list['helices'][org][c]

    rmsd_data_list['helices']['all'][c] = rmsd_data_list['helices']['all'][c]/prot_number
    for c in rmsd_data_list['strands'][org]:
        if c in rmsd_data_list['strands']['all']:
            rmsd_data_list['strands']['all'][c] += rmsd_data_list['strands'][org][c]
        else:
            rmsd_data_list['strands']['all'][c] = rmsd_data_list['strands'][org][c]

    rmsd_data_list['strands']['all'][c] = rmsd_data_list['strands']['all'][c]/prot_number

for c in rmsd_data_list['helices']['all']:
        rmsd_data_list['helices']['all'][c] = rmsd_data_list['helices']['all'][c]/prot_number

for c in rmsd_data_list['strands']['all']:
        rmsd_data_list['strands']['all'][c] = rmsd_data_list['strands']['all'][c]/prot_number

hdiffsum = 0
sdiffsum = 0

for org in name_list:
    hdiff = 0
    sdiff = 0
    for c in rmsd_data_list['helices']['all']:
        if c in rmsd_data_list['helices'][org]:
            hdiff += abs(rmsd_data_list['helices']['all'][c] - rmsd_data_list['helices'][org][c])
        else:
            hdiff += abs(rmsd_data_list['helices']['all'][c])

    for c in rmsd_data_list['strands']['all']:
        if c in rmsd_data_list['strands'][org]:
            sdiff += abs(rmsd_data_list['strands']['all'][c] - rmsd_data_list['strands'][org][c])
        else:
            sdiff += abs(rmsd_data_list['strands']['all'][c])

    hdiffsum += hdiff
    sdiffsum += sdiff

print(f"avgdiff helices {hdiffsum/len(name_list)} strands {sdiffsum/len(name_list)}")

name_list = list(name_list)

print("------------")
if not args.all:
    print("str ", end=" ")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
    print()

    for org1 in name_list:
        print(f"{org1:04s}", end=" ")

        for org2 in name_list:
            hdiff = 0
            for c in rmsd_data_list['helices']['all']:
                if c in rmsd_data_list['helices'][org1]:
                    if c in rmsd_data_list['helices'][org2]:
                        hdiff += abs(rmsd_data_list['helices'][org1][c] - rmsd_data_list['helices'][org2][c])
                    else:
                        hdiff += abs(rmsd_data_list['helices'][org1][c])
                elif c in rmsd_data_list['helices'][org2]:
                    hdiff += abs(rmsd_data_list['helices'][org2][c])
            sdiff = 0
            for c in rmsd_data_list['strands']['all']:
                if c in rmsd_data_list['strands'][org1]:
                    if c in rmsd_data_list['strands'][org2]:
                        sdiff += abs(rmsd_data_list['strands'][org1][c] - rmsd_data_list['strands'][org2][c])
                    else:
                        sdiff += abs(rmsd_data_list['strands'][org1][c])
                elif c in rmsd_data_list['strands'][org2]:
                    sdiff += abs(rmsd_data_list['strands'][org2][c])
            sdiffsum += sdiff
            print(f"{int(sdiff):02d}  ", end=" ")
        print()

print("============")

if not args.all:
    print("hel ", end=" ")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
    print()

    for org1 in name_list:
        print(f"{org1:04s}", end=" ")

        for org2 in name_list:
            hdiff = 0
            for c in rmsd_data_list['helices']['all']:
                if c in rmsd_data_list['helices'][org1]:
                    if c in rmsd_data_list['helices'][org2]:
                        hdiff += abs(rmsd_data_list['helices'][org1][c] - rmsd_data_list['helices'][org2][c])
                    else:
                        hdiff += abs(rmsd_data_list['helices'][org1][c])
                elif c in rmsd_data_list['helices'][org2]:
                    hdiff += abs(rmsd_data_list['helices'][org2][c])
            sdiff = 0
            for c in rmsd_data_list['strands']['all']:
                if c in rmsd_data_list['strands'][org1]:
                    if c in rmsd_data_list['strands'][org2]:
                        sdiff += abs(rmsd_data_list['strands'][org1][c] - rmsd_data_list['strands'][org2][c])
                    else:
                        sdiff += abs(rmsd_data_list['strands'][org1][c])
                elif c in rmsd_data_list['strands'][org2]:
                    sdiff += abs(rmsd_data_list['strands'][org2][c])
            sdiffsum += sdiff
            print(f"{int(hdiff):02d}  ", end=" ")
        print()
print("============")

if not args.all:
    print("helices")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
        c = 0
        while c < 50:
            if str(c) in rmsd_data_list["helices"][org1]:
                print(rmsd_data_list["helices"][org1][str(c)], end=" ")
            else:
                print("0", end=" ")
            c += 1
        print()
    print()

if not args.all:
    print("strands")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
        c = 0
        while c < 50:
            if str(c) in rmsd_data_list["strands"][org1]:
                print(rmsd_data_list["strands"][org1][str(c)], end=" ")
            else:
                print("0", end=" ")
            c += 1
        print()
    print()

if not args.all:
    print("sses")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
        c = 0
        while c < 50:
            if str(c) in rmsd_data_list["sses"][org1]:
                print(rmsd_data_list["sses"][org1][str(c)], end=" ")
            else:
                print("0", end=" ")
            c += 1
        print()
    print()

if not args.all:
    print("len_helices")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
        c = 0
        while c < 50:
            if str(c) in rmsd_data_list["len_helices"][org1]:
                print(rmsd_data_list["len_helices"][org1][str(c)], end=" ")
            else:
                print("0", end=" ")
            c += 1
        print()
    print()

if not args.all:
    print("len_strands")

    for org1 in name_list:
        print(f"{org1:0s}", end=" ")
        c = 0
        while c < 50:
            if str(c) in rmsd_data_list["len_strands"][org1]:
                print(rmsd_data_list["len_strands"][org1][str(c)], end=" ")
            else:
                print("0", end=" ")
            c += 1
        print()
    print()


df = pd.DataFrame(data_list)

order_tag = ["Arabidopsis", "Caenorhabditis", "Candida", "Danio", "Dictyostelium", "Drosophila", "Escherichia", "Glycine",
             "Homo", "Methanocaldococcus", "Mus", "Oryza", "Rattus", "Saccharomyces", "Schizosaccharomyces", "Zea"]
order_tag = ["Arabido..", "Caenorh..", "Candida", "Danio", "Dictyos..", "Drosoph..", "Escheri..", "Glycine",
             "Homo", "Methano..", "Mus", "Oryza", "Rattus", "Sacchar..", "Schizos..", "Zea"]
order_tag = ["Caen", "Dani", "Mus", "Ratt", "Homo", "Dros", "Brug", "Wuch", "Onch", "Stro", "Tric", "Drac", "Schi", "Camp", "Ente",
             "Stap", "Stre", "Haem", "Kleb", "Pseu", "Salm", "Shig", "Neis", "Heli",
             "My_le", "My_tu", "My_ul", "Noca", "Madu", "Spor", "Clad", "Ajel", "Para", "Fons",
             "Leis", "Tr_br", "Tr_cr", "Plas"]
# order_tag = [ "Brug", "Wuch", "Onch", "Stro", "Tric", "Drac", "Schi", "Camp", "Ente", \
#             "Stap", "Stre", "Haem", "Kleb", "Pseu", "Salm", "Shig", "Neis", "Heli",  \
#             "My_le", "My_tu", "My_ul", "Noca", "Madu", "Spor", "Clad", "Ajel", "Para", "Fons", \
#             "Leis", "Tr_br", "Tr_cr", "Plas" ]
# order_tag = ["Dani", "Mus", "Ratt", "Homo", "Dros", "Brug", "Wuch", "Onch", "Stro", "Tric", "Drac", "Caen", "Schi", ]
# order_tag = ["Dani", "Mus", "Ratt", "Homo", "Dros", "Brug", "Onch", "Tric", "Caen" ]

# order_tag = ["Madu", "Spor", "Clad", "Ajel", "Para", "Fons", "Sacc", "Cand", "Schi"]
# order_tag = ["Arab", "Glyc", "Oryz", "Zea"]
# order_tag = ["Camp", "Ente", "Haem", "Heli", "Kleb", "My_le", "My_tu", "My_ul", "Neis", \
#             "Noca", "Pseu", "Salm", "Shig", "Stap", "Stre", "Esch"]
# order_tag = ["Camp", "Ente", "Haem", "Heli", "Kleb", "My_le", "My_tu", "Neis", \
#             "Noca", "Pseu", "Salm", "Shig", "Stap", "Stre", "Esch"]
# order_tag = ['Ence', 'Nema', 'Spiz', '2Spi', 'Sync',  'Lobo', 'Pseu', 'Micr', 'Schi', 'Dani', 'Sacc']
# order_tag = ["Capr", "Odoc", "Biso", "Bos", "Feli", "Mus", "Ratt", "Dani", "Pogo", "Pipr", "Ayth", "Stro",  "Acti", "Dros", "Apis", "Varr", "Tric", "Caen", "Brug", "Capi", "Ling"]
print("----------------------------------------------------------------")


if args.all:
    print("aa", len(args.all))
    sns.violinplot(data=df, y="sses number", x="organism", hue="sses type", split=True, cut=0, bw=0.06, gridsize=200, scale="width", scale_hue=True)
else:
    sns.violinplot(data=df, y="sses number", x="organism", hue="sses type", split=True, cut=0, bw=0.08, gridsize=200, scale="width", scale_hue=True)

print("len(data)", len(data))
h_counter = 0
h_sum = 0
h_odch = 0
s_counter = 0
s_sum = 0
s_odch = 0
for i in data_list:
    if i["sses type"] == "helices":
        h_sum += i["sses number"]
        h_counter += 1
    if i["sses type"] == "strands":
        s_sum += i["sses number"]
        s_counter += 1

for i in data_list:
    if i["sses type"] == "helices":
        h_odch += (abs(i["sses number"]-h_sum / h_counter))**1
    if i["sses type"] == "strands":
        s_odch += (abs(i["sses number"]-s_sum / s_counter))**1

print("h - avg", h_sum / h_counter, h_odch / h_counter)
print("s - avg", s_sum / s_counter, s_odch / s_counter)
print("len(data_list)", len(data_list))
plt.show()
