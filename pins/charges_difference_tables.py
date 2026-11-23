#!/usr/bin/env python3
# tool which generate for given pin directory its standart graphs with up/down part charges sum

import argparse
import os
import json


# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
args = parser.parse_args()

dir_list = []  # list of all pin subdirectories (PIN1,PIN2,...)
data = {}      # charges data data[sse][dirname][filename] = list of res charges - from 0 to end of "sse"
counter = 0

for dirname in os.listdir(args.directory):
  if os.path.isdir(dirname):
      continue
  dir_list.append(dirname)
  # go through all PINs in PIN directory and create data record for them
  for filename in os.listdir(f"{args.directory}/{dirname}"):
    # read PIN json and read residues charges in SSEs
    with open(f"{args.directory}/{dirname}/{filename}", "r") as fd:
        counter += 1
        jdata = json.load(fd)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in data:
                    data[sse] = {}
                if dirname not in data[sse]:
                    data[sse][dirname] = {}
                res_min = 100000
                for i in jdata["charges"][sse]:
                    if int(i) < res_min:
                       res_min = int(i)
                rec = {}
                for i in jdata["charges"][sse]:
                    rec[int(i) - res_min] = jdata["charges"][sse][i] / 100
                data[sse][dirname][filename] = rec
# data are red
print(f"data read ({counter} files), counting avg for each pin type ================================")
print(f"data.keys {data.keys()}")
print(f"data.keys {data['H6'].keys()}")
for pi in data['H6'].keys():
    print("---", pi, len(data['H6'][pi]))

counter = {}
avg = {}  # avg charges per sses/pin class
for sse in data:
    avg[sse] = {}
    counter[sse] = {}
    for dirname in data[sse]:
        avg[sse][dirname] = []
        counter[sse][dirname] = []
        is_in = True
        res_counter = 0
        while is_in:
            is_in = False
            avg_ch = 0
            dom_counter = 0
            for filename in data[sse][dirname]:
                 if len(data[sse][dirname][filename]) > res_counter:
                      dom_counter += 1
                      avg_ch += data[sse][dirname][filename][res_counter]
                      is_in = True
            if is_in is True:
                 avg[sse][dirname].append(avg_ch/dom_counter)
                 counter[sse][dirname].append(dom_counter)
            res_counter += 1
print("avg", avg["H59"]["PIN1"])
print("avg", counter["H59"]["PIN1"])
print(len(counter["H59"]["PIN1"]))

sses_list = list(avg.keys())

if list(avg.keys()) == ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"]:  # small dataset
     sses_dict = {"H2": "1", "H3": "2", "H5": "3", "H6": "4a", "H8": "4b", "H9": "5", "H19": "6a", "H20": "6b", "H21": "7", "H22": "8", "H23": "9a", "H25": "9b", "H26": "10"}  # small dataset
elif list(avg.keys()) == ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]:  # BIG DATASET
     sses_dict = {"H6": "1", "H7": "2", "H9": "3", "H11": "4a", "H13": "4b", "H15": "5", "H44": "6a", "H45": "6b", "H48": "7", "H51": "8", "H53": "9a", "H56": "9b", "H59": "10"}  # BIG DATASET
else:
     print("WARNING: there is no mapping to cannoninal pin names for this dataset")
     sses_dict = {}
     for i in list(avg.keys()):
         sses_dict[i] = i

if True:
    # average charges values per SSE per PIN
    print("1/ average charges values per SSE per PIN ================================")
    for sse in avg:
        if sse in sses_list:
            print(f"sse: {sses_dict[sse]}")
            for dirname in dir_list:
                print(f"avg {dirname}:", end=" ")
                for l in range(len(avg[sse][dirname])):
                    if counter[sse][dirname][l] == counter[sse][dirname][0]:
                        if avg[sse][dirname][l] < 0.1 and avg[sse][dirname][l] > -0.1:
                            print("0", end=" ")
                        else:
                            print(f"{avg[sse][dirname][l]:.1f}", end=" ")
                    else:
                        if avg[sse][dirname][l] < 0.1 and avg[sse][dirname][l] > -0.1:
                            print(f"0({int(counter[sse][dirname][l]/counter[sse][dirname][0]*100)}%)", end=" ")
                        else:
                            print(f"{avg[sse][dirname][l]:.1f}({int(counter[sse][dirname][l]/counter[sse][dirname][0]*100)})", end=" ")
                print()
    print()

aavg = {}

# count avg for all PINS cumulative
for sse in avg:
    aavg[sse] = []
    is_in = True
    res_counter = 0
    while is_in:
        is_in = False
        aavg_ch = 0
        dom_counter = 0
        for dirname in avg[sse]:
             if len(avg[sse][dirname]) > res_counter:
                  dom_counter += 1
                  aavg_ch += avg[sse][dirname][res_counter]
                  is_in = True
        if is_in is True:
            aavg[sse].append(aavg_ch/dom_counter)
            res_counter += 1

# avg for all PINS cumulative
print("2/ avg for all PINS cumulative ================================")
if True:
    for sse in aavg:
        if sse in sses_list:
            print(f"cumulative avg {sses_dict[sse]}:", end=" ")
            for r in aavg[sse]:
                print(f"{r:.1f}", end=" ")
            print()
    print()

# deviation counting
deviation = {}
for sse in data:
    deviation[sse] = {}
    for dirname in data[sse]:
        deviation[sse][dirname] = []
        is_in = True
        res_counter = 0
        while is_in:
            is_in = False
            avg_deviation = 0
            dom_counter = 0
            for filename in data[sse][dirname]:
                 if len(data[sse][dirname][filename]) > res_counter:
                      dom_counter += 1
                      avg_deviation += abs((data[sse][dirname][filename][res_counter] - avg[sse][dirname][res_counter])*counter[sse][dirname][res_counter]/counter[sse][dirname][0])
                      is_in = True
            if is_in is True:
                 deviation[sse][dirname].append(avg_deviation / dom_counter)
            res_counter += 1


if True:
    print("3/ deviation values per SSE per PIN ================================")
    for sse in deviation:
        if sse in sses_list:
            print(f"sse: {sses_dict[sse]}")
            for dirname in dir_list:
                print(f"deviation {dirname} {len(deviation[sse][dirname])}:", end=" ")
                c =0
                for r in deviation[sse][dirname]:
                    print(f"{r:.1f}({counter[sse][dirname][c]/counter[sse][dirname][0]})", end=" ")
                    c += 1
                print()
    print()
#sys.exit()

if True:
    print("4/ deviation PINx vs PINy ================================= ")
    print("avg deviation ", end=" ")
    for dirname in dir_list:
        print(dirname, end=" ")
    print()
    for sse in deviation:
        if sse in sses_list:
            print(f"deviation {sses_dict[sse]}:", end=" ")
            for dirname in dir_list:
                print(f"{sum(deviation[sse][dirname]):.1f}", end=" ")
            print()

dif = {}
for sse in data:
    dif[sse] = {}
    for dirname in data[sse]:
        dif[sse][dirname] = []
        is_in = True
        res_counter = 0
        while is_in:
            is_in = False
            avg_dis = 0
            dom_counter = 0
            for filename in data[sse][dirname]:
                 if len(data[sse][dirname][filename]) > res_counter:
                      dom_counter += 1
                      avg_dis += (data[sse][dirname][filename][res_counter] - aavg[sse][res_counter]) * (data[sse][dirname][filename][res_counter] - aavg[sse][res_counter])
                      is_in = True
            if is_in is True:
                 dif[sse][dirname].append(avg_dis/dom_counter)
            res_counter += 1

if True:
    print()
    print("5/ difference  ===================================================================== ")
    print("avg difference", end=" ")
    for dirname in dir_list:
        print(dirname, end=" ")
    print()
    for sse in deviation:
        if sse in sses_list:
            print(f"{sses_dict[sse]}", end=" | ")
            for dirname in dir_list:
                print(f"{sum(dif[sse][dirname]):.2f}", end=" ")
            print()
# ------------------

dif2 = {}
for sse in data:
    dif2[sse] = {}
    for dirname in data[sse]:
        dif2[sse][dirname] = []
        is_in = True
        res_counter = 0
        while is_in:
            is_in = False
            avg_dis = 0
            dom_counter = 0
            if len(avg[sse][dirname]) > res_counter:
                dom_counter += 1
                avg_dis += abs(avg[sse][dirname][res_counter]-aavg[sse][res_counter])*(counter[sse][dirname][res_counter])/(counter[sse][dirname][0])
                is_in = True
            if is_in is True:
                dif2[sse][dirname].append(avg_dis/dom_counter)
            res_counter += 1


if True:
    print()
    print("6/ difference2  ===================================================================== ")
    print(" avg difference2", end=" ")
    for dirname in dir_list:
        print(dirname, end=" ")
    print()
    for sse in deviation:
        if sse in sses_list:
            print(f"{sses_dict[sse]}", end=" | ")
            for dirname in dir_list:
                print(f"{sum(dif2[sse][dirname]):.2f}", end=" ")
            print()

dir_list2 = ["PIN1", "PIN2", "PIN3", "PIN4", "PIN5", "PIN6", "PIN7", "PIN8"]
# dir_list2 = ["PIN3", "PIN4"]

pin_diff = {}
print()
print("diff2", end=" ")
for dirname2 in dir_list2:
    print(dirname2, end=" ")
print()

for dirname1 in dir_list2:
    print(f"<tr><td>{dirname1}</td>", end=" ")
    for dirname2 in dir_list2:
        s = 0
        c = 0
        for sse in avg:
            for l in range(len(avg[sse][dirname1])):
                if len(avg[sse][dirname2]) > l:
                    s += (abs(avg[sse][dirname1][l] - avg[sse][dirname2][l])) * (counter[sse][dirname1][l] + counter[sse][dirname2][l]) / (counter[sse][dirname1][0] + counter[sse][dirname2][0])
                    c += 1
        print(f"<td>{s:.2f}</td>", end=" ")
    print("</tr>")

