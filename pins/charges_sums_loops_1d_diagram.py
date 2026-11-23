#!/usr/bin/env python3
# tool which generate for given pin directory its standart graphs with up/down part charges sum

import argparse
import time
import os
import json
import svgwrite
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def up_part(ch_data, part1, part2, text = ""):
    up = []
    down = []
    res_min = 10000
    res_max = 0
    for res in ch_data:
        if int(res) < res_min:
            res_min = int(res)
        if int(res) > res_max:
            res_max = int(res)
    res_middle = res_min + (res_max - res_min) / 2
    aux = []
    res = res_min
    while res < res_middle:
        aux.append((ch_data[str(res)]))
        res += 1

    if part1 == "up":
        up = aux
    else:
        down = aux
    aux2 = []
    while res <= res_max:
        aux2.append((ch_data[str(res)]))
        res += 1
    if part2 == "up":
        up += aux2
    else:
        down += aux2

    return up, down


# r_upup_part2(name, data, sses_list, data_all_res, 1, 2, "up", name)
def sses_part(name, data, sses_list, data_all_res, min_index, max_index, part1, part2, text):
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = None
        for i in data[sses_list[min_index]][name]:
            if res_min is None or res_min > int(i):
                res_min = int(i)
        res_min += 1

        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max is None or res_max < int(i):
                res_max = int(i)

        ch_data = {}
        res = res_min-1
        while res <= res_max:
            ch_data[str(res)] = (data_all_res[res])
            res += 1
#        print("ch_data", text, ch_data, "----------------")
        print("ch_data", name, text, ch_data, "----------------")
        return up_part(ch_data, part1, part2, text)
    else:
        return [], []


#    up_part2(name, data, sses_list, data_all_res, 0, 1, "down", name)
def up_part2(name, data, sses_list, data_all_res, min_index, max_index, part, text):
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = 1
        for i in data[sses_list[min_index]][name]:
            if res_min is None or res_min < int(i):
                res_min = int(i)
        res_min += 1

        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max is None or res_max > int(i):
                res_max = int(i)

        ch_data = {}
        res = res_min
        while res < res_max:
            ch_data[str(res)] = (data_all_res[res])
            res += 1
        print("ch_data", name, text, ch_data, "----------------")
        if part == "ud":
            return up_part(ch_data, "up", "down", text)
        elif part == "du":
            return up_part(ch_data, "down", "up", text)
        elif part == "up":
            return up_part(ch_data, part, part, text)
        elif part == "down":
            return up_part(ch_data, part, part, text)
        else:
            print("part je spatne", part)
            exit
    else:
        return [], []


def up_part_start(name, data, sses_list, data_all_res, sindex, part, text):
    if name in data[sses_list[sindex]]:
        res_min = 1000
        for i in data[sses_list[sindex]][name]:
            if res_min is None or res_min > int(i):
                res_min = int(i)

        ch_data = {}
        res = 1
        while res < res_min:
            ch_data[str(res)] = data_all_res[res]
            res += 1
        print("ch_data", name, text, ch_data, "----------------")

        return up_part(ch_data, part, part, text)
    else:
        return [], []


def up_part_end(name, data, sses_list, data_all_res, sindex, part, text):
    if name in data[sses_list[sindex]]:
        res_max = 1
        for i in data[sses_list[sindex]][name]:
            if res_max is None or res_max < int(i):
                res_max = int(i)

        ch_data = {}
        res = res_max + 1
        while res in data_all_res:
            ch_data[str(res)] = data_all_res[res]
            res += 1
        print("ch_data", name, text, ch_data, "----------------")
        return up_part(ch_data, part, part, text)
    else:
        return [], []


def aall(name, data, sses_list, data_all_res, sindex, part, text):
    if name in data[sses_list[sindex]]:
        ch_data = {}
        res = 0
        while res in data_all_res:
            ch_data[str(res)] = data_all_res[res]
            res += 1
            
        slist = []
        ssum = 0
        for i in data_all_res:
            ssum += data_all_res[i]
#            slist.append(i)
            
        print("ch1_data", name, text, ssum, "----------------")
#        print("ch1_data", name, text, data_all_res, "----------------")
        
        return up_part(ch_data, part, part, text)
    else:
        return [], []





def graph_loops(data, sses_order, fn_order, avg, first, data_all_res, mapping_names):
    # generate diagram and record
    name = fn_order[0]

    c_sses = 0
    for sses in sses_order:
         c_sses += 1

    #      sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"] # small dataset(45 pcs)
    sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]  # BIG DATASET
    sum_up = []
    sum_down = []

    r_up, r_down = up_part_start(name, data, sses_list, data_all_res, 0, "up", name)
    sum_up += r_up
    sum_down += r_down
    print(f" st0: {sum(r_up)} ", end="|")

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 0, 1, "down", name)
    sum_up += r_up
    sum_down += r_down
#    print(f"{r_up}(0-1)", end="| ")

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 1, 2, "up", name)
    sum_up += r_up
    sum_down += r_down
    print(f"1-2: {sum(r_up)} ", end="|")

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 2, 3, "down", name)
    sum_up += r_up
    sum_down += r_down
#    print(f"{r_up}(2-3)", end="| ")

#    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 3, 4, "up", name)
#    sum_up += r_up
#    sum_down += r_down

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 4, 5, "up", name)
    sum_up += r_up
    sum_down += r_down
    print(f"3-4: {sum(r_up)} ", end="|")

    # big loop
    big_loop_up, big_loop_down = up_part2(name, data, sses_list, data_all_res, 5, 6, "down", name)
#    sum_up += r_up
#    sum_down += r_down
#    pom_up += sum(r_up)
#    pom_down += sum(r_down)

#    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 6, 7, "up", name)
#    sum_up += r_up
#    sum_down += r_down

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 7, 8, "up", name)
    sum_up += r_up
    sum_down += r_down
    print(f"7-8: {sum(r_up)} ", end="|")

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 8, 9, "down", name)
    sum_up += r_up
    sum_down += r_down
#    print(f"{r_up}", end="| ")

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 9, 10, "up", name)
    sum_up += r_up
    sum_down += r_down
    print(f"9-10: {sum(r_up)} ", end="|")

#    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 10, 11, "up", name)
#    sum_up += r_up
#    sum_down += r_down

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 11, 12, "down", name)
    sum_up += r_up
    sum_down += r_down
#    print(f"{r_up}", end="| ")

    r_up, r_down = up_part_end(name, data, sses_list, data_all_res, 12, "up", name)
    sum_up += r_up
    sum_down += r_down
    print(f"e12: {sum(r_up)} ", end="||")

    sname = name.split("-")
    # print loops only charges
    res = [mapping_names[name], int(sum(sum_up)), int(sum(sum_down)), sname[1], sum(big_loop_down), int(min(list(data[sses_list[6]][name].keys()))) - int(max(list(data[sses_list[5]][name].keys())))]
    print("res je", res[0],"++", res[1], res[2], name )
    return res


def graph_full(data, sses_order, fn_order, avg, first, data_all_res, mapping_names):
    # generate diagram and record
    name = fn_order[0]

    c_sses = 0
    for sses in sses_order:
         c_sses += 1

    #      sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"] # small dataset(45 pcs)
    sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]  # BIG DATASET
    sum_up = []
    sum_down = []
    vypis_sum_up = []
    vypis_sum_down = []

    # start 
    r_up, r_down = up_part_start(name, data, sses_list, data_all_res, 0, "up", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up),"start"])
    vypis_sum_down.append([sum(r_down),"start"])
#    print(f" st0: {sum(r_up)} ", end="|")

    # SSE 1
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 0, 0, "up", "down", "SSE1")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE1"])
    vypis_sum_down.append([sum(r_down), "SSE1"])
#    print(f"| sse1: {sum(r_up)}/{sum(r_down)} ", end="||")

    # SSE 1 <-> 2 connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 0, 1, "down", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"{r_up}(0-1)", end="| ")

    # SSE 2
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 1, 1, "down", "up",  "SSE2")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE2"])
    vypis_sum_down.append([sum(r_down), "SSE2"])
#    print(f"| sse2: {sum(r_up)}/{sum(r_down)} ", end="||")

    # SSE 2 <-> 3 connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 1, 2, "up", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"2-3: {sum(r_up)} ", end="|")

    # SSE 3
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 2, 2, "up", "down", "SSE3")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE3"])
    vypis_sum_down.append([sum(r_down), "SSE3"])
#    print(f"| sse1: {sum(r_up)}/{sum(r_down)} ", end="||")

    # SSE 3 <-> 4a connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 2, 3, "down", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"{r_up}(2-3)", end="| ")

    # SSE 4a
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 3, 3, "down", "down", "SSE4a")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE4a"])
    vypis_sum_down.append([sum(r_down), "SSE4a"])
#    print(f"| sse1: {sum(r_up)}/{sum(r_down)} ", end="||")

    # SSE 4a <-> 4b connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 3, 4, "du", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))

    # SSE 4b
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 4, 4, "up", "up", "SSE4b")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE4b"])
    vypis_sum_down.append([sum(r_down), "SSE4b"])
#    print(f"| sse1: {sum(r_up)}/{sum(r_down)} ", end="||")

    # SSE 4b <-> 5 connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 4, 5, "up", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"3-4: {sum(r_up)} ", end="|")

    # SSE 5
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 5, 5, "up", "down", "SSE5")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE5"])
    vypis_sum_down.append([sum(r_down), "SSE5"])

    # big loop 5 <-> 6a
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 5, 6, "down", name)
    big_loop_up = r_up
    big_loop_down = r_down
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))

    # SSE 6a
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 6, 6, "down", "down", "SSE6a")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE6a"])
    vypis_sum_down.append([sum(r_down), "SSE6a"])

    # SSE 6a <-> 6b connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 6, 7, "du", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))

    # SSE 6b
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 7, 7, "up", "up", "SSE6b")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE6b"])
    vypis_sum_down.append([sum(r_down), "SSE6b"])

    # SSE 6b <-> 7 connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 7, 8, "up", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"7-8: {sum(r_up)} ", end="|")

    # SSE 7
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 8, 8, "up", "down", "SSE7")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE7"])
    vypis_sum_down.append([sum(r_down), "SSE7"])

    # SSE 7 <-> 8 connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 8, 9, "down", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"{r_up}", end="| ")

    # SSE 8
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 9, 9, "down", "up", "SSE8")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE8"])
    vypis_sum_down.append([sum(r_down), "SSE8"])

    # SSE 8 <-> 9a connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 9, 10, "up", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"9-10: {sum(r_up)} ", end="|")

    # SSE 9a
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 10, 10, "up", "up", "SSE9a")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE9a"])
    vypis_sum_down.append([sum(r_down), "SSE9a"])

    # SSE 9a <-> 9b connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 10, 11, "ud", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))

    # SSE 9b
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 11, 11, "down", "down", "SSE9b")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE9b"])
    vypis_sum_down.append([sum(r_down), "SSE9b"])

    # SSE 9b <-> 10 connection
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 11, 12, "down", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"{r_up}", end="| ")

    # SSE 10
    r_up, r_down = sses_part(name, data, sses_list, data_all_res, 12, 12, "down", "up", "SSE10")
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append([sum(r_up), "SSE10"])
    vypis_sum_down.append([sum(r_down), "SSE10"])

    # end
    r_up, r_down = up_part_end(name, data, sses_list, data_all_res, 12, "up", name)
    sum_up += r_up
    sum_down += r_down
    vypis_sum_up.append(sum(r_up))
    vypis_sum_down.append(sum(r_down))
#    print(f"e12: {sum(r_up)} ", end="||")
    
    aall(name, data, sses_list, data_all_res, 12, "up", name)
    

    sname = name.split("-")
    # print loops only charges
    res = [mapping_names[name], int(sum(sum_up)), int(sum(sum_down)), sname[1], sum(big_loop_down), int(min(list(data[sses_list[6]][name].keys()))) - int(max(list(data[sses_list[5]][name].keys())))]
    print("res je", res[0],"++", res[1], res[2], name )
    print("sum_up", name,  vypis_sum_up)
    print("sum_down", name, vypis_sum_down)
    return res


def avg(vlist):
    vcounter = 0
    vsum = 0
    for i in vlist:
        vcounter += 1
        vsum += i
    return(vsum /vcounter)

# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument('--charges_dir', required=True, action="store",
                    help="charges data dir")
parser.add_argument("--name", help="name")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = os.path.basename(os.path.normpath(args.directory))
new_dir = f"charges-{name}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)
print("new_dir", new_dir, " - created")

data = {}      # charges data
fn_order = []
data_all_res = {}
pin_mapping = {}
for dirname in os.listdir(args.directory):
  if os.path.isdir(dirname):
      continue
  print("pin directory name:", dirname[0:4], ": ", dirname)
  for filename in os.listdir(f"{args.directory}/{dirname}"):
    pin_mapping[filename] = dirname[0:4]
    # go through all PINs in PIN directory and create record for them

    # read PIN json and read residues charges in SSEs
    print("read", f"{args.directory}/{dirname}/{filename}")
    with open(f"{args.directory}/{dirname}/{filename}", "r") as fd:
        jdata = json.load(fd)
        fn_order.append(filename)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in data:
                    data[sse] = {}
                data[sse][filename] = jdata["charges"][sse]

    # read PIN json and read residues charges in parts outside SSEs
    sfilename = filename.split("-")
    ac_filename = f"{sfilename[1]}.alphacharges"
    with open(f"{args.charges_dir}/{ac_filename}", "rb") as fd:
        ac_res = {}
        r = fd.read()

        line = ""
        start = False

        for c in r:
            if chr(c) != "\n":
                # read separate character and put them to to lines
                line = line + chr(c)
            else:
                if "ATOM" in line:
                    # the first line contains the rest of old stuff at the beginning, need to be moved
                    pos = line.find("ATOM")
                    line = line[pos:]
                if line.startswith("ATOM") and start:
                    if int(line[21:26]) not in ac_res:
                        ac_res[int(line[21:26])] = []
                    ac_res[int(line[21:26])].append(float(line[55:62]))

                line = ""
            if line.startswith("TER"):
                # read only stuff between "TER" signs
                if not start:
                    start = True
                else:
                    break

    # data aggregation
    data_all_res[filename] = {}
    for r in ac_res:
        s = 0
        for ss in ac_res[r]:
            s += ss
        data_all_res[filename][r] = s*100
    for sse in data:
        if filename in data[sse]:
            for res in data[sse][filename]:
                data[sse][filename][res] = data_all_res[filename][int(res)]

pin_mapping_order = sorted(pin_mapping.items(), key=lambda x: x)
fn_order = []

for i in pin_mapping:
    fn_order.append(i)
sses_order = list(data.keys())

# for all PINs generate the output html record and diagram
records = []

for f in fn_order:
    # generate diagram and recors
    r = graph_full(data, sses_order, [f], 2, True, data_all_res[f], pin_mapping)
    if r is not False:
        records.append(r)

const = 2
# do diagram
width = 1400
height = 1400

# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-upper_loops-down_loops222-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))

csum = {}
ccounter = {}


for sline in records:
        if (sline[0] == "pin1") or (sline[0] == "PIN1"):
            color = "red"
            c = 0
        elif (sline[0] == "pin2") or (sline[0] == "PIN2"):
            color = "green"
            c = 1
#            dwg.add(dwg.text(f"{sline[2]}/{sline[1]}", insert=((width/2+int(sline[2])/2-1)*const, (height/2-int(sline[1])/2+1)*const), color="black", font_size=20, style="text-anchor:center"))

        elif (sline[0] == "pin3") or (sline[0] == "PIN3"):
            color = "blue"
            c = 2
        elif (sline[0] == "pin4") or (sline[0] == "PIN4"):
            color = "pink"
            c = 3
        elif (sline[0] == "pin5") or (sline[0] == "PIN5"):
            color = "magenta"
#            dwg.add(dwg.text(f"___{sline[2]}/{sline[1]}", insert=((width/2+int(sline[2])/2-1)*const, (height/2-int(sline[1])/2+1)*const), color="black", font_size=8, style="text-anchor:center"))
            c = 4
        elif (sline[0] == "pin6") or (sline[0] == "PIN6"):
            color = "cyan"
            c = 5
        elif (sline[0] == "pin7") or (sline[0] == "PIN7"):
            color = "brown"
            c = 6
        elif (sline[0] == "pin8") or (sline[0] == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        csum[c] += sline[2]
        ccounter[c] += 1
#        if c == 4:
        dwg.add(dwg.line(start=((width/2+int(sline[2])/2-2)*const, (height/2-int(sline[1])/2+2)*const), end=((width/2+int(sline[2])/2+3)*const, (height/2-int(sline[1])/2-3)*const), stroke_width=4, stroke=color))




#dwg.add(dwg.line(start=((width/2+500)*const, (height/2-500)*const), end=((width/2-500)*const, (height/2+500)*const), stroke="black", stroke_width=2))

#dwg.add(dwg.text(f"{0}", insert=((width/2+3)*const, (height/2-3)*const), color="white", font_size=12, style="text-anchor:center"))

#diagonaly
dwg.add(dwg.line(start=((width/2-650)*const, (height/2-500)*const), end=((width/2+350)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-600)*const, (height/2-500)*const), end=((width/2+400)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-550)*const, (height/2-500)*const), end=((width/2+450)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-500)*const), end=((width/2+500)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-450)*const, (height/2-500)*const), end=((width/2+550)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-400)*const, (height/2-500)*const), end=((width/2+600)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-350)*const, (height/2-500)*const), end=((width/2+650)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-300)*const, (height/2-500)*const), end=((width/2+700)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-250)*const, (height/2-500)*const), end=((width/2+750)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-200)*const, (height/2-500)*const), end=((width/2+800)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-150)*const, (height/2-500)*const), end=((width/2+850)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-100)*const, (height/2-500)*const), end=((width/2+900)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2-50)*const, (height/2-500)*const), end=((width/2+950)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2+1000)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2+50)*const, (height/2-500)*const), end=((width/2+1050)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
dwg.add(dwg.line(start=((width/2+100)*const, (height/2-500)*const), end=((width/2+1100)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))


#hlavni kriz
dwg.add(dwg.line(start=((width/2-500)*const, (height/2)*const), end=((width/2+1000)*const, (height/2)*const), stroke="black", stroke_width=2))
dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2)*const, (height/2+500)*const), stroke="black", stroke_width=2))
#dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2)*const, (height/2+500)*const), stroke="black", stroke_width=2))

#vedlejsi osy
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-100/2)*const), end=((width/2+1000)*const, (height/2-100/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-200/2)*const), end=((width/2+1000)*const, (height/2-200/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-300/2)*const), end=((width/2+1000)*const, (height/2-300/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-400/2)*const), end=((width/2+1000)*const, (height/2-400/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-500/2)*const), end=((width/2+1000)*const, (height/2-500/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-600/2)*const), end=((width/2+1000)*const, (height/2-600/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+100/2)*const), end=((width/2+1000)*const, (height/2+100/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+200/2)*const), end=((width/2+1000)*const, (height/2+200/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+300/2)*const), end=((width/2+1000)*const, (height/2+300/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+400/2)*const), end=((width/2+1000)*const, (height/2+400/2)*const), stroke="grey"))

dwg.add(dwg.line(start=((width/2-600/2)*const, (height/2-500)*const), end=((width/2-600/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500/2)*const, (height/2-500)*const), end=((width/2-500/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-400/2)*const, (height/2-500)*const), end=((width/2-400/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-300/2)*const, (height/2-500)*const), end=((width/2-300/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-200/2)*const, (height/2-500)*const), end=((width/2-200/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-100/2)*const, (height/2-500)*const), end=((width/2-100/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-0/2)*const, (height/2-500)*const), end=((width/2-0/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+100/2)*const, (height/2-500)*const), end=((width/2+100/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+200/2)*const, (height/2-500)*const), end=((width/2+200/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+300/2)*const, (height/2-500)*const), end=((width/2+300/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+400/2)*const, (height/2-500)*const), end=((width/2+400/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+500/2)*const, (height/2-500)*const), end=((width/2+500/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+600/2)*const, (height/2-500)*const), end=((width/2+600/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+700/2)*const, (height/2-500)*const), end=((width/2+700/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+800/2)*const, (height/2-500)*const), end=((width/2+800/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+900/2)*const, (height/2-500)*const), end=((width/2+900/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+1000/2)*const, (height/2-500)*const), end=((width/2+1000/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+1100/2)*const, (height/2-500)*const), end=((width/2+1100/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+1200/2)*const, (height/2-500)*const), end=((width/2+1200/2)*const, (height/2+500)*const), stroke="grey"))

dwg.add(dwg.text(f"{1}", insert=((width/2+3)*const, (height/2-100/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{2}", insert=((width/2+3)*const, (height/2-200/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{3}", insert=((width/2+3)*const, (height/2-300/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{4}", insert=((width/2+3)*const, (height/2-400/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
#dwg.add(dwg.text(f"{5}", insert=((width/2+3)*const, (height/2-500/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
#dwg.add(dwg.text(f"{6}", insert=((width/2+3)*const, (height/2-600/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{1}", insert=((width/2+3)*const, (height/2+100/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-2}", insert=((width/2+3)*const, (height/2+200/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-4}", insert=((width/2+3)*const, (height/2+400/2-3)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-3}", insert=((width/2+3)*const, (height/2+300/2-3)*const), color="white", font_size=40, style="text-anchor:center"))

#dwg.add(dwg.text("charge of", insert=((width/2+2)*const, (height/2-300/2-3-183)*const), font_size=46, style="text-anchor:center"))
#dwg.add(dwg.text("non central loop side", insert=((width/2+2)*const, (height/2-300/2-3-183+27)*const), font_size=46, style="text-anchor:center"))

dwg.add(dwg.text(f"{-6}", insert=((width/2-600/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-5}", insert=((width/2-500/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-4}", insert=((width/2-400/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-3}", insert=((width/2-300/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-2}", insert=((width/2-200/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{-1}", insert=((width/2-100/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{0}", insert=((width/2+0/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{1}", insert=((width/2+100/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{2}", insert=((width/2+200/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{3}", insert=((width/2+300/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{4}", insert=((width/2+400/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{5}", insert=((width/2+500/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{6}", insert=((width/2+600/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{7}", insert=((width/2+700/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{8}", insert=((width/2+800/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{9}", insert=((width/2+900/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{10}", insert=((width/2+1000/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{11}", insert=((width/2+1100/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
dwg.add(dwg.text(f"{12}", insert=((width/2+1200/2+3)*const, (height/2-5)*const), color="white", font_size=40, style="text-anchor:center"))
#dwg.add(dwg.text("charge of", insert=((width/2+400/2+2+300)*const, (height/2+9+11)*const), color="white", font_size=46, style="text-anchor:center"))
#dwg.add(dwg.text("central loop side", insert=((width/2+400/2+2+300)*const, (height/2-3+9+38)*const), color="white", font_size=46, style="text-anchor:center"))
pcRotate = 'rotate(45,%s, %s)' % ((width/2+600/2+3)*const, (height/2+200/2-3)*const)
pcRotate = 'rotate(45,%s, %s)' % ((width/2+400/2+3)*const, (height/2+400/2-3)*const)

#dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))
if True:
    y=-180
    x= 350
    dwg.add(dwg.rect(insert=((width/2+185+x)*const, (height/2-177+y)*const), size=(33*const,330*const), fill="white"))
#    dwg.add(dwg.rect(insert=((width/2+185+x)*const, (height/2-177+y)*const), size=(110*const,330*const), fill="white"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+y)*const), stroke_width=24, stroke="red"))
    dwg.add(dwg.text("PIN1", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+40+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+40+y)*const), stroke_width=24, stroke="green"))
    dwg.add(dwg.text("PIN2", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+40+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+80+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+80+y)*const), stroke_width=24, stroke="blue"))
    dwg.add(dwg.text("PIN3", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+80+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+120+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+120+y)*const), stroke_width=24, stroke="pink"))
    dwg.add(dwg.text("PIN4", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+120+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+160+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+160+y)*const), stroke_width=24, stroke="magenta"))
    dwg.add(dwg.text("PIN5", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+160+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+200+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+200+y)*const), stroke_width=24, stroke="cyan"))
    dwg.add(dwg.text("PIN6", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+200+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+240+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+240+y)*const), stroke_width=24, stroke="brown"))
    dwg.add(dwg.text("PIN7", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+240+y)*const), font_size=64, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-5+400/2+x)*const, (height/2-300/2+5+280+y)*const), end=((width/2+400/2+5+x)*const, (height/2-300/2-5+280+y)*const), stroke_width=24, stroke="orange"))
    dwg.add(dwg.text("PIN8", insert=((width/2+400/2+2+12+x)*const, (height/2-300/2+9+280+y)*const), font_size=64, style="text-anchor:center"))
dwg.save()


# for all PINs generate the output html record and diagram
records = []

for f in fn_order:
    # generate diagram and recors
    r = graph_loops(data, sses_order, [f], 2, True, data_all_res[f], pin_mapping)
    if r is not False:
        records.append(r)

const = 2
# do diagram
width = 1400
height = 1400

# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-upper_loops-down_loops-{name}-2.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))

csum = {}
ccounter = {}


for sline in records:
        if (sline[0] == "pin1") or (sline[0] == "PIN1"):
            color = "red"
            c = 0
        elif (sline[0] == "pin2") or (sline[0] == "PIN2"):
            color = "green"
            c = 1
#            dwg.add(dwg.text(f"{sline[2]}/{sline[1]}", insert=((width/2+int(sline[2])/2-1)*const, (height/2-int(sline[1])/2+1)*const), color="black", font_size=20, style="text-anchor:center"))

        elif (sline[0] == "pin3") or (sline[0] == "PIN3"):
            color = "blue"
            c = 2
        elif (sline[0] == "pin4") or (sline[0] == "PIN4"):
            color = "pink"
            c = 3
        elif (sline[0] == "pin5") or (sline[0] == "PIN5"):
            color = "magenta"
#            dwg.add(dwg.text(f"___{sline[2]}/{sline[1]}", insert=((width/2+int(sline[2])/2-1)*const, (height/2-int(sline[1])/2+1)*const), color="black", font_size=8, style="text-anchor:center"))
            c = 4
        elif (sline[0] == "pin6") or (sline[0] == "PIN6"):
            color = "cyan"
            c = 5
        elif (sline[0] == "pin7") or (sline[0] == "PIN7"):
            color = "brown"
            c = 6
        elif (sline[0] == "pin8") or (sline[0] == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        csum[c] += sline[2]
        ccounter[c] += 1
#        if c == 4:
        dwg.add(dwg.line(start=((width/2+int(sline[2])/2-1)*const, (height/2-int(sline[1])/2+1)*const), end=((width/2+int(sline[2])/2+2)*const, (height/2-int(sline[1])/2-2)*const), stroke_width=8, stroke=color))




#dwg.add(dwg.line(start=((width/2+500)*const, (height/2-500)*const), end=((width/2-500)*const, (height/2+500)*const), stroke="black", stroke_width=2))

#dwg.add(dwg.text(f"{0}", insert=((width/2+3)*const, (height/2-3)*const), color="white", font_size=12, style="text-anchor:center"))

#diagonaly
#dwg.add(dwg.line(start=((width/2-650)*const, (height/2-500)*const), end=((width/2+350)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-600)*const, (height/2-500)*const), end=((width/2+400)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-550)*const, (height/2-500)*const), end=((width/2+450)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-500)*const, (height/2-500)*const), end=((width/2+500)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-450)*const, (height/2-500)*const), end=((width/2+550)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-400)*const, (height/2-500)*const), end=((width/2+600)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-350)*const, (height/2-500)*const), end=((width/2+650)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-300)*const, (height/2-500)*const), end=((width/2+700)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-250)*const, (height/2-500)*const), end=((width/2+750)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-200)*const, (height/2-500)*const), end=((width/2+800)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-150)*const, (height/2-500)*const), end=((width/2+850)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-100)*const, (height/2-500)*const), end=((width/2+900)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2-50)*const, (height/2-500)*const), end=((width/2+950)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2+1000)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2+50)*const, (height/2-500)*const), end=((width/2+1050)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
#dwg.add(dwg.line(start=((width/2+100)*const, (height/2-500)*const), end=((width/2+1100)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))


#hlavni kriz
dwg.add(dwg.line(start=((width/2-500)*const, (height/2)*const), end=((width/2+1000)*const, (height/2)*const), stroke="black", stroke_width=2))
dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2)*const, (height/2+500)*const), stroke="black", stroke_width=2))
#dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2)*const, (height/2+500)*const), stroke="black", stroke_width=2))

#vedlejsi osy
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-100/2)*const), end=((width/2+1000)*const, (height/2-100/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-200/2)*const), end=((width/2+1000)*const, (height/2-200/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-300/2)*const), end=((width/2+1000)*const, (height/2-300/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-400/2)*const), end=((width/2+1000)*const, (height/2-400/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-500/2)*const), end=((width/2+1000)*const, (height/2-500/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2-600/2)*const), end=((width/2+1000)*const, (height/2-600/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+100/2)*const), end=((width/2+1000)*const, (height/2+100/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+200/2)*const), end=((width/2+1000)*const, (height/2+200/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+300/2)*const), end=((width/2+1000)*const, (height/2+300/2)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500)*const, (height/2+400/2)*const), end=((width/2+1000)*const, (height/2+400/2)*const), stroke="grey"))

dwg.add(dwg.line(start=((width/2-600/2)*const, (height/2-500)*const), end=((width/2-600/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-500/2)*const, (height/2-500)*const), end=((width/2-500/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-400/2)*const, (height/2-500)*const), end=((width/2-400/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-300/2)*const, (height/2-500)*const), end=((width/2-300/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-200/2)*const, (height/2-500)*const), end=((width/2-200/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-100/2)*const, (height/2-500)*const), end=((width/2-100/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2-0/2)*const, (height/2-500)*const), end=((width/2-0/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+100/2)*const, (height/2-500)*const), end=((width/2+100/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+200/2)*const, (height/2-500)*const), end=((width/2+200/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+300/2)*const, (height/2-500)*const), end=((width/2+300/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+400/2)*const, (height/2-500)*const), end=((width/2+400/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+500/2)*const, (height/2-500)*const), end=((width/2+500/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+600/2)*const, (height/2-500)*const), end=((width/2+600/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+700/2)*const, (height/2-500)*const), end=((width/2+700/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+800/2)*const, (height/2-500)*const), end=((width/2+800/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+900/2)*const, (height/2-500)*const), end=((width/2+900/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+1000/2)*const, (height/2-500)*const), end=((width/2+1000/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+1100/2)*const, (height/2-500)*const), end=((width/2+1100/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.line(start=((width/2+1200/2)*const, (height/2-500)*const), end=((width/2+1200/2)*const, (height/2+500)*const), stroke="grey"))

dwg.add(dwg.text(f"{1}", insert=((width/2+3)*const, (height/2-100/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{2}", insert=((width/2+3)*const, (height/2-200/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{3}", insert=((width/2+3)*const, (height/2-300/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{4}", insert=((width/2+3)*const, (height/2-400/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{5}", insert=((width/2+3)*const, (height/2-500/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
#dwg.add(dwg.text(f"{6}", insert=((width/2+3)*const, (height/2-600/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{1}", insert=((width/2+3)*const, (height/2+100/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-2}", insert=((width/2+3)*const, (height/2+200/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-4}", insert=((width/2+3)*const, (height/2+400/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-3}", insert=((width/2+3)*const, (height/2+300/2-3)*const), color="white", font_size=16, style="text-anchor:center"))

#dwg.add(dwg.text("charge of small loops", insert=((width/2+2)*const, (height/2-100/2-3-133+14)*const), font_size=20, style="text-anchor:center"))
#dwg.add(dwg.text("on non central loop side", insert=((width/2+2)*const, (height/2-100/2-3-133+25)*const), font_size=20, style="text-anchor:center"))

dwg.add(dwg.text(f"{-6}", insert=((width/2-600/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-5}", insert=((width/2-500/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-4}", insert=((width/2-400/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-3}", insert=((width/2-300/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-2}", insert=((width/2-200/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{-1}", insert=((width/2-100/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{0}", insert=((width/2+0/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{1}", insert=((width/2+100/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{2}", insert=((width/2+200/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{3}", insert=((width/2+300/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{4}", insert=((width/2+400/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{5}", insert=((width/2+500/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{6}", insert=((width/2+600/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{7}", insert=((width/2+700/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{8}", insert=((width/2+800/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{9}", insert=((width/2+900/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{10}", insert=((width/2+1000/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{11}", insert=((width/2+1100/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
dwg.add(dwg.text(f"{12}", insert=((width/2+1200/2+3)*const, (height/2-5)*const), color="white", font_size=16, style="text-anchor:center"))
#dwg.add(dwg.text("charge of small loops", insert=((width/2+200/2+2+80)*const, (height/2+9)*const), color="white", font_size=20, style="text-anchor:center"))
#dwg.add(dwg.text("on central loop side", insert=((width/2+200/2+2+80)*const, (height/2+9+11)*const), color="white", font_size=20, style="text-anchor:center"))
pcRotate = 'rotate(45,%s, %s)' % ((width/2+600/2+3)*const, (height/2+200/2-3)*const)
pcRotate = 'rotate(45,%s, %s)' % ((width/2+400/2+3)*const, (height/2+400/2-3)*const)

#dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))
v1 = "full"
if True:
    y=-25
    x=25
    dwg.add(dwg.rect(insert=((width/2-1+400/2-5+x)*const, (height/2-300/2+1-7+y)*const), size=(35*const,80 *const), fill="white"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+y)*const), stroke_width=8, stroke="red"))
    dwg.add(dwg.text("PIN1", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+10+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+10+y)*const), stroke_width=8, stroke="green"))
    dwg.add(dwg.text("PIN2", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+10+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+20+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+20+y)*const), stroke_width=8, stroke="blue"))
    dwg.add(dwg.text("PIN3", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+20+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+30+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+30+y)*const), stroke_width=8, stroke="pink"))
    dwg.add(dwg.text("PIN4", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+30+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+40+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+40+y)*const), stroke_width=8, stroke="magenta"))
    dwg.add(dwg.text("PIN5", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+40+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+50+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+50+y)*const), stroke_width=8, stroke="cyan"))
    dwg.add(dwg.text("PIN6", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+50+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+60+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+60+y)*const), stroke_width=8, stroke="brown"))
    dwg.add(dwg.text("PIN7", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+60+y)*const), font_size=16, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-1+400/2+x)*const, (height/2-300/2+1+70+y)*const), end=((width/2+400/2+2+x)*const, (height/2-300/2-2+70+y)*const), stroke_width=8, stroke="orange"))
    dwg.add(dwg.text("PIN8", insert=((width/2+400/2+2+8+x)*const, (height/2-300/2+9-7+70+y)*const), font_size=16, style="text-anchor:center"))


print("len", len(records))


dwg.save()
width = 1200
height = 200
data_list = {}
data_list_u = {}
data_list_d = {}
data_list_bl = {}
data_list_longlength = {}
data_listd = []
data_listu = []

const = 2

dmin = []

for sline in records:
    if sline[0] not in data_list_d:
         data_list_d[sline[0]] = []
    data_list_d[sline[0]].append(sline[2])
    if sline[0] not in data_list_u:
         data_list_u[sline[0]] = []
    data_list_u[sline[0]].append(sline[1])
    if sline[0] not in data_list_bl:
         data_list_bl[sline[0]] = []
    data_list_bl[sline[0]].append(sline[4])
    if sline[0] not in data_list_longlength:
         data_list_longlength[sline[0]] = []
    data_list_longlength[sline[0]].append(sline[5])


for sline in records:
    if sline[2] < 8:
        data_listd.append({"type": sline[0], "down": sline[2]/100, "up": sline[1]/100, })

for sline in records:
    if sline[1] < 8:
        data_listu.append({"type": sline[0], "down": sline[1]/100, "up": sline[1]/100, })

# draw down small loops 1d diagram
# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-down-small-loops-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))
dwg.add(dwg.line(start=(width/2*const, (30)*const), end=(width/2*const, 30*const+7*40), stroke_width=4, stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=(width/2*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(200)/2)*const, (30)*const), end=((width/2+int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{2}", insert=((width/2+int(200)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(400)/2)*const, (30)*const), end=((width/2+int(400)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{4}", insert=((width/2+int(400)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))

for pin in data_list_u:
        data_list_d[pin].sort()
        l5 = int(len(data_list_d[pin])/20)
        l10 = int(len(data_list_d[pin])/10)
        if (pin == "pin1") or (pin == "PIN1"):
            color = "red"
            c = 0
        elif (pin == "pin2") or (pin == "PIN2"):
            color = "green"
            c = 1
        elif (pin == "pin3") or (pin == "PIN3"):
            color = "blue"
            c = 2
        elif (pin == "pin4") or (pin == "PIN4"):
            color = "pink"
            c = 3
        elif (pin == "pin5") or (pin == "PIN5"):
            color = "magenta"
            c = 4
        elif (pin == "pin6") or (pin == "PIN6"):
            color = "cyan"
            c = 5
        elif (pin == "pin7") or (pin == "PIN7"):
            color = "brown"
            c = 6
        elif (pin == "pin8") or (pin == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        dwg.add(dwg.line(start=((width/2+int(data_list_d[pin][0])/2-3)*const, (30)*const+c*40), end=((width/2+int(data_list_d[pin][-1])/2+3)*const, (30)*const+c*40), stroke_width=1, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(data_list_d[pin][l5])/2-3)*const, (30)*const+c*40), end=((width/2+int(data_list_d[pin][-1-l5])/2+3)*const, (30)*const+c*40), stroke_width=4, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(data_list_d[pin][l10])/2-3)*const, (30)*const+c*40), end=((width/2+int(data_list_d[pin][-1-l10])/2+3)*const, (30)*const+c*40), stroke_width=9, stroke=color))
dwg.save()

# draw up smapp loops 1d diagram
# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-up-small-loops-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))
dwg.add(dwg.line(start=((width/2-int(200)/2)*const, (30)*const), end=((width/2-int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-2}", insert=((width/2-int(200)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=(width/2*const, (30)*const), end=(width/2*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=(width/2*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(200)/2)*const, (30)*const), end=((width/2+int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{2}", insert=((width/2+int(200)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(400)/2)*const, (30)*const), end=((width/2+int(400)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{4}", insert=((width/2+int(400)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))

for pin in data_list_u:
        data_list_u[pin].sort()
        l5 = int(len(data_list_u[pin])/20)
        l10 = int(len(data_list_u[pin])/10)

        if (pin == "pin1") or (pin == "PIN1"):
            color = "red"
            c = 0
        elif (pin == "pin2") or (pin == "PIN2"):
            color = "green"
            c = 1
        elif (pin == "pin3") or (pin == "PIN3"):
            color = "blue"
            c = 2
        elif (pin == "pin4") or (pin == "PIN4"):
            color = "pink"
            c = 3
        elif (pin == "pin5") or (pin == "PIN5"):
            color = "magenta"
            c = 4
        elif (pin == "pin6") or (pin == "PIN6"):
            color = "cyan"
            c = 5
        elif (pin == "pin7") or (pin == "PIN7"):
            color = "brown"
            c = 6
        elif (pin == "pin8") or (pin == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        dwg.add(dwg.line(start=((width/2+int(data_list_u[pin][0])/2-3)*const, (30)*const+c*40), end=((width/2+int(data_list_u[pin][-1])/2+3)*const, (30)*const+c*40), stroke_width=1, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(data_list_u[pin][l5])/2-3)*const, (30)*const+c*40), end=((width/2+int(data_list_u[pin][-1-l5])/2+3)*const, (30)*const+c*40), stroke_width=4, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(data_list_u[pin][l10])/2-3)*const, (30)*const+c*40), end=((width/2+int(data_list_u[pin][-1-l10])/2+3)*const, (30)*const+c*40), stroke_width=9, stroke=color))
dwg.save()

# draw big loops 1d diagram
# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-big-loop-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))

dwg.add(dwg.line(start=((width/2-int(400)/2)*const, (30)*const), end=((width/2-int(400)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-4}", insert=((width/2-int(400)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2-int(200)/2)*const, (30)*const), end=((width/2-int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-2}", insert=((width/2-int(200)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=(width/2*const, (30)*const), end=(width/2*const, 30*const+7*40), stroke_width=4, stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=(width/2*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(200)/2)*const, (30)*const), end=((width/2+int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{2}", insert=((width/2+int(200)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(400)/2)*const, (30)*const), end=((width/2+int(400)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{4}", insert=((width/2+int(400)/2)*const, (30)*const), color="grey", font_size=20, style="text-anchor:center"))

for pin in data_list_bl:
    for d in data_list_bl[pin]:
        data_list_bl[pin].sort()
        l5 = int(len(data_list_bl[pin])/20)
        l10 = int(len(data_list_bl[pin])/10)
#        print("pin dd", pin,"l je ", l5, l10, data_list_[pin][0], data_list_d[pin][l5], data_list_d[pin][l10], data_list_d[pin][-1-l10], data_list_d[pin][-1-l5], data_list_d[pin][-1], data_list_d[pin])
        if (pin == "pin1") or (pin == "PIN1"):
            color = "red"
            c = 0
        elif (pin == "pin2") or (pin == "PIN2"):
            color = "green"
            c = 1
        elif (pin == "pin3") or (pin == "PIN3"):
            color = "blue"
            c = 2
        elif (pin == "pin4") or (pin == "PIN4"):
            color = "pink"
            c = 3
        elif (pin == "pin5") or (pin == "PIN5"):
            color = "magenta"
            c = 4
        elif (pin == "pin6") or (pin == "PIN6"):
            color = "cyan"
            c = 5
        elif (pin == "pin7") or (pin == "PIN7"):
            color = "brown"
            c = 6
        elif (pin == "pin8") or (pin == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        dwg.add(dwg.line(start=((width/2+int(d)/2)*const-5, (30)*const+c*40-5), end=((width/2+int(d)/2)*const+5, (30)*const+c*40+5), stroke_width=2, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(d)/2)*const+5, (30)*const+c*40-5), end=((width/2+int(d)/2)*const-5, (30)*const+c*40+5), stroke_width=2, stroke=color))
dwg.save()

### violin graph big loop sum charge
constx = 0.5
consty = 1
width = width
height=height *2+40
k = 3
dwg = svgwrite.Drawing(f"{new_dir}/diagram-big-loop2-charge-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="white"))
print("rozmery ", width, height)
vyska = 50
o=300
dwg.add(dwg.text(f"average charge of central loop", insert=((10+int(100)/k+240)*constx*2-75, (30)*consty-10+vyska-27), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(1800)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(1800)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-18}", insert=((-17+width/k-int(1800)/k)*constx+o-4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(1600)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(1600)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-16}", insert=((-17+width/k-int(1600)/k)*constx+o-4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(1400)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(1400)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-14}", insert=((-17+width/k-int(1400)/k)*constx+o-4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(1200)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(1200)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-12}", insert=((-17+width/k-int(1200)/k)*constx+o-4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(1000)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(1000)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-10}", insert=((-17+width/k-int(1000)/k)*constx+o-4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(800)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(800)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-8}", insert=((-17+width/k-int(800)/k)*constx+o+2, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(600)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(600)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-6}", insert=((-17+width/k-int(600)/k)*constx+o+2, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k-int(400)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(400)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-4}", insert=((-17+width/k-int(400)/k)*constx+o+2, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))

dwg.add(dwg.line(start=((width/k-int(200)/k)*constx+o, (30)*consty+vyska), end=((width/k-int(200)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-2}", insert=((-17+width/k-int(200)/k)*constx+o+2, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))

dwg.add(dwg.line(start=(width/k*constx+o, (30)*consty+vyska), end=(width/k*constx+o, 30*consty+8*40+vyska), stroke_width=4, stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=(width/k*constx+o-3, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))

dwg.add(dwg.line(start=((width/k+int(200)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(200)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{2}", insert=((-17+width/k+int(200)/k)*constx+o+4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(400)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(400)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{4}", insert=((-17+width/k+int(400)/k)*constx+o+4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(600)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(600)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{6}", insert=((-17+width/k+int(600)/k)*constx+o+4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(800)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(800)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{8}", insert=((-17+width/k+int(800)/k)*constx+o+4, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(1000)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(1000)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{10}", insert=((-17+width/k+int(1000)/k)*constx+o, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(1200)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(1200)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{12}", insert=((-17+width/k+int(1200)/k)*constx+o, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(1400)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(1400)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{14}", insert=((-17+width/k+int(1400)/k)*constx+o, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(1600)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(1600)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{16}", insert=((-17+width/k+int(1600)/k)*constx+o, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/k+int(1800)/k)*constx+o, (30)*consty+vyska), end=((width/k+int(1800)/k)*constx+o, 30*consty+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{18}", insert=((-17+width/k+int(1800)/k)*constx+o, (30)*consty-10+vyska), color="grey", font_size=20, style="text-anchor:center"))

r = 80
#const = 1
for pin in data_list_bl:
    r =int(len(data_list_bl[pin])/4)
    data_list_bl[pin].sort()
    ll_min = min(data_list_bl[pin])
    ll_max = max(data_list_bl[pin])
    ll_avg = avg(data_list_bl[pin])
#    print("dlb", ll_min, data_list_bl[pin])
    last = 0
    for i in range(r):
        int_min = ll_min + (ll_max - ll_min)*i/r
        int_max = ll_min + (ll_max - ll_min)*(i+1)/r + 1
#        print(f"int_min ll_min={ll_min}, ll_max={ll_max}, ll_max-ll_min={ll_max-ll_min}, i/r={i/r}, celkem={(ll_max - ll_min)*i/r}")
        counter = 0
        for c in data_list_bl[pin]:
            if int_min <= c and c < int_max:
                counter += 1
#        print("im", int_min, int_max, counter, data_list_bl[pin])
        if (pin == "pin1") or (pin == "PIN1"):
            color = "red"
            c = 0
        elif (pin == "pin2") or (pin == "PIN2"):
            color = "green"
            c = 1
        elif (pin == "pin3") or (pin == "PIN3"):
            color = "blue"
            c = 2
        elif (pin == "pin4") or (pin == "PIN4"):
            color = "pink"
            c = 3
        elif (pin == "pin5") or (pin == "PIN5"):
            color = "magenta"
            c = 4
        elif (pin == "pin6") or (pin == "PIN6"):
            color = "cyan"
            c = 5
        elif (pin == "pin7") or (pin == "PIN7"):
            color = "brown"
            c = 6
        elif (pin == "pin8") or (pin == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        c = c + 0.5
        counter = counter*1.4*r/len(data_list_bl[pin])
#        print("counter", counter,"int_min", int_min, "int_max", int_max, "width", width/k*const, const)
        print("sourad", [int((width/k+int_min/k)*constx+o), int((30-last)*consty+c*40)])
        dwg.add(dwg.polygon([[(width/k+int_min/k)*constx+o, (30-last)*consty+c*40+vyska], [(width/k+int_min/k)*constx+o, (30+last)*consty+c*40+vyska], [(width/k+int_max/k)*constx+o, (30+counter)*consty+c*40+vyska], [(width/k+int_max/k)*constx+o, (30-counter)*consty+c*40+vyska]], fill=color, stroke=color))
        last = counter
    print("stred", [int((width/k+(ll_avg+5)/k)*constx+o), int((30-5)*consty+c*40)])
    dwg.add(dwg.line(start=((width/k+(ll_avg+5)/k)*constx+o, 30*consty+c*40+vyska), end=((width/k+(ll_avg-5)/k)*constx+o, 30*consty+c*40+vyska), fill="black", stroke="black", stroke_width=10))

dwg.add(dwg.text(f"PIN1", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+0*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN2", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+1*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN3", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+2*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN4", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+3*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN5", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+4*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN6", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+5*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN7", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+6*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN8", insert=((-17+width/k-int(800)/k)*constx+o+8-220, int((30-5)*consty+7*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))

dwg.save()
heigh = height /2



# draw big loops 1d diagram
# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-big-loop-size-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))

dwg.add(dwg.line(start=((width/2-int(400)/2)*const, (30)*const), end=((width/2-int(400)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-400}", insert=((width/2-int(400)/2)*const, (30)*const-10), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2-int(200)/2)*const, (30)*const), end=((width/2-int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{-200}", insert=((width/2-int(200)/2)*const, (30)*const-10), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=(width/2*const, (30)*const), end=(width/2*const, 30*const+7*40), stroke_width=4, stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=(width/2*const, (30)*const-10), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(200)/2)*const, (30)*const), end=((width/2+int(200)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{200}", insert=((width/2+int(200)/2)*const, (30)*const-10), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+int(400)/2)*const, (30)*const), end=((width/2+int(400)/2)*const, 30*const+7*40), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{400}", insert=((width/2+int(400)/2)*const, (30)*const-10), color="grey", font_size=20, style="text-anchor:center"))

for pin in data_list_longlength:
    for d in data_list_longlength[pin]:
        if (pin == "pin1") or (pin == "PIN1"):
            color = "red"
            c = 0
        elif (pin == "pin2") or (pin == "PIN2"):
            color = "green"
            c = 1
        elif (pin == "pin3") or (pin == "PIN3"):
            color = "blue"
            c = 2
        elif (pin == "pin4") or (pin == "PIN4"):
            color = "pink"
            c = 3
        elif (pin == "pin5") or (pin == "PIN5"):
            color = "magenta"
            c = 4
        elif (pin == "pin6") or (pin == "PIN6"):
            color = "cyan"
            c = 5
        elif (pin == "pin7") or (pin == "PIN7"):
            color = "brown"
            c = 6
        elif (pin == "pin8") or (pin == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        dwg.add(dwg.line(start=((width/2+int(d)/2)*const-5, (30)*const+c*40-5), end=((width/2+int(d)/2)*const+5, (30)*const+c*40+5), stroke_width=2, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(d)/2)*const+5, (30)*const+c*40-5), end=((width/2+int(d)/2)*const-5, (30)*const+c*40+5), stroke_width=2, stroke=color))
dwg.save()


const = 1
# draw histogram big loop size
dwg = svgwrite.Drawing(f"{new_dir}/diagram-big-loop-size2-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))
k = 1.5
odsaz = 100
vyska =30
dwg.add(dwg.text(f"size of central loop", insert=((10+int(100)/k)*const-5+odsaz, (30)*const-10+vyska-27), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=(10*const+odsaz, (30)*const+vyska), end=(10*const+odsaz, 30*const+8*40+vyska), stroke_width=4, stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=(10*const-5+odsaz, (30)*const-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((10+int(100)/k)*const+odsaz, (30)*const+vyska), end=((10+int(100)/k)*const+odsaz, 30*const+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{100}", insert=((-2+int(100)/k)*const-5+odsaz, (30)*const-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((10+int(200)/k)*const+odsaz, (30)*const+vyska), end=((10+int(200)/k)*const+odsaz, 30*const+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{200}", insert=((-2+int(200)/k)*const-5+odsaz, (30)*const-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((10+int(300)/k)*const+odsaz, (30)*const+vyska), end=((10+int(300)/k)*const+odsaz, 30*const+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{300}", insert=((-2+int(300)/k)*const-5+odsaz, (30)*const-10+vyska), color="grey", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((10+int(400)/k)*const+odsaz, (30)*const+vyska), end=((10+int(400)/k)*const+odsaz, 30*const+8*40+vyska), stroke_width=2, stroke="grey"))
dwg.add(dwg.text(f"{400}", insert=((-2+int(400)/k)*const-5+odsaz, (30)*const-10+vyska), color="grey", font_size=20, style="text-anchor:center"))

r = 80

for pin in data_list_longlength:
    r =int(len(data_list_longlength[pin])/4)
    data_list_longlength[pin].sort()
    ll_min = min(data_list_longlength[pin])
    ll_max = max(data_list_longlength[pin])
    ll_avg = avg(data_list_longlength[pin])
    print("ll_min", ll_min, ll_avg, ll_max), data_list_longlength[pin]
    last = 0
    for i in range(r):
        int_min = ll_min + (ll_max - ll_min)*i/r
        int_max = ll_min + (ll_max - ll_min)*(i+1)/r + 1
        counter = 0
        for c in data_list_longlength[pin]:
            if int_min <= c and c < int_max:
                counter += 1
        print("im", int_min, int_max, counter, data_list_longlength[pin])
        if (pin == "pin1") or (pin == "PIN1"):
            color = "red"
            c = 0
        elif (pin == "pin2") or (pin == "PIN2"):
            color = "green"
            c = 1
        elif (pin == "pin3") or (pin == "PIN3"):
            color = "blue"
            c = 2
        elif (pin == "pin4") or (pin == "PIN4"):
            color = "pink"
            c = 3
        elif (pin == "pin5") or (pin == "PIN5"):
            color = "magenta"
            c = 4
        elif (pin == "pin6") or (pin == "PIN6"):
            color = "cyan"
            c = 5
        elif (pin == "pin7") or (pin == "PIN7"):
            color = "brown"
            c = 6
        elif (pin == "pin8") or (pin == "PIN8"):
            color = "orange"
            c = 7
        if c not in csum:
            csum[c] = 0
            ccounter[c] = 0
        counter = counter*1*r/len(data_list_longlength[pin])
#        print("counter", counter, const)
#        print("sourad", [(10+int_min)*const, (30-last)*const+c*40], [(10+int_min)*const, (30+last)*const+c*40], [(10+int_max)*const, (30+counter)*const+c*40], [(10+int_max)*const, (30-counter)*const+c*40])
        dwg.add(dwg.polygon([[(10+int_min/k)*const+odsaz, (30-last)*const+c*40+20+vyska], [(10+int_min/k)*const+odsaz, (30+last)*const+c*40+20+vyska], [(10+int_max/k)*const+odsaz, (30+counter)*const+c*40+20+vyska], [(10+int_max/k)*const+odsaz, (30-counter)*const+c*40+20+vyska]], fill=color, stroke=color))
        last = counter
    print("stred", [int((10+(ll_avg+5)/k)*const), int((30-5)*const+c*40)])
    dwg.add(dwg.line(start=((10+(ll_avg+2)/k)*const+odsaz, 30*const+c*40+20+vyska), end=((10+(ll_avg-2)/k)*const+odsaz, 30*const+c*40+20+vyska), fill="black", stroke="black", stroke_width=10))


pom = -160
dwg.add(dwg.text(f"PIN1", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+0*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN2", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+1*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN3", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+2*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN4", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+3*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN5", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+4*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN6", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+5*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN7", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+6*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))
dwg.add(dwg.text(f"PIN8", insert=((-17+width/k-int(800)/k)*const+8-40+pom, int((30-5)*const+7*40+20+vyska) ), color="black", font_size=20, style="text-anchor:center"))


dwg.save()



print("data_list_bl", data_list_bl.keys())
print("DONE")

#df = pd.DataFrame(data_list_bl)
#sns.violinplot(data=df, y="up", x="type", split=True, cut=0, bw=.2)
#plt.show()

