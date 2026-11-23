#!/usr/bin/env python3
# tool which generate for given pin directory its standart graphs with up/down part charges sum

import argparse
import time
import os
import json
import svgwrite


def up_part(ch_data, part1, part2, text):

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
        aux.append(ch_data[str(res)])
        res += 1

    if part1 == "up":
        up = aux
    else:
        down = aux
    aux2 = []
    while res <= res_max:
        aux2.append(ch_data[str(res)])
        res += 1
    if part2 == "up":
        up += aux2
    else:
        down += aux2

    return up, down


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
            ch_data[str(res)] = data_all_res[res]
            res += 1
        return up_part(ch_data, part, part, text)
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
        return up_part(ch_data, part, part, text)
    else:
        return [], []


def dolu(x, y, res_list, sses, opacity, text="", mmax=0, avg=False, protein_name=""):
    res_min = None
    res_max = None
    for i in res_list:
        if res_min is None or res_min > int(i):
            res_min = int(i)
        if res_max is None or res_max < int(i):
            res_max = int(i)


def graph_full(data, sses_order, fn_order, avg, first, data_all_res, mapping_names, name):
    # generate diagram and record

    # sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"]    # small dataset(45 pcs)
    sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]  # BIG DATASET

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 0, 1, "down", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 1, 2, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 2, 3, "down", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 3, 4, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 4, 5, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 5, 6, "down", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 6, 7, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 7, 8, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 8, 9, "down", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 9, 10, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 10, 11, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 11, 12, "up", name)

    r_up, r_down = up_part_end(name, data, sses_list, data_all_res, 12, "up", name)

    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 5, 6, "down", name)
    # print loops only charges
    return [mapping_names[name], int(sum(r_up)), int(sum(r_down)), [mapping_names[name]]]


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

# print("=========================")
# for all PINs generate the output html record and diagram
records = []

for f in fn_order:
    # generate diagram and recors
    r = graph_full(data, sses_order, [f], 2, True, data_all_res[f], pin_mapping, f)
    if r is not False:
        records.append(r)

const = 2
# do diagram
width = 1200
height = 200

# draw axes
dwg = svgwrite.Drawing(f"{new_dir}/diagram-up_down-{name}.svg", size=(width*const, height*const))
dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))

dwg.add(dwg.line(start=((width/2-600/2)*const, (height/2-500)*const), end=((width/2-600/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{-6}", insert=((width/2-600/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2-400/2)*const, (height/2-500)*const), end=((width/2-400/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{-4}", insert=((width/2-400/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2-200/2)*const, (height/2-500)*const), end=((width/2-200/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{-2}", insert=((width/2-200/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2-0/2)*const, (height/2-500)*const), end=((width/2-0/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{0}", insert=((width/2+0/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+200/2)*const, (height/2-500)*const), end=((width/2+200/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{2}", insert=((width/2+200/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+400/2)*const, (height/2-500)*const), end=((width/2+400/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{4}", insert=((width/2+400/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+600/2)*const, (height/2-500)*const), end=((width/2+600/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{6}", insert=((width/2+600/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
dwg.add(dwg.line(start=((width/2+800/2)*const, (height/2-500)*const), end=((width/2+800/2)*const, (height/2+500)*const), stroke="grey"))
dwg.add(dwg.text(f"{8}", insert=((width/2+800/2+3)*const, (height-10)*const), color="white", font_size=20, style="text-anchor:center"))
pcRotate = 'rotate(45,%s, %s)' % ((width/2+600/2+3)*const, (height/2+200/2-3)*const)
pcRotate = 'rotate(45,%s, %s)' % ((width/2+400/2+3)*const, (height/2+400/2-3)*const)

csum = {}
ccounter = {}
for sline in records:
        if (sline[0] == "pin1") or (sline[0] == "PIN1"):
            color = "red"
            c = 0
        elif (sline[0] == "pin2") or (sline[0] == "PIN2"):
            color = "green"
            c = 1
        elif (sline[0] == "pin3") or (sline[0] == "PIN3"):
            color = "blue"
            c = 2
        elif (sline[0] == "pin4") or (sline[0] == "PIN4"):
            color = "pink"
            c = 3
        elif (sline[0] == "pin5") or (sline[0] == "PIN5"):
            color = "magenta"
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
        dwg.add(dwg.line(start=((width/2+int(sline[2])/2-3)*const, (30-int(sline[1])/2-3)*const+c*40), end=((width/2+int(sline[2])/2+3)*const, (30-int(sline[1])/2+5)*const+c*40), stroke_width=2, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(sline[2])/2+3)*const, (30-int(sline[1])/2-3)*const+c*40), end=((width/2+int(sline[2])/2-3)*const, (30-int(sline[1])/2+5)*const+c*40), stroke_width=2, stroke=color))

for i in range(8):
     if i in ccounter:
         dwg.add(dwg.line(start=((width/2+int(csum[i]/ccounter[i])/2-3)*const, (31)*const+i*40), end=((width/2+int(csum[i]/ccounter[i])/2+3)*const, (31)*const+i*40), stroke_width=14, stroke="black"))
         dwg.add(dwg.text(f"PIN{i+1}", insert=((width/2+int(csum[i]/ccounter[i])/2-3)*const+5, (37)*const+i*40+10), color="white", font_size=20, style="text-anchor:center"))

dwg.save()


