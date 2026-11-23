#!/usr/bin/env python3
# tool which generate for given pin directory its standart graphs with up/down part charges sum

import argparse
import time
import os
import json
import svgwrite
import sys


HEADER_TEMPLATE = """
<html>
    <head>
        <title>2DProt status page</title>
        <style>
            table, td, th {
                border: 1px solid black;
            }

            table {
                border-collapse: collapse;
                width: 90%;
            }

            td {
                padding: .5em;
            }
        </style>
    </head>
    <body>
    <table>
"""


def draw_avg_sse_diagram(pin, fn_order, int_set, name, data_all_res, length):
    # print sse
    width = 20
    if length == "short":
        height = 70
    elif length == "half":
        height = 110
    else:
        height = 220
    res_max = 0
    lengths = []
    for fn in int_set[pin]:
        if (len(int_set[pin][fn]) > 1):
            print("do not use graph option for multii intervals")
            sys.exit()
        start, end = int_set[pin][fn][0]
        lengths.append(end-start+1)
    lengths.sort()
    res_max = lengths[int(len(lengths)/10*9.5)]
    print("lenths je", lengths)
    print("res_max je", res_max)

    opacity = 0.2
    ratio = (height-20)/res_max

    dwg = svgwrite.Drawing(f"{new_dir}/{pin}/charges-{pin}-{name}-avg.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
    for fn in int_set[pin]:
            start, end = int_set[pin][fn][0]
            counter = start
            while counter <= end:
                if data_all_res[fn][counter] < -40:
                    color = "red"
                elif data_all_res[fn][counter] > 40:
                    color = "blue"
                else:
                    color = "white"
                shapes = dwg.add(dwg.g(id='shapes', fill=color))
                shapes.add(dwg.rect(insert=(5, 10+ratio*(counter-start)), size=(10, ratio), fill=color, stroke=color, style=f"opacity:{opacity}"))
                counter += 1
    dwg.save()


def record(pin_order, fn_order, int_set, data_all_res, output_file_name, header=False, graph=False, main_title="", length="short"):
    if header:
        st = "w"
    else:
        st = "a"
    with open(output_file_name, st) as fd:
        # at first do a header with pin names
        if header:
            fd.write(HEADER_TEMPLATE)
        # the number of pin domains in the given section
        fd.write(f"    <tr><td rowspan=4>{main_title}</td>\n")
        fd.write("    <td>pin name, number of domains</td>\n")
        for pin in pin_order:
            fd.write(f"        <td><strong>{pin}</strong> <small>({len(fn_order[pin])} dom.)</small></td>\n")
        fd.write("    </tr>\n")

        # length of the  domains in the given section
        fd.write("    <tr>\n")
        fd.write("    <td>avg length of loop</td>\n")
        for pin in pin_order:
            os.makedirs(f"{new_dir}/{pin}", exist_ok=True)
            len_set = []
            for fn in fn_order[pin]:
                ll = 0
                if fn in int_set[pin]:
                    for interval in int_set[pin][fn]:
                        ll += int(interval[1]) - int(interval[0]) + 1
                len_set.append(ll)
            len_set.sort()
            fd.write(f"        <td><small>{sum(len_set)/len(len_set):.1f} res.</small></td>\n")
        fd.write("    </tr>\n")

        # average charge
        fd.write("    <tr>\n")
        fd.write("    <td>avg charge of loop </td>\n")
        for pin in pin_order:
            print("record pro pin", pin, len(int_set[pin]))
            ch_set = []
            ch_set2 = []
            for fn in int_set[pin]:
                print("pin", pin, "file", fn, "interval", int_set[pin][fn], end=": ")
                av_ch = 0
                for interval in int_set[pin][fn]:
                    start, end = interval
                    counter = start
                    while counter <= end:
                        av_ch += data_all_res[fn][counter]
                        print(int(data_all_res[fn][counter]), end=", ")
                        counter += 1
                ch_set.append(av_ch)
                ch_set2.append(round(av_ch/100))
                print(" - ", int(av_ch/100))
                ch_set.sort()
            fd.write(f"        <td><strong>{sum(ch_set)/(len(ch_set)*100):.1f}<strong></td>\n")
        fd.write("    </tr>\n")

        # draw avg graph
        fd.write("    <tr>\n")
        fd.write("    <td>graph</td>\n")
        for pin in pin_order:
            if graph:
                fd.write(f'        <td bgcolor="grey">start<br><img src="{pin}/charges-{pin}-{graph}-avg.svg"><br>end</td>\n')
                draw_avg_sse_diagram(pin, fn_order, int_set, graph, data_all_res, length)
        fd.write("    </tr>\n")

        # finish record table
        fd.write("</body>\n")
        fd.write("</html>\n")


def mmax(slist):
    mmax = 0
    for s in slist:
        if int(s) > mmax:
            mmax = int(s)
    return mmax


def mmin(slist):
    mmin = 10000
    for s in slist:
        if int(s) < mmin:
            mmin = int(s)
    return mmin


def record_main_loop(pin_order, fn_order, sses_order, data_all_res, data, output_file_name):
 #  sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"] # small dataset
    sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]   # BIG DATASET

    # between H9 and H19
    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[0]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, header=True, graph="SSE1", main_title="SSE1", length="long")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[1]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE2", main_title="SSE2", length="long")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[2]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE3", main_title="SSE3", length="long")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[3]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE4a", main_title="SSE4a", length="half")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[4]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE4b", main_title="SSE4b", length="half")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[5]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE5", main_title="SSE5", length="long")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[6]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE6a", main_title="SSE6a", length="half")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[7]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE6b", main_title="SSE6b", length="half")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[7]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE7", main_title="SSE7", length="long")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[8]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE8", main_title="SSE8", length="long")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[9]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE9a", main_title="SSE9a", length="half")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[10]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE9b", main_title="SSE9b", length="half")

    int_set = {}    # set of intervals which are interesting for us
    for pin in pin_order:
        if pin not in int_set:
            int_set[pin] = {}
        for f_name in fn_order[pin]:
            k = data[sses_list[11]][pin][f_name]
            int_set[pin][f_name] = [(min(k), max(k))]
    record(pin_order, fn_order, int_set, data_all_res, output_file_name, graph="SSE10", main_title="SSE10", length="long")


# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains subdirs with protein layouts (e.g. subdir pin1, pins,...)")
parser.add_argument('--charges_dir', required=True, action="store",
                    help="charges data dir")
parser.add_argument('--org_name_dir', required=False, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument("--name", help="name")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = os.path.basename(os.path.normpath(args.directory))
new_dir = f"charges-{name}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)
print("new_dir", new_dir)
output_file_name = f"{new_dir}/index.html"

data = {}          # charges data for pins->files->sses
aux_data = {}      # organism name
pin_order = []     # order of pin list
fn_order = {}      # order of protein list in pin directories
data_all_res = {}  # charges of all residues in proteins (not only sses)
for pin_dir in os.listdir(args.directory):
    pin_order.append(pin_dir)
    fn_order[pin_dir] = []
    for filename in os.listdir(f"{args.directory}/{pin_dir}"):
        # go through all PINs in PIN directory and create html record and diagram for them
        print("filename", pin_dir, filename)

        # read PIN json and read residues charges in SSEs
        with open(f"{args.directory}/{pin_dir}/{filename}", "r") as fd:
            jdata = json.load(fd)
            fn_order[pin_dir].append(filename)
            for sse in jdata["h_type"]:
                if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                    if sse not in data:
                        data[sse] = {}
                    if pin_dir not in data[sse]:
                        data[sse][pin_dir] = {}
                    aux = {}
                    for i in jdata["charges"][sse]:
                        aux[int(i)] = jdata["charges"][sse][i]

                    data[sse][pin_dir][filename] = aux
#                    print("data_all_res22", jdata["charges"][sse])


        # and (if org_name_dir is set) organism name as well
        if args.org_name_dir:
            sfilename = filename.split("-")
            with open(f"{args.org_name_dir}/{sfilename[1]}-F1.cif") as fd:
                line = fd.readline()
                while line:
                    if line.startswith("_ma_target_ref_db_details.organism_scientific"):
                        aux_data[filename] = line[55:]
                        break
                    line = fd.readline()
        else:
            aux_data[filename] = "--"

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
                        print("data_all_res", data_all_res[filename][int(res)])


aux_data_order = sorted(aux_data.items(), key=lambda x: x[1])
sses_order = list(data.keys())

print("=========================")
record_main_loop(pin_order, fn_order, sses_order, data_all_res, data, output_file_name)
