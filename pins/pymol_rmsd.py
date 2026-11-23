#!/usr/bin/env python3
# go through unannotated channels and write which channel is the nearest one

import argparse
import pymol
from pymol.cgo import *
from pymol import cmd
import json
import os


def read_file_list(filename):
    fl = []
    with open(filename, "r") as fd:
        line = fd.readline()
        while line:
            print("line", line)
            sline = line.split(" ")
            fl.append([sline[0], sline[1], sline[2][:-1]])
            line = fd.readline()
    return fl


def generate_axis_set(ch_filename):
    axis = []
    with open(ch_filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.split(" ")
            axis.append([float(sline[0]), float(sline[1]), float(sline[2])])
            line = fd.readline()
    return axis


def generate_cgo_input(filename, radius, color, directory, ligand_rotation):
    s_color = color.split(",")
    print("color", color)
    print("s_color", s_color)
    c = (float(s_color[0]), float(s_color[1]), float(s_color[2]))

    with open(ligand_rotation, "r") as rf:
        r_data = json.load(rf)

    with open(filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.strip().split(",")
            for rec in sline:
                srec = rec.split("'")
                if len(srec) < 2:
                    continue
                fn = srec[1]
                ch_filename = f"{directory}/{fn}"
                axis = generate_axis_set(ch_filename)
                linelist = [BEGIN, LINES, COLOR, c[0], c[1], c[2]]
                for i in axis:
                    fix = list(numpy.array(numpy.dot(numpy.linalg.inv(r_data["rotation"]),
                          numpy.array([i[0], i[1], i[2]]) - numpy.array([r_data["translation"][0][0], r_data["translation"][1][0], r_data["translation"][2][0]]))))
                    linelist.extend([VERTEX, fix[0], fix[1], fix[2]])
                linelist.append(END)
                print("processed file :", fn)
                cmd.load_cgo(linelist, fn[2:])
            line = fd.readline()


# read arguments
parser = argparse.ArgumentParser(
    description='generate statistic regarding channel with channel names from the input directory')
parser.add_argument('--pdb_dir1', required=True, action="store",
                    help="protein domain name, which should be shown with the channels")
parser.add_argument('--pdb_dir2', required=True, action="store",
                    help="protein domain name, which should be shown with the channels")
args = parser.parse_args()

pymol.finish_launching()

csum = 0
ccount = 0
cc1 = 0
cc2 = 0
data = {}
fn_list = []

for filename in os.listdir(f"{args.pdb_dir1}"):
    fn_list.append(filename)
    data[filename] = {}
    print("filename", filename)
    cc1 += 1
    if cc1 == 200:
        break

    cc2 = 0
    fn_list2 = []
    for filename2 in os.listdir(f"{args.pdb_dir2}"):
        cc2 += 1
        fn_list2.append(filename2)

        cmd.load(f"{args.pdb_dir1}/{filename}")
        name1 = os.path.basename(f"{args.pdb_dir1}/{filename}")[:-4]

        cmd.load(f"{args.pdb_dir2}/{filename2}")
        name2 = os.path.basename(f"{args.pdb_dir2}/{filename2}")[:-4]

        if name1 == name2:
            data[filename][filename2] = 0
            continue

        ppp = cmd.align(name1, name2)
        ccount += 1
        csum += ppp[0]
        data[filename][filename2] = ppp[0]
        print(f"rmsd {name1} vs {name2} - {ppp[0]:.1f}")
        cmd.remove(name1)
        cmd.remove(name2)
        if cc2 == 200:
            break

cmd.quit()

print()
print("html table of rmsds ======")
print("<table>")
print("<tr>", end=" ")
for fn in fn_list2:
    print("<td>fn</td>", end=" ")
print("</tr>")

for fn in fn_list:
    print(f"<tr><td>{fn}</td>", end=" ")
    for fn2 in fn_list2:
        if data[fn][fn2] < 1:
            print("<td bgcolor='#00DD00'>%.1f</td>" % data[fn][fn2], end=" ")
        elif data[fn][fn2] < 2:
            print("<td bgcolor='#00BB00'>%.1f</td>" % data[fn][fn2], end=" ")
        elif data[fn][fn2] < 3:
            print("<td bgcolor='#009900'>%.1f</td>" % data[fn][fn2], end=" ")
        elif data[fn][fn2] < 4:
            print("<td bgcolor='#007700'>%.1f</td>" % data[fn][fn2], end=" ")
        elif data[fn][fn2] < 5:
            print("<td bgcolor='#005500'>%.1f</td>" % data[fn][fn2], end=" ")
        else:
            print("<td bgcolor='#FF0000'>%.1f</td>" % data[fn][fn2], end=" ")
    print("</tr>")
print("</table>")
