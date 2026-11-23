#!/usr/bin/env python3
# go through unannotated channels and write which channel is the nearest one

import argparse
import pymol
from pymol.cgo import *
from pymol import cmd
import json
import numpy
import sys
sys.path.append("../..")
import utils_links
import utils_common
import utils_pdbe_api


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
                for i in axis: #[-3:-1]:
                    fix = list(numpy.array(numpy.dot(numpy.linalg.inv(r_data["rotation"]),
                          numpy.array([i[0], i[1], i[2]]) - numpy.array([r_data["translation"][0][0], r_data["translation"][1][0], r_data["translation"][2][0]]))))
                    linelist.extend([VERTEX, fix[0], fix[1], fix[2]])
                linelist.append(END)
                cmd.load_cgo(linelist, fn[2:])
            line = fd.readline()


# read arguments
parser = argparse.ArgumentParser(
    description='generate statistic regarding channel with channel names from the input directory')
parser.add_argument('--directory', required=True, action="store",
                    help="channel files directory")
parser.add_argument('--file_list', required=True, action="store",
                    help="file containing protein names and chain numbers of hanis which would be shown")
parser.add_argument('--protein_domain', required=False, action="store",
                    help="protein domain name, which should be shown with the channels")
parser.add_argument('--pd2', required=False, action="store",
                    help="second protein domain name, which should be shown with the channels")
args = parser.parse_args()


pymol.finish_launching()
if args.protein_domain:
    pcd = args.protein_domain
else:
    pcd = None
    counter = 0
    for file_name in os.listdir(args.directory):
        if file_name.startswith("gnuplot-") and file_name.endswith(".txt") and "-sses-" not in file_name:
            s = file_name.split("-")
            pcd = s[1]
            counter += 1
            if counter == 1:
                break
    if pcd is None:
        print("protein_domain have to be set ", file_name)
        sys.exit()


protein = pcd[:4]
chain = pcd[4:-2]
domain = pcd[-2:]

sequence = utils_pdbe_api.get_domain_ranges(protein, chain, domain)
all_seq = ""
for borders in sequence:
    if all_seq == "":
        all_seq = f"{borders[0]}-{borders[1]}"
    else:
        all_seq = f"{all_seq},{borders[0]}-{borders[1]}"
cmd.fetch(protein)
cmd.remove('solvent')
cmd.extract(f"{protein}_{chain}{domain}", f"chain {chain} and resi {all_seq}")


link_ligands_rotation = utils_links.link_ligands_rotation(protein, chain, domain)
ligand_rotation = f"ligands_data-{protein}{chain}{domain}-rotation.json"
ret = utils_common.curl_file(link_ligands_rotation, "./", ligand_rotation, option="--insecure")
if not ret:
    sys.exit()


file_list = read_file_list(args.file_list)
for i in file_list:
    file_name = i[0]
    radius = i[1]
    color = i[2]
    generate_cgo_input(file_name, radius, color, args.directory, ligand_rotation)

