#!/usr/bin/env python3
# tool which generates diagram like:
# http://147.251.21.23/pub/2023_12_05_piny/up_and_down/diagram-up2.png
import argparse
import os
import json

const = 2

# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--dir', required=True, action="store",
                    help="list of all layouts'")
parser.add_argument('--pdbs', required=True, action="store",
                    help="list of pdbs files of all layouts (generated*/cif)")
parser.add_argument('--sse', required=True, action="store",
                    help="sse name")
args = parser.parse_args()


aa_list = {}
bbb = {}

for fl in os.listdir(args.dir):
    with open(f"{args.dir}/{fl}", "r") as fd:
        print("fll", fl)
        data = json.load(fd)
        if (args.sse) in data["helices_residues"]:
            se = data["helices_residues"][args.sse]

    sfl = fl.split("-")
    name = sfl[1]+"-F1.pdb"
    aa_list[sfl[1]] = {}

    with open(f"{args.pdbs}/{name}") as fd:
        line = fd.readline()
        while (line):
            if line.startswith("ATOM"):
                res = int(line[22:26])
                if ((res > int(se[0]) and res < int(se[1]))) or ((res > int(se[1])) and (res < int(se[0]))):
                    if res not in bbb:
                        bbb[res] = {}
                    amin = line[17:20]
                    if amin == "ARG":
                        amin_sc = "R"
                    elif amin == "HIS":
                        amin_sc = "H"
                    elif amin == "LYS":
                        amin_sc = "K"
                    elif amin == "ASP":
                        amin_sc = "D"
                    elif amin == "GLU":
                        amin_sc = "G"
                    elif amin == "SER":
                        amin_sc = "S"
                    elif amin == "THR":
                        amin_sc = "T"
                    elif amin == "ASN":
                        amin_sc = "N"
                    elif amin == "GLN":
                        amin_sc = "Q"
                    elif amin == "CYS":
                        amin_sc = "C"
                    elif amin == "SEC":
                        amin_sc = "U"
                    elif amin == "GLY":
                        amin_sc = "G"
                    elif amin == "PRO":
                        amin_sc = "P"
                    elif amin == "ALA":
                        amin_sc = "A"
                    elif amin == "VAL":
                        amin_sc = "V"
                    elif amin == "ILE":
                        amin_sc = "I"
                    elif amin == "LEU":
                        amin_sc = "L"
                    elif amin == "MET":
                        amin_sc = "M"
                    elif amin == "PHE":
                        amin_sc = "F"
                    elif amin == "TYR":
                        amin_sc = "Y"
                    elif amin == "TRP":
                        amin_sc = "W"

                    else:
                        amin_sc = amin
                    if res not in aa_list[sfl[1]]:
                        if amin_sc not in bbb[res]:
                            bbb[res][amin_sc] = 1
                        else:
                            bbb[res][amin_sc] += 1
                    aa_list[sfl[1]][res] = amin_sc

            line = fd.readline()

for p in aa_list:
    for i in aa_list[p]:
        print(f"{aa_list[p][i]}", end="")
    print()

