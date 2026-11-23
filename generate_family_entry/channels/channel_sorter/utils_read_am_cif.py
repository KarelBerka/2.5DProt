#!/usr/bin/env python3

import json
import numpy as np

def filter_empty_strings_from_list(res):
    # remove empty strings from res list of strings
    r = []
    for i in res:
        if i != "":
            r.append(i)
    return r


def read_all_res_am(cif_filename):
    # return coords of all atoms

    res_list = {}
    pos = {}

    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.Cartn_x")):
                pos["x"] = counter - start
            if (line.startswith("_atom_site.Cartn_y")):
                pos["y"] = counter - start
            if (line.startswith("_atom_site.Cartn_z")):
                pos["z"] = counter - start
            if (line.startswith("_atom_site.id")):
                pos["id"] = counter - start
            if (line.startswith("_atom_site.B_iso_or_equiv")):
                pos["am"] = counter - start
            if (line.startswith("_atom_site.label_seq_id")):
                pos["res"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

#    print("pos", pos)
    with open(cif_filename, "r") as fcif:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input ATOM record to parts divided by spaces,
            # there can be more spaces between each record
            # ATOM 1    N N   . MET A 1 1   ? 39.525  21.464  -49.802 1.0 50.28 ? 1   MET A N   1 A0A017SFB8 UNP 1   M
            # ATOM 2    C CA  . MET A 1 1   ? 40.350  20.496  -49.049 1.0 50.28 ? 1   MET A CA  1 A0A017SFB8 UNP 1   M
            #
#            print("ss", ss)
            if line.startswith("ATOM"):
                atom_id = ss[pos["id"]]
                x = float(ss[pos["x"]])
                y = float(ss[pos["y"]])
                z = float(ss[pos["z"]])
                am = float(ss[pos["am"]])
                res = int(ss[pos["res"]])
                res_list[res] = am
            line = fcif.readline()

    return res_list


