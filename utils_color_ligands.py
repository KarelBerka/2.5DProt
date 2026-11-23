#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import json
from pathlib import Path
import time


chemical_data_link = "../ligand_colors/chemie/data_chemical_new.json"
chemical_colors_link = "../ligand_colors/chemie/colors_chemical_new.json"
functional_data_link = "../ligand_colors/funkce/data_function_new.json"
functional_colors_link = "../ligand_colors/funkce/colors_function.json"


def set_ligand_chemical_color(ecmp):
    empty = []
    ssum = 0
    t=time.time()
    with open(chemical_data_link, "r") as f:
        data = json.load(f)
        for lig in data:
            ssum += 1
            if (lig["AA_derivatives"] == 0) and (lig["NA_derivatives"]["pyrimidine"] == 0) and (lig["NA_derivatives"]["pyridine"] == 0) and (lig["NA_derivatives"]["purine"] == 0) \
                and (lig["saccharides_derivatives"]["pyran"] == 0) and (lig["saccharides_derivatives"]["furan"] == 0) and (lig["lipids"] == 0) and (lig["metallo-organic"] == 0) \
                and (lig["organic_compounds"]["carbonyl"] == 0) and (lig["organic_compounds"]["ether"] == 0) and (lig["organic_compounds"]["carboxyl"] == 0) \
                and (lig["organic_compounds"]["alcohol"] == 0) and (lig["organic_compounds"]["cyclohexane_derivatives"] == 0):
                    empty.append(lig["id"])
            if lig["id"] == ecmp:
                with open(chemical_colors_link, "r") as f:
                    c_data = json.load(f)
                    if lig["AA_derivatives"]:
                        return c_data["AA_derivatives"]
                    elif lig["NA_derivatives"]["purine"]:
                        return c_data["purine"]
                    elif lig["NA_derivatives"]["pyrimidine"]:
                        return c_data["pyrimidine"]
                    elif lig["NA_derivatives"]["pyridine"]:
                        return c_data["pyridine"]
                    elif lig["saccharides_derivatives"]["pyran"]:
                        return c_data["pyran"]
                    elif lig["saccharides_derivatives"]["furan"]:
                        return c_data["furan"]
                    elif lig["lipids"]:
                        return c_data["lipids"]
                    elif lig["inorganic_compounds"]:
                        return c_data["inorganic_compounds"]
                    elif lig["metallo-organic"]:
                        return c_data["metallo-organic"]
                    elif lig["organic_compounds"]["carboxyl"]:
                        return c_data["carboxyl"]
                    elif lig["organic_compounds"]["carbonyl"]:
                        return c_data["carbonyl"]
                    elif lig["organic_compounds"]["ether"]:
                        return c_data["ether"]
                    elif lig["organic_compounds"]["alcohol"]:
                        return c_data["alcohol"]
                    elif lig["organic_compounds"]["cyclohexane_derivatives"]:
                        return c_data["cyclohexane_derivatives"]
                    break
    return "black"


def set_ligand_chemical_type(ecmp):
    empty = []
    ssum = 0
    t=time.time()
    with open(chemical_data_link, "r") as f:
        data = json.load(f)
        for lig in data:
            ssum += 1
            if (lig["AA_derivatives"] == 0) and (lig["NA_derivatives"]["pyrimidine"] == 0) and (lig["NA_derivatives"]["pyridine"] == 0) and (lig["NA_derivatives"]["purine"] == 0) \
                and (lig["saccharides_derivatives"]["pyran"] == 0) and (lig["saccharides_derivatives"]["furan"] == 0) and (lig["lipids"] == 0) and (lig["metallo-organic"] == 0) \
                and (lig["organic_compounds"]["carbonyl"] == 0) and (lig["organic_compounds"]["ether"] == 0) and (lig["organic_compounds"]["carboxyl"] == 0) \
                and (lig["organic_compounds"]["alcohol"] == 0) and (lig["organic_compounds"]["cyclohexane_derivatives"] == 0):
                    empty.append(lig["id"])
            if lig["id"] == ecmp:
                with open(chemical_data_link, "r") as f:
                    c_data = json.load(f)
                    if lig["AA_derivatives"]:
                        return "AA derivatives"
                    elif lig["NA_derivatives"]["purine"]:
                        return "purine"
                    elif lig["NA_derivatives"]["pyrimidine"]:
                        return "pyrimidine"
                    elif lig["NA_derivatives"]["pyridine"]:
                        return "pyridine"
                    elif lig["saccharides_derivatives"]["pyran"]:
                        return "pyran"
                    elif lig["saccharides_derivatives"]["furan"]:
                        return "furan"
                    elif lig["lipids"]:
                        return "lipids"
                    elif lig["inorganic_compounds"]:
                        return "inorganic compounds"
                    elif lig["metallo-organic"]:
                        return "metallo organic"
                    elif lig["organic_compounds"]["carboxyl"]:
                        return "carboxyl"
                    elif lig["organic_compounds"]["carbonyl"]:
                        return "carbonyl"
                    elif lig["organic_compounds"]["ether"]:
                        return "ether"
                    elif lig["organic_compounds"]["alcohol"]:
                        return "alcohol"
                    elif lig["organic_compounds"]["cyclohexane_derivatives"]:
                        return "cyclohexane derivatives"
                    break
                return "not sorted"
    return "not sorted"


def set_ligand_functional_color(ecmp):
    empty = []
    ssum = 0
    t=time.time()
    with open(functional_data_link, "r") as f:
        data = json.load(f)
        for lig in data:
            ssum += 1
            if (lig["drug"] == 0) and (lig["metabolite"] == 0) and (lig["cofactor"] == 0) and (lig["reactant"] == 0) \
                and (lig["NMR_solvent"] == 0):
                    empty.append(lig["id"])
            if lig["id"] == ecmp:
                with open(functional_colors_link, "r") as f:
                    c_data = json.load(f)
                    if lig["cofactor"]:
                        return c_data["cofactor"]
                    elif lig["drug"]:
                        return c_data["drug"]
                    elif lig["metabolite"]:
                        return c_data["metabolite"]
                    elif lig["reactant"]:
                        return c_data["reactant"]
                    elif lig["NMR_solvent"]:
                        return c_data["NMR_solvent"]
                    break
    return "black"


def set_ligand_functional_type(ecmp):
    empty = []
    ssum = 0
    t=time.time()
    with open(functional_data_link, "r") as f:
        data = json.load(f)
        for lig in data:
            ssum += 1
            if (lig["drug"] == 0) and (lig["metabolite"] == 0) and (lig["cofactor"] == 0) and (lig["reactant"] == 0) \
                and (lig["NMR_solvent"] == 0):
                    empty.append(lig["id"])
            if lig["id"] == ecmp:
                with open("../ligand_colors/example_colors_2.json", "r") as f:
                    c_data = json.load(f)
                    if lig["cofactor"]:
                        return "cofactor"
                    elif lig["drug"]:
                        return "drug"
                    elif lig["metabolite"]:
                        return "metabolite"
                    elif lig["reactant"]:
                        return "reactant"
                    elif lig["NMR_solvent"]:
                        return "NMR solvent"
                    break
                return "not sorted"
    return "not sorted"


def get_rainbow_colors(sses,res):
    if len(sses) <= 0:
        return []
    r = {}
    min = 10000
    max = 0
    for sse in sses:
       if sse in res:
          if res[sse][0] < min:
             min = res[sse][0]
          if res[sse][1] > max:
             max = res[sse][1]
    i = 0
    coef=(max-min)
    if coef == 0:
       coef = 1

    for sse in sses:
       if sse in res:
          i = i + 1
          hr=int((res[sse][0]+res[sse][1])/2)
          color_s='s'+str(int((hr-min)*999/coef)).zfill(3)
          color_c = pymol.cmd.get_color_tuple(color_s)
          color_p = (int(color_c[0]*255), \
                   int(color_c[1]*255), \
                   int(color_c[2]*255))
          color_sse = '#%02x%02x%02x' % color_p
          r[sse]=(color_sse)
       else:
          r[sse]=('black')
    return r

