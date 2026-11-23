#!/usr/bin/env python3

import json
import numpy as np


def add_ligand_record(ligands, lig_name, counter, x, y, z, cc, res, chain):
    # auxiliary function for generate_ligand_json
    # add record in generate_ligand_json format
    if lig_name != "HOH":
        # HOH is ommited
        # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
        # HETATM 3835 FE   HEM A   1      17.140   3.115  15.066  1.00 14.14          FE
        # HETATM 8238  S   SO4 A2001      10.885 -15.746 -14.404  1.00 47.84           S  
        d = {"chain": chain,
             "auth_chain": chain,
             "n_atoms": counter,
             "center": [x/counter, y/counter, z/counter],
             "entity": cc,
             "entity_type": "?",
             "entity_comp": lig_name,
             "entity_description": "",
             "res": res
             }
        ligands[cc] = (d)


def filter_empty_strings_from_list(res):
    # remove empty strings from res list of strings
    r = []
    for i in res:
        if i != "":
            r.append(i)
    return r


def find_fe_hem_coords_auth(cif_filename, chain_id=""):
    # find coords of HEM ligand if it is present
    # otherwise it returns 0
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
            if (line.startswith("_atom_site.type_symbol")):
                pos["type"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 6:
                break
    ret = 0
    with open(cif_filename, "r") as fcif:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input HETATMrecord to parts divided by spaces,
            # there can be more spaces between each record
            # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
            # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

            if line.startswith("HETATM"):
                lig_name = ss[pos["name"]]
                lig_chain = ss[pos["chain"]]
                ligand_type = ss[pos["type"]]

                # ligand can have several types of names HEM, HEC, HEB
                if ((chain_id == lig_chain) or (chain_id == "")) and ((lig_name == "HEM") or (lig_name == "HEC") or (lig_name == "HEB") )  and ligand_type == "FE":
                    x = float(ss[pos["x"]])
                    y = float(ss[pos["y"]])
                    z = float(ss[pos["z"]])
                    ret = [x, y, z]
                    break
            line = fcif.readline()

    return ret


def find_all_fe_hem_coords(cif_filename, chain_id=""):
    # find coords of HEM ligand if it is present
    # otherwise it returns 0
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
            if (line.startswith("_atom_site.type_symbol")):
                pos["type"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 5:
                break
    ret = []
    with open(cif_filename, "r") as fcif:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input HETATMrecord to parts divided by spaces,
            # there can be more spaces between each record
            # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
            # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

            if line.startswith("HETATM"):
                lig_name = ss[pos["name"]]
                ligand_type = ss[pos["type"]]

                # ligand can have several types of names HEM, HEC, HEB
                if ((lig_name == "HEM") or (lig_name == "HEC") or (lig_name == "HEB") )  and ligand_type == "FE":
                    x = float(ss[pos["x"]])
                    y = float(ss[pos["y"]])
                    z = float(ss[pos["z"]])
                    ret.append([x, y, z])
            line = fcif.readline()
    return ret


def find_all_ligands_sizes(cif_filename, ommit_list):
    # return size of all ligands
    lig_list = {}
    pos = {}

    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start
            if (line.startswith("_atom_site.label_entity_id")):
                pos["entity"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    hetatm_part = False
    with open(cif_filename, "r") as fcif:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input HETATMrecord to parts divided by spaces,
            # there can be more spaces between each record
            # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
            # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

            if line.startswith("HETATM"):
                hetatm_part = True
                lig_name = ss[pos["name"]]
                lig_chain = ss[pos["chain"]]
                lig_entity = ss[pos["entity"]]

                if lig_name not in ommit_list:
                    lig_id = lig_entity+lig_chain
                    if lig_id not in lig_list:
                        lig_list[lig_id] = {"name": lig_name, "size": 0}
                    lig_list[lig_id]["size"] += 1
            line = fcif.readline()

            if line.startswith("ATOM") and (hetatm_part):
                return lig_list
    return lig_list


def find_all_ligands_coords(cif_filename, chain_id, ommit_list, wanted_list=""):
    # return coords of all ligands
    lig_list = {}
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
            if (line.startswith("_atom_site.label_entity_id")):
                pos["res"] = counter - start
            if (line.startswith("_atom_site.label_asym_id")):
                pos["res_label"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    with open(cif_filename, "r") as fcif:
        line = fcif.readline()

        while line:
            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input HETATMrecord to parts divided by spaces,
            # there can be more spaces between each record
            # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
            # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

            if line.startswith("HETATM"):
                lig_name = ss[pos["name"]]
                lig_chain = ss[pos["chain"]]

                if ((chain_id == lig_chain) or (chain_id == "ALL")) and lig_name not in ommit_list:
                    if (wanted_list != "") and lig_name not in wanted_list:
                        line = fcif.readline()
                        continue
                    lig_number = f"{ss[pos['res']]}_{ss[pos['res_label']]}"
                    try:
                        x = float(ss[pos["x"]])
                    except:
                        x = float(ss[pos["x"]][1:])
                    try:
                        y = float(ss[pos["y"]])
                    except:
                        y = float(ss[pos["y"]][1:])
                    try:
                        z = float(ss[pos["z"]])
                    except:
                        z = float(ss[pos["z"]][1:])
                    if lig_number not in lig_list:
                        lig_list[lig_number] = {"name": lig_name, "coords": []}
                    lig_list[lig_number]["coords"].append([x, y, z])
            line = fcif.readline()

    return lig_list

def find_all_ligands_data(cif_filename, chain_id, ommit_list="", wanted_list=""):
    # return coords of all ligands
    lig_list = {}
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
            if (line.startswith("_atom_site.label_entity_id")):
                pos["res"] = counter - start
            if (line.startswith("_atom_site.label_asym_id")):
                pos["res_label"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    with open(cif_filename, "r") as fcif:
        line = fcif.readline()

        while line:
            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input HETATMrecord to parts divided by spaces,
            # there can be more spaces between each record
            # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
            # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

            if line.startswith("HETATM"):
                lig_name = ss[pos["name"]]
                lig_chain = ss[pos["chain"]]

                if ((chain_id == lig_chain) or (chain_id == "ALL")) and lig_name not in ommit_list:
                    if (wanted_list != "") and lig_name not in wanted_list:
                        line = fcif.readline()
                        continue
                    lig_number = f"{ss[pos['res']]}_{ss[pos['res_label']]}"
                    try:
                        x = float(ss[pos["x"]])
                    except:
                        x = float(ss[pos["x"]][1:])
                    y = float(ss[pos["y"]])
                    z = float(ss[pos["z"]])
                    if lig_number not in lig_list:
                        lig_list[lig_number] = {"name": lig_name, \
                                                "coords": [], \
                                                "chain": ss[pos['res_label']], \
                                                "auth_chain": lig_chain}

                    lig_list[lig_number]["coords"].append([x, y, z])
            line = fcif.readline()

    return lig_list



def remove_ligands(input_cif_filename, output_cif_filename, omit_list = []):
    # return coords of all ligands
    lig_list = {}
    pos = {}

    with open(input_cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    with open(input_cif_filename, "r") as incif, open(output_cif_filename, "w") as outcif:
        line = incif.readline()
        while line:
            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input HETATMrecord to parts divided by spaces,
            # there can be more spaces between each record
            # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
            # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

            if line.startswith("HETATM"):
                lig_name = ss[pos["name"]]
                lig_chain = ss[pos["chain"]]

                if lig_name in omit_list:
                    outcif.write(line)
            else:
                outcif.write(line)
            line = incif.readline()

    return lig_list


def generate_ligand_json(cif_filename, ligand_filename):
    # return ligand record (defined in add_ligand_record function
    ligands = {}

    pos = {}

    # at first read the record positions:
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
            if (line.startswith("_atom_site.label_entity_id")):
                pos["res"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_comp_id")):
                pos["name"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 6:
                break

    with open(cif_filename, "r") as fcif:
        lig_name = ""
        lig_number = -1
        x, y, z, counter = 0, 0, 0, 0
        cc = 1

        line = fcif.readline()
        lig_chain = ""
        lig_number = 0

        while line:
            if line.startswith("HETATM"):
                s = line.split(" ")
                ss = filter_empty_strings_from_list(s)
                # split input HETATMrecord to parts divided by spaces,
                # there can be more spaces between each record
                # line HETATM 3767 C CHA HEM . . . B 2 -15.56 -19.721 -10.646 1 33.1 ? CHA HEM 508 A 1 508 ? ? ? ?
                # line HETATM 3768 C CHB HEM . . . B 2 -16.452 -22.296 -14.515 1 29.45 ? CHB HEM 508 A 1 508 ? ? ? ?

                if lig_number != ss[pos["res"]] or lig_chain != ss[pos["chain"]]:
                    if lig_name != "":
                        # one ligand data is already red, add relevant record
                        add_ligand_record(ligands, lig_name, counter, x, y, z, cc, int(lig_number), lig_chain)

                        x, y, z, counter = 0, 0, 0, 0
                        cc += 1
                    lig_number = ss[pos["res"]]
                    lig_name = ss[pos["name"]]
                    lig_chain = ss[pos["chain"]]

                x += float(ss[pos["x"]])
                y += float(ss[pos["y"]])
                z += float(ss[pos["z"]])

                counter += 1
            line = fcif.readline()

        if lig_name != "":
            # last ligand data are already red, add relevant record
            add_ligand_record(ligands, lig_name, counter, x, y, z, cc, int(lig_number), lig_chain)

    fd = open(ligand_filename, "w")
    json.dump(ligands, fd, separators=(',', ':'))

    return True


def read_all_atoms_coords(cif_filename, resi=None, chain=None):
    # return coords of all atoms

    atom_list = {}
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
            if (line.startswith("_atom_site.label_seq_id")):
                pos["label_res"] = counter - start
            if (line.startswith("_atom_site.label_asym_id")):
                pos["label_chain"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

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
            if line.startswith("ATOM")  and not line.startswith("ATOMS"):
                atom_id = ss[pos["id"]]
                x = float(ss[pos["x"]])
                y = float(ss[pos["y"]])
                z = float(ss[pos["z"]])
                res = int(ss[pos["label_res"]])
                a_chain = ss[pos["label_chain"]]
                if resi == None or (resi[0]<=res and res<=resi[1]):
                  if chain == None or (chain == a_chain):
                    atom_list[atom_id] = [x, y, z]
            line = fcif.readline()

    return atom_list


def read_all_res(cif_filename, resi=None):
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
            if (line.startswith("_atom_site.label_comp_id")):
                pos["label_res_id"] = counter - start
            if (line.startswith("_atom_site.label_seq_id")):
                pos["label_res"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

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
            if line.startswith("ATOM"):
                res_id = ss[pos["label_res_id"]]
                res = int(ss[pos["label_res"]])
                if resi == None or (resi[0]<=res and res<=resi[1]):
                    if res not in res_list:
                        res_list[res] = res_id
            line = fcif.readline()
    return res_list



def find_all_label_chains(cif_filename, with_helices = False):
    # return chain list

    if with_helices:
        # we want to return only chains with nhelices - so find in which chain helices are
        helices = []
        with open(cif_filename, "r") as fcif:
            line = fcif.readline()

        ch_h_list = []
        pos = {}
        with open(cif_filename, "r") as fcif:
            counter = 0
            start = 0
            line = fcif.readline()

            while line:
                if line.startswith("loop_"):
                    start = counter+1
                if (line.startswith("_struct_conf.end_auth_asym_id")):
                    pos["chain"] = counter - start

                counter += 1
                line = fcif.readline()
                if len(pos) == 7:
                    break
        with open(cif_filename, "r") as fcif:
            line = fcif.readline()

            while line:
                s = line.split(" ")
                ss = filter_empty_strings_from_list(s)
                # split input HELX_P record to parts divided by spaces,
                # there can be more spaces between two records
                # HELX_P HELX_P37 AE1 PHE B 251 ? SER B 258 ? PHE B 271 SER B 278 1 ? 8  
                # HELX_P HELX_P38 AE2 SER B 271 ? HIS B 304 ? SER B 291 HIS B 324 1 ? 34 
                if line.startswith("HELX"):
                    chain = ss[pos["chain"]]
                    if chain not in ch_h_list:
                        ch_h_list.append(chain)
                line = fcif.readline()

    ch_list = []
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    ch_list = []
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
            if line.startswith("ATOM"):
                chain = ss[pos["chain"]]
                if chain not in ch_list:
                    if with_helices:
                        # if we need to ouput only chains with helices (omit ligand only chains)
                        if chain in ch_h_list:
                            ch_list.append(chain)
                    else:
                        # we will show all chains even without helices
                        ch_list.append(chain)

            line = fcif.readline()

    return ch_list


def find_all_auth_chains(cif_filename, with_helices = False):
    # return chain list

    if with_helices:
        # we want to return only chains with nhelices - so find in which chain helices are
        helices = []
        with open(cif_filename, "r") as fcif:
            line = fcif.readline()

        ch_h_list = []
        pos = {}
        with open(cif_filename, "r") as fcif:
            counter = 0
            start = 0
            line = fcif.readline()

            while line:
                if line.startswith("loop_"):
                    start = counter+1
                if (line.startswith("_struct_conf.end_auth_asym_id")):
                    pos["chain"] = counter - start

                counter += 1
                line = fcif.readline()
                if len(pos) == 7:
                    break
        with open(cif_filename, "r") as fcif:
            line = fcif.readline()

            while line:
                s = line.split(" ")
                ss = filter_empty_strings_from_list(s)
                # split input HELX_P record to parts divided by spaces,
                # there can be more spaces between two records
                # HELX_P HELX_P37 AE1 PHE B 251 ? SER B 258 ? PHE B 271 SER B 278 1 ? 8  
                # HELX_P HELX_P38 AE2 SER B 271 ? HIS B 304 ? SER B 291 HIS B 324 1 ? 34 
                if line.startswith("HELX"):
                    chain = ss[pos["chain"]]
                    if chain not in ch_h_list:
                        ch_h_list.append(chain)
                line = fcif.readline()

    ch_list = []
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    ch_list = []
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
            if line.startswith("ATOM") and not line.startswith("ATOMS"):
                chain = ss[pos["chain"]]
                if chain not in ch_list:
                    if with_helices:
                        # if we need to ouput only chains with helices (omit ligand only chains)
                        if chain in ch_h_list:
                            ch_list.append(chain)
                    else:
                        # we will show all chains even without helices
                        ch_list.append(chain)

            line = fcif.readline()

    return ch_list


def change_atoms_coords(cif_filename, atom_list, new_filename):
    # change coords of all ATOMS (not ligands)
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

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    with open(cif_filename, "r") as fcif, open(new_filename, "w") as fnew:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input ATOM record to parts divided by spaces,
            # there can be more spaces between each record
            # ATOM 1    N N   . MET A 1 1   ? 39.525  21.464  -49.802 1.0 50.28 ? 1   MET A N   1 A0A017SFB8 UNP 1   M
            # ATOM 2    C CA  . MET A 1 1   ? 40.350  20.496  -49.049 1.0 50.28 ? 1   MET A CA  1 A0A017SFB8 UNP 1   M
            #
            if line.startswith("ATOM"):
                atom_id = ss[pos["id"]]
                if atom_id in atom_list:

                    ll = len(ss[pos["x"]])
                    if ll == 3:
                        line = line.replace(ss[pos["x"]], f"{atom_list[atom_id][0]:3.2f}", 1)
                    elif ll == 4:
                        line = line.replace(ss[pos["x"]], f"{atom_list[atom_id][0]:4.2f}", 1)
                    elif ll == 5:
                        line = line.replace(ss[pos["x"]], f"{atom_list[atom_id][0]:5.2f}", 1)
                    elif ll == 6:
                        line = line.replace(ss[pos["x"]], f"{atom_list[atom_id][0]:6.2f}", 1)
                    elif ll == 7:
                        line = line.replace(ss[pos["x"]], f"{atom_list[atom_id][0]:7.2f}", 1)

                    ll = len(ss[pos["y"]])
                    if ll == 3:
                        line = line.replace(ss[pos["y"]], f"{atom_list[atom_id][1]:3.3f}", 1)
                    elif ll == 4:
                        line = line.replace(ss[pos["y"]], f"{atom_list[atom_id][1]:4.3f}", 1)
                    elif ll == 5:
                        line = line.replace(ss[pos["y"]], f"{atom_list[atom_id][1]:5.3f}", 1)
                    elif ll == 6:
                        line = line.replace(ss[pos["y"]], f"{atom_list[atom_id][1]:6.3f}", 1)
                    elif ll == 7:
                        line = line.replace(ss[pos["y"]], f"{atom_list[atom_id][1]:7.3f}", 1)

                    ll = len(ss[pos["z"]])
                    if ll == 3:
                        line = line.replace(ss[pos["z"]], f"{atom_list[atom_id][2]:5.2f}")
                    elif ll == 4:
                        line = line.replace(ss[pos["z"]], f"{atom_list[atom_id][2]:4.3f}")
                    elif ll == 5:
                        line = line.replace(ss[pos["z"]], f"{atom_list[atom_id][2]:5.3f}")
                    elif ll == 6:
                        line = line.replace(ss[pos["z"]], f"{atom_list[atom_id][2]:6.3f}")
                    elif ll == 7:
                        line = line.replace(ss[pos["z"]], f"{atom_list[atom_id][2]:7.3f}")
            fnew.write(line)
            line = fcif.readline()

    return atom_list


def transform_label_to_auth_asym_id(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_seq_id")):
                pos["label_asym_id"] = counter - start
            if (line.startswith("_atom_site.auth_seq_id")):
                pos["auth_asym_id"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()

    diff = {}
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
            if line.startswith("ATOM"):
              if (chain == ss[pos["chain"]]):
                label_asym_id = ss[pos["label_asym_id"]]
                auth_asym_id = ss[pos["auth_asym_id"]]
#                diff[int(auth_asym_id)] = int(label_asym_id)
                diff[int(label_asym_id)] = int(auth_asym_id)
#                return int(auth_asym_id) - int(label_asym_id)
            line = fcif.readline()
    return diff


def transform_auth_to_label_asym_id(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_seq_id")):
                pos["label_asym_id"] = counter - start
            if (line.startswith("_atom_site.auth_seq_id")):
                pos["auth_asym_id"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()

    diff = {}
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
            if line.startswith("ATOM"):
              if (chain == ss[pos["chain"]]):
                label_asym_id = ss[pos["label_asym_id"]]
                auth_asym_id = ss[pos["auth_asym_id"]]
#                diff[int(auth_asym_id)] = int(label_asym_id)
                diff[int(auth_asym_id)] = int(label_asym_id)
#                return int(auth_asym_id) - int(label_asym_id)
            line = fcif.readline()
    return diff


def transform_auth_asym_id_to_label_asym_id(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_asym_id")):
                pos["label_chain"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["auth_chain"] = counter - start

            counter += 1
            line = fcif.readline()

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
            if line.startswith("ATOM") and not line.startswith("ATOMS") :
              if (chain == ss[pos["auth_chain"]]):
                  return ss[pos["label_chain"]]
            line = fcif.readline()
    return None


def transform_label_asym_id_to_auth_asym_id(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_asym_id")):
                pos["auth_chain"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["label_chain"] = counter - start

            counter += 1
            line = fcif.readline()

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
            if line.startswith("ATOM"):
              if (chain == ss[pos["auth_chain"]]):
                  return ss[pos["label_chain"]]
            line = fcif.readline()
    return None


def transform_auth_id_name_to_label_id_name(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_asym_id")):
                pos["auth_chain"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["label_chain"] = counter - start

            counter += 1
            line = fcif.readline()

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
            if line.startswith("ATOM"):
              if (chain == ss[pos["label_chain"]]):
                  return ss[pos["auth_chain"]]
            line = fcif.readline()
    return None


def find_auth_asym_id_end(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.auth_seq_id")):
                pos["auth_asym_id"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()

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
            if line.startswith("ATOM"):
              if (chain == ss[pos["chain"]]):
                
                auth_asym_id = ss[pos["auth_asym_id"]]
            line = fcif.readline()
    return int(auth_asym_id)


def find_auth_asym_id_len(cif_filename, chain):
    pos = {}
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.auth_seq_id")):
                pos["auth_seq_id"] = counter - start
            if (line.startswith("_atom_site.label_seq_id")):
                pos["label_seq_id"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()

    auth_asym_start = -1
#    auth_asym_end = -1
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
            if line.startswith("ATOM"):
              if (chain == ss[pos["chain"]]):
#                print("line", line[:-1])
                
                if auth_asym_start == -1:
                        # pymol does not prodce auth_seq_id in its save output - here is sufficient to use label as well
                    auth_asym_start = ss[pos["auth_seq_id"]]
                auth_asym_end = ss[pos["auth_seq_id"]]

            line = fcif.readline()
    return int(auth_asym_end) - int(auth_asym_start)


def find_label_asym_id_len(cif_filename, chain):
    pos = {}
    print("cif_filename, chain", cif_filename, chain)
    with open(cif_filename, "r") as fcif:
        counter = 0
        start = 0
        line = fcif.readline()

        while line:
            if line.startswith("loop_"):
                start = counter+1
            if (line.startswith("_atom_site.label_seq_id")):
                pos["label_asym_id"] = counter - start
            if (line.startswith("_atom_site.label_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()

    label_asym_start = -1
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
            if line.startswith("ATOM"):
              if (chain == ss[pos["chain"]]):
                if label_asym_start == -1:
                    print("line", line[:-1])

                    label_asym_start = ss[pos["label_asym_id"]]
                label_asym_end = ss[pos["label_asym_id"]]
            line = fcif.readline()
    return int(label_asym_end) - int(label_asym_start)


def change_atoms_coords_with_rotation_file_label(cif_filename, rotation_filename, new_filename, chain=None, all_ligands = False):
    # change coords of all ATOMS (not ligands)
    # read rotation files - to rotate ligands
    with open(f"{rotation_filename}", "r") as fr:
        r_data = json.load(fr)
        r = np.linalg.inv(r_data["rotation"])
        t = r_data["translation"]

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
            if (line.startswith("_atom_site.label_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    with open(cif_filename, "r") as fcif, open(new_filename, "w") as fnew:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input ATOM record to parts divided by spaces,
            # there can be more spaces between each record
            # ATOM 1    N N   . MET A 1 1   ? 39.525  21.464  -49.802 1.0 50.28 ? 1   MET A N   1 A0A017SFB8 UNP 1   M
            # ATOM 2    C CA  . MET A 1 1   ? 40.350  20.496  -49.049 1.0 50.28 ? 1   MET A CA  1 A0A017SFB8 UNP 1   M
            #
            if (line.startswith("ATOM") or line.startswith("HETATM")) and ("MISSING" not in line):
                atom_id = ss[pos["id"]]
                x = float(ss[pos["x"]])
                y = float(ss[pos["y"]])
                z = float(ss[pos["z"]])
                sv = [x,y,z]
                atom_list = [r[0][0]*sv[0]+r[1][0]*sv[1]+r[2][0]*sv[2]+t[0][0], \
                     r[0][1]*sv[0]+r[1][1]*sv[1]+r[2][1]*sv[2]+t[1][0], \
                     r[0][2]*sv[0]+r[1][2]*sv[1]+r[2][2]*sv[2]+t[2][0]]
                ll = len(ss[pos["x"]])
                line = line.replace(ss[pos["x"]], f"{str(atom_list[0])[:ll]}", 1)
                ll = len(ss[pos["y"]])
                line = line.replace(ss[pos["y"]], f"{str(atom_list[1])[:ll]}", 1)
                ll = len(ss[pos["z"]])
                line = line.replace(ss[pos["z"]], f"{str(atom_list[2])[:ll]}", 1)
                if chain== None or (chain == ss[pos["chain"]]):
                    # change coords
                    fnew.write(line)
                elif (all_ligands) and line.startswith("HETATM"):
                    # sometimes we want to show all ligands
                    fnew.write(line)
            else:
                # just rewrite the line
                fnew.write(line)
            line = fcif.readline()

    return atom_list


def change_atoms_coords_with_rotation_file(cif_filename, rotation_filename, new_filename, chain=None, all_ligands = False):
    # change coords of all ATOMS (not ligands)
    # read rotation files - to rotate ligands
    with open(f"{rotation_filename}", "r") as fr:
        r_data = json.load(fr)
        r = np.linalg.inv(r_data["rotation"])
        t = r_data["translation"]

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
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    with open(cif_filename, "r") as fcif, open(new_filename, "w") as fnew:
        line = fcif.readline()

        while line:

            s = line.split(" ")
            ss = filter_empty_strings_from_list(s)
            # split input ATOM record to parts divided by spaces,
            # there can be more spaces between each record
            # ATOM 1    N N   . MET A 1 1   ? 39.525  21.464  -49.802 1.0 50.28 ? 1   MET A N   1 A0A017SFB8 UNP 1   M
            # ATOM 2    C CA  . MET A 1 1   ? 40.350  20.496  -49.049 1.0 50.28 ? 1   MET A CA  1 A0A017SFB8 UNP 1   M
            #
            if (line.startswith("ATOM") or line.startswith("HETATM")) and ("MISSING" not in line):
                atom_id = ss[pos["id"]]
                x = float(ss[pos["x"]])
                y = float(ss[pos["y"]])
                z = float(ss[pos["z"]])
                sv = [x,y,z]
                atom_list = [r[0][0]*sv[0]+r[1][0]*sv[1]+r[2][0]*sv[2]+t[0][0], \
                     r[0][1]*sv[0]+r[1][1]*sv[1]+r[2][1]*sv[2]+t[1][0], \
                     r[0][2]*sv[0]+r[1][2]*sv[1]+r[2][2]*sv[2]+t[2][0]]
                ll = len(ss[pos["x"]])
                line = line.replace(ss[pos["x"]], f"{str(atom_list[0])[:ll]}", 1)
                ll = len(ss[pos["y"]])
                line = line.replace(ss[pos["y"]], f"{str(atom_list[1])[:ll]}", 1)
                ll = len(ss[pos["z"]])
                line = line.replace(ss[pos["z"]], f"{str(atom_list[2])[:ll]}", 1)
                if chain== None or (chain == ss[pos["chain"]]):
                    # change coords
                    fnew.write(line)
                elif (all_ligands) and line.startswith("HETATM"):
                    # sometimes we want to show all ligands
                    fnew.write(line)
            else:
                # just rewrite the line
                fnew.write(line)
            line = fcif.readline()

    return atom_list


def find_res_avg_coord(cif_filename, chain, res):
    # change coords of all ATOMS (not ligands)
    # read rotation files - to rotate ligands
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
#            if (line.startswith("_atom_site.auth_seq_id")):
#                pos["res"] = counter - start
#            if (line.startswith("_atom_site.auth_asym_id")):
#                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_seq_id")):
                pos["res"] = counter - start
            if (line.startswith("_atom_site.label_asym_id")):
                pos["chain"] = counter - start

            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    sum_x, sum_y, sum_z, counter = 0, 0, 0, 0
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
            if line.startswith("ATOM") and not line.startswith("ATOM O1 OF") and not line.startswith("ATOMS"):
                atom_id = ss[pos["id"]]
                x = float(ss[pos["x"]])
                y = float(ss[pos["y"]])
                z = float(ss[pos["z"]])
                l_chain = ss[pos["chain"]]
                l_res = int(ss[pos["res"]])
                if l_chain == chain and l_res == res:
                    sum_x += x
                    sum_y += y
                    sum_z += z
                    counter += 1
            line = fcif.readline()
    return [sum_x/counter, sum_y/counter, sum_z/counter]


def find_res_avg_coord_auth(cif_filename, chain, res):
    # change coords of all ATOMS (not ligands)
    # read rotation files - to rotate ligands
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
            if (line.startswith("_atom_site.auth_seq_id")):
                pos["res"] = counter - start
            if (line.startswith("_atom_site.auth_asym_id")):
                pos["chain"] = counter - start
            if (line.startswith("_atom_site.label_seq_id")):
                pos["l_res"] = counter - start
            counter += 1
            line = fcif.readline()
            if len(pos) == 7:
                break

    sum_x, sum_y, sum_z, counter = 0, 0, 0, 0
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
            if line.startswith("ATOM") and not line.startswith("ATOM O1 OF"):
                atom_id = ss[pos["id"]]
                x = float(ss[pos["x"]])
                y = float(ss[pos["y"]])
                z = float(ss[pos["z"]])
                a_chain = ss[pos["chain"]]
                a_res = int(ss[pos["res"]])
                l_res = int(ss[pos["l_res"]])
                if a_chain == chain and l_res == res:
                    sum_x += x
                    sum_y += y
                    sum_z += z
                    counter += 1
            line = fcif.readline()
    return [sum_x/counter, sum_y/counter, sum_z/counter]


