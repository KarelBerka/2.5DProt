#!/usr/bin/env python3
import json


def is_pdb(pdb_filename):
    # test the first line of the file.
    # pdb file starts with keyword HEADER (cif file do not have it)
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        if line.startswith("HEADER"):
            return True
        else:
            return False


def find_res_coords(pdb_filename, res_number, chain_id):
    # for the given res find its average cooordinates
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        a_counter = 0
        a_x = 0
        a_y = 0
        a_z = 0
        while line:
            if line.startswith("ATOM "):
                # read ATOM tag
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                # ATOM     32  N  AARG A  -3      11.281  86.699  94.383  0.50 35.88           N  
                res = int(line[22:26])
                if res == res_number and chain_id == line[21]:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    a_counter += 1
                    a_x += x
                    a_y += y
                    a_z += z
            line = fpdb.readline()
    return a_x/a_counter, a_y/a_counter, a_z/a_counter


def find_res_atom(pdb_filename):
    # find coordinates of all atoms of all residues in pdb file
    alist = []
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        while line:
            if line.startswith("ATOM "):
                #  read ATOM tag
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                # ATOM     32  N  AARG A  -3      11.281  86.699  94.383  0.50 35.88           N  
                    res = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    alist.append([x, y, z, res])
            line = fpdb.readline()
    return alist


def find_all_res_coords(pdb_filename, chain_id):
    # find average cooordinates for all residues
    res_coords = {}
    res_list = {}
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        while line:
            if line.startswith("ATOM ") and chain_id == line[21]:
                res = int(line[22:26])
                if res not in res_list:
                    res_list[res] = line[17:20]
            line = fpdb.readline()

    for res in res_list:
        res_coords[res] = {"coords": find_res_coords(pdb_filename, res, chain_id), "aa": res_list[res]}
    return res_coords


def find_all_ligands_coords(pdb_filename, chain_id, ommit_list):
    # find average cooordinates for all ligands
    lig_list = {}
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        while line:
            if line.startswith("HETATM") and chain_id == line[21] and line[17:20] not in ommit_list:
                lig_number = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if lig_number not in lig_list:
                    lig_list[lig_number] = {"name": line[16:20], "coords": []}
                lig_list[lig_number]["coords"].append([x, y, z])
            line = fpdb.readline()

    return lig_list


def find_fe_hem_coords(pdb_filename, chain_id):
    # find average cooordinates for HEM if it is present otherwise it returns 0
    ret = 0
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        while line:
            if line.startswith("HETATM") and chain_id == line[21] and line[12:14] == "FE":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ret = [x, y, z]
                break
            line = fpdb.readline()
    return ret


def detect_sses(pdb_filename, output_filename, protein):
    # read pdb file - detect sses
    # read helices and sheets information
    # and create *detected file (like SecStrAnnotator output)
    # pdb file documentation
    # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    sses = []
    chains = set()
    beta_connectivity = []
    counter = 0
    start_st = ""
    prev_st = ""
    with open(pdb_filename, "r") as fpdb:
        line = fpdb.readline()
        while line:
            if line.startswith("HELIX"):
                #  read helix tag - example
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                # HELIX    1  HA GLY A   86  GLY A   94  1                                   9   
                # HELIX    1 AA1 MET A  132  LYS A  139  1
                # HELIX   10  10 LEU A  172  THR A  187  1
                sses.append({"label": "H"+str(int(line[7:10])),
                                "chain_id": line[19],    # there is only 1 digit to chain label -> problem for big proteins
                                "start": int(line[21:25]),
                                "end": int(line[33:37])})
                if line[19] not in chains:
                    chains.add(line[19])
                counter += 1

            if line.startswith("SHEET"):
                #  read sheet tag - example
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                # SHEET    1   A 5 THR A 107  ARG A 110  0
                # SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
                name = "E"+str(int(counter))
                sheet_c = int(line[7:10])
                orientation = int(line[38:40])
                start_st = name
                start_bc = orientation

                if prev_st != "":
                    if sheet_c == 1:
                        if start_st != "" and start_bc != 0:
                            beta_connectivity.append([prev_st, start_st, orientation])
                    else:
                        if orientation != 0:
                            beta_connectivity.append([prev_st, name, orientation])
                prev_st = name
                sses.append({"label": name,
                                "chain_id": line[21],    # there is only 1 digit to chain label -> problem for big proteins
                                "start": int(line[22:26]),
                                "end": int(line[33:37])})
                if line[21] not in chains:
                    chains.add(line[21])
                counter += 1

            line = fpdb.readline()

    if start_st != "" and start_bc != 0:
        # if the last sheet have connected the first strand with the last
        # mark the connection to beta connectivity tag
        beta_connectivity.append([prev_st, start_st, orientation])

    for s in sses:
        s["start_vector"] = find_res_coords(pdb_filename, s["start"], s["chain_id"])
        s["end_vector"] = find_res_coords(pdb_filename, s["end"], s["chain_id"])

    aux = {"secondary_structure_elements": sses,
           "beta_connectivity": beta_connectivity}
    data = {protein: aux}

    fd = open(output_filename, "w")
    json.dump(data, fd, separators=(',', ':'))

    return chains


def add_ligand_record(ligands, lig_name, counter, x, y, z, cc, res, chain):
    # do 2DProt input format ligand redord
    # TOTO: where is this function used
    if lig_name != "HOH":
        #  HOH is ommited
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
             "res": res,
             }
        ligands[cc] = (d)


def generate_ligand_json(pdb_filename, ligand_filename):
    # read ligand data from pdb file and generate ligand json
    ligands = {}
    with open(pdb_filename, "r") as fpdb:
        lig_name = ""
        lig_number = -1
        lig_chain = ""
        x, y, z, counter = 0, 0, 0, 0
        cc = 1

        line = fpdb.readline()
        while line:
            if line.startswith("HETATM"):
                #  read HETATOM tag
                # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
                # HETATM 8237 MG    MG A1001      13.872  -2.555 -29.045  1.00 27.36          MG 
                if lig_number != line[22:26] or lig_chain != line[21]:
                    if lig_name != "":
                        # one ligand data are already red, add relevant record
                        add_ligand_record(ligands, lig_name, counter, x, y, z, cc, int(lig_number), lig_chain)

                        x, y, z, counter = 0, 0, 0, 0
                        cc += 1
                    lig_number = line[22:26]
                    lig_name = line[17:20]
                    lig_chain = line[21]

                x += float(line[30:38])
                y += float(line[38:46])
                z += float(line[46:54])
                counter += 1
            line = fpdb.readline()

        if lig_name != "":
            # last ligand data are already red, add relevant record
            add_ligand_record(ligands, lig_name, counter, x, y, z, cc, int(lig_number), lig_chain)

    fd = open(ligand_filename, "w")
    json.dump(ligands, fd, separators=(',', ':'))

    return
