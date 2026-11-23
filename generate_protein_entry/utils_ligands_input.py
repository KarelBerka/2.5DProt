#!/usr/bin/env python3

import utils_read_cif
import json
import sys

sys.path.append("..")
import utils_links
import utils_common


def avg_3d_position(point_list):
    x, y, z = 0, 0, 0
    for point in point_list:
        x += point[0]
        y += point[1]
        z += point[2]
    counter = len(point_list)
    return [x/counter, y/counter, z/counter]


def find_description_and_type_from_ligand_cif(cif_ligand_file):
    ed_found = False
    type_found = False
    with open(cif_ligand_file, "r") as fr:
        line = fr.readline()
        while line:
            if line.startswith("_chem_comp.name "):
                ed_found = True
                entity_description = line[16:-1]
                # remove obsolete " or ' signs from start and end of string
                if (entity_description.startswith('"')) or (entity_description.startswith("'")):
                    entity_description = entity_description[1:]
                if (entity_description.endswith('"')) or (entity_description.endswith("'")):
                    entity_description = entity_description[:-1]
            if line.startswith("_chem_comp.type "):
                l_type = line[16:-1]
                type_found = True
            line = fr.readline()
        if not ed_found:
            entity_description = ""
        if not type_found:
            l_type = "UNKNOWN"
    return entity_description, l_type


def create_ligands_input(cif_file, annotation_dir, output_file, coords=False):
    # create ligands input file
    # ligands_data-{prot}.json format
    # "A":{"chain":"A",                                   *
    #      "auth_chain":"A",                              *
    #      "n_atoms":1297,                                *
    #      "center":[34.983,11.433,9.605],                *
    #      "entity":"1",                                  *
    #      "entity_type":"polymer",                       *
    #      "entity_polymer_type":"polypeptide(L)",        *
    #      "entity_comp":null,                            *
    #      "entity_description":"T4 LYSOZYME",            --
    #      "residues":null}                               *

    # read all ligands which are in protein cif
    lig_data = utils_read_cif.find_all_ligands_data(cif_file, "ALL", "HOH", "")
    lig_transf_data = {}
    for lig in lig_data:
        # for all ligands parse its data and transform it to ligands_data-{prot}.json format
        cif_ligand_link = utils_links.link_ligand_info(lig_data[lig]["name"])
        cif_ligand_file = f"{annotation_dir}/{lig_data[lig]['name']}.cif"
        # we have to download {ligand}.cif to find out ligand description and type
        ret = utils_common.download_file(cif_ligand_link, annotation_dir, cif_ligand_file)
        if not ret:
             # {ligand}.cif is not present so ommit this columns
            entity_description = ""
            l_type = ""
        else:
            # find out _chem_comp.name and type columns to fill entity description and type
            entity_description, l_type = find_description_and_type_from_ligand_cif(cif_ligand_file)
        # if entity description contains RNA or DNA - just name the ligand this way
        if "RNA" in entity_description:
            name = "RNA"
        elif "DNA" in entity_description:
            name = "DNA"
        else:
            name = lig_data[lig]["name"]
        # add whole ligand record to ligands list
        ligand = {"chain": lig_data[lig]["chain"],
                  "auth_chain": lig_data[lig]["auth_chain"],
                  "n_atoms": len(lig_data[lig]["coords"]),
                  "center": avg_3d_position(lig_data[lig]["coords"]),
                  "entity": "1",                         # not filled
                  "entity_type": l_type,
                  "entity_polymer_type": "",             # not filled
                  "entity_comp": name,
                  "entity_description": entity_description,
                  "residues": None}                      # not filled
        if coords:
            ligand["coords"] = lig_data[lig]["coords"]
        lig_transf_data[lig_data[lig]["chain"]] = ligand

    # dump the result ligand list to output file
    with open(output_file, "w") as fw:
        json.dump(lig_transf_data, fw)
    return
