#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import utils_common
import copy

AF_URL = "https://alphafold.ebi.ac.uk/files"

def download_protein_cif(directory, protein):
    pdbe_protein_url = f"{AF_URL}/AF-{protein}-F1-model_v6.cif"
    dest_file = f"{directory}/{protein}.cif"
    ret = utils_common.get_url_to_file(pdbe_protein_url, dest_file)
    return ret


def download_protein_pdb(directory, protein):
    pdbe_protein_url = f"{AF_URL}/AF-{protein}-F1-model_v6.pdb"
    dest_file = f"{directory}/{protein}.pdb"
    ret = utils_common.get_url_to_file(pdbe_protein_url, dest_file)
    return ret


