#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

LIGANDS_TIMESTAMP = "2022-07-12"

def link_ligands(protein):
    return f"https://2dprots.ncbr.muni.cz/pub/residue-summary/{LIGANDS_TIMESTAMP}/{protein[1:3]}/{protein}.json"

def link_pdb(protein):
    return f"https://files.rcsb.org/download/{protein}.pdb"

def link_cif(protein):
    return f"https://files.rcsb.org/download/{protein}.cif"

def link_ligand_info(ligand):
    return f"https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{ligand}.cif"

def link_ligands_rotation(protein, chain, domain):
    return f"https://2dprots.ncbr.muni.cz/pub/domain-rotation/{LIGANDS_TIMESTAMP}/{protein[1:3]}/{protein}{chain}{domain}-rotation.json"

def link_channels(protein):
    return f"https://channelsdb2.biodata.ceitec.cz/api/download/pdb/{protein}/json"
#    return f"https://webchem.ncbr.muni.cz/API/ChannelsDB/PDB/{protein}"  old version

def link_channels_af(protein):
    return f"https://channelsdb2.biodata.ceitec.cz/api/download/alphafill/{protein}/json"
