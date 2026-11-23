#!/usr/bin/env python3
# tool which generate for given protein its directory with
# * protein images
# and image per chain


import os
import sys
import time
import json
import shutil

import utils_common
import utils_cath
import utils_pdbe_api
import utils_read_cif

import argparse
from pathlib import Path
#import logging
#import shutil
#import subprocess
#import copy

#import utils_af_api
#sys.path.append('../../../generate_protein_entry/')
#import utils_annotation

#import utils_read_pdb
#import utils_read_cif
#import utils_sses_buffer

#import utils_generate_svg
#import utils_links

FAMILY = "1.10.630.10"
#FAMILY = "1.20.1070.10"

def download_cif_from_pdb_to_dir(cif_dir, protein):
    # we want to download read relevant ligand cif files
    ret = utils_pdbe_api.download_protein_cif(cif_dir, protein)
    if ret is None:
        print(f"ERROR: cif file {protein} can't be found")
        sys.exit(1)

##############################################################

# read arguments
parser = argparse.ArgumentParser(
    description='for a HEM protein list generate set of tunnels - using combination of MOLE calls and sort the tunnels to standart tunnel types')
parser.add_argument('--cath', required=False, action="store",
                    help="generate sorted protein list for given protein")
parser.add_argument('--prot', required=False, action="store",
                    help="generate sorted protein list for given protein")
args = parser.parse_args()

# generate new directory - with timestamp in name
# create new family diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
new_dir = "pdb_uniprot_mapping-" + FAMILY + timestamp
os.mkdir(new_dir)
pdb_jsons_dir = f"{new_dir}/pdb_jsons"
os.mkdir(pdb_jsons_dir)


# at first find all domains in cath family FAMILY ----------
# wget last cath dbs there
if args.cath:
    cf_name = Path(args.cath).name
    path = os.path.join(new_dir, cf_name)
    shutil.copyfile(args.cath, path)
    domain_list = set()

    with open(path, 'r') as f:
        for row in f.read().split('\n'):
            if row == '' or row[0] == '#':
                continue
            row_split = row.split(" ")
            if row_split[2] == FAMILY:
                domain_list.add(row_split[0])
if args.prot:
    domain_list = [args.prot]
else:
    domain_list = utils_cath.cath_domains_of_given_family(new_dir, FAMILY, leave=True)

print(f"We have {len(domain_list)} domains of {FAMILY} family")

protein_list = set()
protein_mapping = {}
uniprot_list = set()

transformation = {}     # we have to transform chain_id to struct_asym_id

# find all UniProt proteins belonging to FAMILY -----------
# for each pdb from cath dbs find mapping to UniProt
counter = 0
for domain in domain_list:
    protein_name = domain[0:4]
    if True:
        protein_link = f"https://www.ebi.ac.uk/pdbe/api/mappings/{protein_name}"
        protein_file = f"{pdb_jsons_dir}/{protein_name}.json"
        try:
            if not os.path.isfile(protein_file):
                utils_common.get_url_to_file(protein_link, protein_file)
            with open(protein_file, "r") as fd:
                data= json.load(fd)
            if len(data[protein_name]["UniProt"]) >= 1:
                mapping = list(data[protein_name]["UniProt"].keys())[0]
                # transform domain chain_id to struct_asym_id
                chain_id = domain[4:-2]
                for inst in data[protein_name]["UniProt"][mapping]["mappings"]:
                    if inst["chain_id"] == chain_id:
                        print("inst", protein_name,  inst)
                        transformation[domain] = f"{protein_name}{inst['chain_id']}{domain[-2:]}"
            else:
                mapping = "-"
        except:
            mapping = "?"

        print(counter, protein_name, mapping)
        counter = counter + 1
        protein_list.add(protein_name)
        protein_mapping[protein_name] = mapping
        uniprot_list.add(mapping)
#    if counter > 20:
#        break

# write list of CATH proteins -> UniProt proteins mapping
output_filename=f"{new_dir}/list_all_cath_proteins"
with open(output_filename, "w") as f:
    for rec in protein_mapping:
        f.write(f"{rec} {protein_mapping[rec]}\n")

# find FAMILY UniProts -> domains mapping (domains from cath already have mapping
#uniprot_list=["P14779"]
print(uniprot_list)

# at first find uniprot ->pdb mapping P08684 -> 1tqn, 1w0e, ...
uniprot_txt_dir = f"{new_dir}/uniprot_txt"
os.mkdir(uniprot_txt_dir)
counter = 0
uniprot_mapping = {}
for uniprot in uniprot_list:
    uniprot_link = f"https://rest.uniprot.org/uniprotkb/{uniprot}.txt"
    uniprot_file = f"{uniprot_txt_dir}/{uniprot}.txt"
    try:
        uniprot_mapping[uniprot] = []
        utils_common.get_url_to_file(uniprot_link, uniprot_file)
        print("p_f", uniprot_file)
        with open(uniprot_file, "r") as fd:
            line = fd.readline()
            while line:
                if line.startswith("DR   PDBsum; "):
                     pdb = line[13:17]
                     uniprot_mapping[uniprot].append(pdb.lower())
                line = fd.readline()
    except:
        # no protein in this uniprot group
        uniprot_mapping[uniprot] = []

print("uniprot_mapping", uniprot_mapping)
for i in uniprot_mapping:
    print("::::", i, uniprot_mapping[i])

#print("transformation", transformation)
for i in transformation:
    print("::::", i, transformation[i])
# now we have to find Uniprot ->pdb -> chains as well
cifs_dir = f"{new_dir}/cifs"
output_filename=f"{new_dir}/list_all_domains"
os.mkdir(cifs_dir)
with open(output_filename, "w") as f_all:
  for up in uniprot_mapping:
    # for each uniprot group go through all its proteins
    output_filename=f"{new_dir}/list_{up}"
    with open(output_filename, "w") as f:
        for p in uniprot_mapping[up]:
            # find all domains of this protein
            # for CATH pdbs we have chains already computed
            # so we have to test whether p is not in this mapping
            found = False
            for d in domain_list:
                if d.startswith(p):
                     if d in transformation:
                         print("->", d, "->", transformation[d])
                         f.write(f"{transformation[d]}\n")
                         f_all.write(f"{up} {transformation[d]}\n")
                     else:
                         print("**",d)
                         f.write(f"{d}\n")
                         f_all.write(f"{up} {d}\n")
                     found = True
            if not found:
                # the protein is not parsed by CATH yet, we have to find all chains - one chain one cyp
                # TODO: do not add all chains some of them are ligands only - do not put them
                # to the list - break consensus computation

                download_cif_from_pdb_to_dir(cifs_dir, p)
                cifs_file = f"{cifs_dir}/{p}.cif"
                a_chains = utils_read_cif.find_all_auth_chains(cifs_file, with_helices=True)
                for ch in a_chains:
                    f.write(f"{p}{ch}\n")
                    f_all.write(f"{up} {p}{ch}\n")
