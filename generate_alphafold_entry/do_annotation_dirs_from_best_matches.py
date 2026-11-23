#!/usr/bin/env python3
# tool which generate for given directory containing best_structures
# (https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/#api-SIFTS-GetBestStructures)
# files (e.g. https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/A0A1D8PCG7)
# creates directory af - containing all "relevant" (coverage > 80%) af domain detected files
#                   pdb - containing all "relevant" (coverage > 80%) pdb domain detected files
# if there is a problem in af file the process is stopped, if it is in pdb file - the protein
# is removed from both directories and the process continue

import argparse
import os
import json
import utils_annotation
import utils_pdbe_api
import sys
import shutil


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--organisms_dir', action="store", required=True,
                    help="directory containing all relevant best_matches files")
parser.add_argument('--af_dir', action="store", required=True,
                    help="generate database for alfa fold cif structures")
parser.add_argument('--output_dir', action="store", required=True,
                    help="directory for the output files")

args = parser.parse_args()

os.mkdir(args.output_dir)
os.mkdir(f"{args.output_dir}/af")
os.mkdir(f"{args.output_dir}/pdb")

counter = 0
for org_file in os.listdir(args.organisms_dir):
    with open(f"{args.organisms_dir}/{org_file}", "r") as fd:
        counter += 1
        # print(counter, org_file)
        data = json.load(fd)
        for rec in data[org_file]:
            if rec["coverage"] > 0.8:   # if the coverage is more then 80%
                # add org_file annotation to af
                # we have to detect only unp_start - unp_end part
                if not os.path.isfile(f"{args.output_dir}/af/{org_file}.cif"):
                    os.symlink(f"{args.af_dir}/{org_file[0]}/{org_file[1]}/{org_file[2]}/{org_file[3]}/AF-{org_file}-F1-model_v4.cif", f"{args.output_dir}/af/{org_file}.cif", )
                af_rec = [org_file, "A", f"{rec['unp_start']}:{rec['unp_end']}"]
                ret = utils_annotation.detect_list(f"{args.output_dir}/af/", f"{args.output_dir}/af/", [af_rec], args.output_dir, whole_protein="cath_alpha_fold")

                if ret != 1:
                    sys.exit()

                # add pdb_id_start_end annotation to pdb
                utils_pdbe_api.download_protein_cif(f"{args.output_dir}/pdb/", rec["pdb_id"])
                try:
                    chain_id = utils_pdbe_api.transform_uniprot_chain_id_to_struct_asym_id(rec["pdb_id"], org_file, rec["chain_id"])
                except:
                    continue
                pdbe_rec = [rec["pdb_id"], chain_id, [rec["start"], rec["end"]]]
                ret = utils_annotation.detect_list(f"{args.output_dir}/pdb/", f"{args.output_dir}/pdb/", [pdbe_rec], args.output_dir, whole_protein="pdbe_best_structures")
                if ret != 1:
                    os.remove(f"{args.output_dir}/af/{org_file}-detected.sses.json")
                else:
                    shutil.move(f'{args.output_dir}/pdb/{rec["pdb_id"]}-detected.sses.json', f'{args.output_dir}/pdb/{pdbe_rec[0]}-{pdbe_rec[1]}-{pdbe_rec[2][0]}-{pdbe_rec[2][1]}-detected.sses.json')

                break

