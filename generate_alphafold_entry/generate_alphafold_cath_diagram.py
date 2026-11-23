#!/usr/bin/env python3
# tool which generate for a cath directory with files per organism
# e.g. methanocaldococcus_jannaschii
# and alphafold dbs of proteins
# diagrams for all domains
# e.g.: ./generate_alphafold_cath_diagram.py --cath_dir cath/2022-07-11.14_52_00/by_organism/methanocaldococcus_jannaschii/ --pdb_dir data/Methanocaldococcus/ --debug


import logging
import sys
import argparse
import time
import os
import subprocess
import shutil

import utils_common
import utils_annotation
import utils_alphafold_confidence

sys.path.append('..')
import utils_generate_svg


def start_logger(debug, info, directory, protein):
    root = logging.getLogger()
    if debug:
        root.setLevel(logging.DEBUG)
    else:
        root.setLevel(logging.INFO)
    logfile = f"{directory}/main-{protein}.log"
    fh = logging.FileHandler(logfile)
    if debug:
        fh.setLevel(level=logging.DEBUG)
    else:
        fh.setLevel(level=logging.INFO)
    root.addHandler(fh)

    return root


# def simlink_cath_alphafold_list(cif_dir, cath_dir):
def simlink_cath_alphafold_list(cath_cif_dir, pdb_dir, output_dir):
    af_list = []
    # af_location="https://alphafold.ebi.ac.uk/files"
    # cath location https://liveuclac-my.sharepoint.com/personal/ucbtnb4_ucl_ac_uk/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fucbtnb4%5Fucl%5Fac%5Fuk%2FDocuments%2FCATH%2DAlphaFold%2D21organisms&ga=1
    # alphafold list could be only simlinked, not to have several copies
    found = False
    for filename in os.listdir(cath_cif_dir):
        if filename.startswith("af_"):
            af_split = filename.split("_")
            if len(af_split) != 4:
                print("af_split ma moc casti", af_split)
                logging.error("af_split ma moc casti", af_split)
                sys.exit()
            af = af_split[1]
            for pdb_filename in os.listdir(pdb_dir):
                if (pdb_filename.startswith(f"AF-{af}-")) and (pdb_filename.endswith(f".pdb")):
                    os.symlink(f"../../{pdb_dir}/{pdb_filename}", f"{output_dir}/{filename}.pdb")
                    os.symlink(f"../../{pdb_dir}/{pdb_filename[:-4]}.cif", f"{output_dir}/{filename}.cif")
                    found = True
                    af_list.append([f"{filename}","A",f"{af_split[2]}:{af_split[3]}"])
                    break
            if not found:
                print(f"Alpha fold entry {af} is not found in {cath_cif_dir}")
                logging.error(f"Alpha fold entry {af} is not found in {cath_cif_dir}")
    return af_list


def generate_cl(cif_dir, af_list, chain, working_dir, annotation_dir, layout_dir):
    for pdi in af_list:
        protein, domain, interval = pdi
        print("protein", protein)
        print("domain", domain)
        print("interval", interval)
        layout_file = f"{layout_dir}/layout-{protein}.json"
        annotation_file = f"{annotation_dir}/{protein}-detected.sses.json"
        print("layout_file", layout_file)
        print("annotation", annotation_file)
        command = f"./generate_domain_layout_from_scratch.py {annotation_file} {chain} {layout_file} {protein} alpha_fold"
        logging.debug(f"DEBUG: Generating layout for {protein}")
        logging.info(f"INFO: Command {command}")
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()
        if (r != 0) or (not os.path.isfile(layout_file)):
            log_fn = f"{working_dir}/generate_protein_layout_fs-fail.txt"
            utils_common.error_report(working_dir, command, None, err, log_fn)
            logging.error(f"generate_domain_layout_from_scratch.py failed - report save to {log_fn}")
            sys.exit()
            return False
    return True


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--cath_dir', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--pdb_dir', required=True, action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')
args = parser.parse_args()


# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
full_name = args.cath_dir
while full_name[-1] == "/":
    full_name = full_name[:-1]
name = os.path.basename(full_name)
print("fn", full_name)
print("name", name)
new_dir = f"generated-{name}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)

# start logger
logger = start_logger(args.debug, args.info, new_dir, name)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# we want to download read relevant cif files
cif_dir = f"{new_dir}/cif"
os.mkdir(cif_dir)
af_list = simlink_cath_alphafold_list(args.cath_dir, args.pdb_dir, cif_dir)

# Annotate download protein structure
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)

ret = utils_annotation.detect_list(annotation_dir, cif_dir, af_list, new_dir, whole_protein="cath_alpha_fold")
if (ret == 0):
    sys.exit()
last_time = utils_common.time_log("1 detected computed", last_time)

# now we have to go through all chains
# found out all domain ranges which are already their diagrams
# and aggregate this intervals
# count layoutn of chain  "A" of imput alpha fold enntity
generate_cl(cif_dir, af_list, "A", new_dir, annotation_dir, new_dir)
last_time = utils_common.time_log("2 layout computed", last_time)
for af in af_list:
    r = utils_alphafold_confidence.switch_color_to_confidence(f"{new_dir}/layout-{af[0]}.json", f"{cif_dir}/{af[0]}.pdb")

# draw single image per chain and whole together
for af in af_list:
    svg_args = argparse.Namespace(
        add_descriptions="1",
        data=False,
        debug=False,
        irotation=None,
        otemplate=None,
        layouts=[f"{new_dir}/layout-{af[0]}.json"],
        opac='0',
        orotation=None,
        itemplate=None,
        iligands=None,
        ligands=None,
        oligands=None,
        output=f"{new_dir}/image-{af[0]}.svg"
    )
    utils_generate_svg.compute(svg_args)
last_time = utils_common.time_log("3 image created", last_time)

if (not args.debug):
    shutil.rmtree(f"{cif_dir}")
    shutil.rmtree(f"{annotation_dir}")
last_time = utils_common.time_log("4 OK", start_time)
