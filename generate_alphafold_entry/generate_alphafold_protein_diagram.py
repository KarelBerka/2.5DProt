#!/usr/bin/env python3
# tool which generate for given protein its directory with
# * protein images
# and image per chain


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
import utils_af_api

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


def generate_cl(cif_file, protein, chain, working_dir, annotation_dir, layout_dir):
    layout_file = f"{layout_dir}/layout-{protein}.json"
    annotation_file = f"{annotation_dir}/{protein}-detected.sses.json"
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
    else:
        return True


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--af', required=True, action="store",
                    help="generate diagram database for given alphafold entry")
parser.add_argument('--input_cif', action="store",
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
new_dir = f"generated-{args.af}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)


# start logger
logger = start_logger(args.debug, args.info, new_dir, args.af)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# we want to download read relevant cif files
if args.input_cif:
    utils_common.symlink_alphafold_list(annotation_dir, args.input_cif, [args.af])
else:
    ret = utils_af_api.download_protein_cif(annotation_dir, args.af)
    if ret is None:
            print(f"ERROR: cif file {protein} can't be found in af dbs")
            sys.exit(1)
    ret = utils_af_api.download_protein_pdb(annotation_dir, args.af)
    if ret is None:
            print(f"ERROR: cif file {protein} can't be found in af dbs")
            sys.exit(1)

cif_file = f"{annotation_dir}/{args.af}.cif"

# Annotate download protein structure
#annotation_dir = f"{new_dir}/annotation"
#os.mkdir(annotation_dir)
ret = utils_annotation.detect_list(annotation_dir, annotation_dir, [args.af], new_dir, whole_protein=True)
print("ret je", ret)
if (ret == 0):
    sys.exit()
last_time = utils_common.time_log("1 detected computed", last_time)

# now we have to go through all chains
# found out all domain ranges which are already their diagrams
# and aggregate this intervals
# count layoutn of chain  "A" of imput alpha fold enntity
cif_file = f"{annotation_dir}/{args.af}.cif"
generate_cl(cif_file, args.af, "A", new_dir, annotation_dir, new_dir)
last_time = utils_common.time_log("2 layout computed", last_time)

r = utils_alphafold_confidence.switch_color_to_confidence(f"{new_dir}/layout-{args.af}.json", f"{annotation_dir}/{args.af}.pdb")

# draw single image per chain and whole together
svg_args = argparse.Namespace(
    add_descriptions="1",
    data=False,
    debug=False,
    irotation=None,
    otemplate=None,
    layouts=[f"{new_dir}/layout-{args.af}.json"],
    opac='0',
    orotation=None,
    itemplate=None,
    iligands=None,
    ligands=None,
    oligands=None,
    transparency=1,
    groups=False,
    output=f"{new_dir}/image-{args.af}.svg"
)
utils_generate_svg.compute(svg_args)
last_time = utils_common.time_log("3 image created", last_time)

if (not args.debug):
    shutil.rmtree(f"{annotation_dir}")
last_time = utils_common.time_log("4 OK", start_time)
