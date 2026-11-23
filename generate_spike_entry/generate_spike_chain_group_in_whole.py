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

import utils_stats
import utils_common
import utils_annotation
import utils_transform_layout_sses_names
import utils_color_protein

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


def generate_layout(cif_file, chain_layout_file, chain_domain_layout_files, working_dir):
    command = f"./generate_chain_layout.py {cif_file} {chain_layout_file} {chain_domain_layout_files}"
    logging.debug(f"DEBUG: generating chain layout {chain_layout_file}")
    logging.debug(f"DEBUG: command {command}")
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0) or (not os.path.isfile(chain_layout_file)):
        log_fn = f"{working_dir}/generate_chain_layout_fs-fail.txt"
        utils_common.error_report(working_dir, command, None, err, log_fn)
        logging.error(f"generate_domain_layout.py failed - report save to {log_fn}")
        sys.exit()
        return False
    else: 
        return True


def generate_whole_protein_from_scratch(cif_file, protein, working_dir, annotation_file, stats_file):
    layout_file = f"{working_dir}/layout-{protein[0]}{protein[1]}.json"
    command = f"./generate_covid_layout.py {annotation_file} ALL ALL {layout_file} {stats_file} None"
    print("command", command)
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
    description='generate diagram for a spike protein as a whole')
parser.add_argument('--pc_list', required=True, action="store",
                    help="generate diagram database for given protein-chain entry")
parser.add_argument('--input_cif', action="store",
                    help="which will be put the result alphafold directory")
parser.add_argument('--start_layout', action="store",
                    help="define start layout")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('--consensus', action="store",
                    help="consensus files location")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')

args = parser.parse_args()


# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
new_dir = f"generated-{args.pc_list}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)

# start logger
logger = start_logger(args.debug, args.info, new_dir, args.pc_list)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

p_list = []
pc_list = []
with open(args.pc_list) as fd:
    ln = fd.readline()
    while ln:
        s = str(ln).split()
        p_list.append(s[0])
        pc_list.append((s[0],s[1],"ALL", f"{s[0]}{s[1]}00"))
        ln = fd.readline()

# we want to download read relevant cif files
cif_dir = f"{new_dir}/annotation"
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)
for p in p_list:
    cif_file = f"{cif_dir}/{p}.cif"
#    print("pc", pc)
#    if not (os.path.isfile(cif_file)):
    shutil.copyfile(f"{args.input_cif}/{p}.cif", cif_file)
    shutil.copyfile(f"{args.input_cif}/{p}.pdb", f"{cif_dir}/{p}.pdb")


# Detect download protein structure
#for pc in pc_list:
shutil.copyfile(f"{args.consensus}/consensus.cif", f"{annotation_dir}/consensus.cif")
shutil.copyfile(f"{args.consensus}/consensus.sses.json", f"{annotation_dir}/consensus-template.sses.json")
shutil.copyfile(f"{args.consensus}/diagram.json", f"{annotation_dir}/diagram.json")

ret = utils_annotation.annotate_with_given_template_list(annotation_dir, pc_list, new_dir, domain_boundaries = False)
#ret = utils_annotation.detect_chain_list(annotation_dir, cif_dir, pc_list, new_dir)
if (ret == 0):
        sys.exit()
last_time = utils_common.time_log("1 detected computed", last_time)

diagram_file = f"{annotation_dir}/diagram.json"

logging.info("INFO: Starting save statistics and counting average protein")
stats_file = new_dir + "/sses-stats.json"
avg_protein = new_dir + "/avg-protein.json"

ret = utils_stats.copy_statistic(annotation_dir, diagram_file, stats_file)
if ret != []:
    logging.critical(f"save_statistics critically failed - report save to {log_fn}")
    sys.exit()

#detect_file = f"{annotation_dir}/{args.protein}-detected.sses.json"
#utils_color_protein.switch_color_to_chain(detect_file)

for p in pc_list:
    generate_whole_protein_from_scratch(cif_file, p, new_dir, f"{annotation_dir}/{p[0]}{p[1]}00-annotated.sses.json", stats_file)

last_time = utils_common.time_log("2 layout ready ", start_time)


# draw single image per chain and whole together
svg_args = argparse.Namespace(
    add_descriptions="0",
    data=True,
    debug=False,
    irotation=None,
    otemplate=None,
    layouts=[f"{new_dir}"],
    opac='0',
    orotation=None,
    itemplate=None,
    output=f"{new_dir}/image-{args.pc_list}.svg"
)

utils_generate_svg.compute(svg_args)
svg_args.add_descriptions="1"
for pc in pc_list:
    svg_args.layouts = [f"{new_dir}/layout-{pc[0]}{pc[1]}.json"]
    svg_args.output = f"{new_dir}/image-{pc[0]}{pc[1]}.svg"
    utils_generate_svg.compute(svg_args)


last_time = utils_common.time_log("4 image ready", last_time)
last_time = utils_common.time_log("5 OK", start_time)

sys.exit()
