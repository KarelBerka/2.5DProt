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


def generate_whole_protein_from_scratch(cif_file, protein, working_dir, annotation_file):
    layout_file = f"{working_dir}/layout-{protein}.json"
    command = f"./generate_covid_layout.py {annotation_file} ALL ALL {layout_file} None None"
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
parser.add_argument('--pg', required=True, action="store",
                    help="file with list of protein groups")
parser.add_argument('--input_cif', action="store",
                    help="which will be put the result alphafold directory")
parser.add_argument('--start_layout', action="store",
                    help="define start layout")
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
new_dir = f"generated-{args.pg}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)

p_list = []
with open(args.pg, "r") as fd:
    ln = fd.readline()
    while ln:
        p_list.append(ln[:-1])
        ln = fd.readline()

print("p_list", p_list)

# start logger
logger = start_logger(args.debug, args.info, new_dir, args.pg)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# we want to download read relevant cif files
cif_dir = f"{new_dir}/cif"
os.mkdir(cif_dir)
for protein in p_list:
    cif_file = f"{cif_dir}/{protein}.cif"
    if not (os.path.isfile(cif_file)):
        shutil.copyfile(f"{args.input_cif}/{protein}.cif", f"{cif_dir}/{protein}.cif")
        shutil.copyfile(f"{args.input_cif}/{protein}.pdb", f"{cif_dir}/{protein}.pdb")
sys.exit()

# Detect download protein structure
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)
annot_file = f"{annotation_dir}/{args.protein}-detected.sses.json"
if not (os.path.isfile(annot_file)):
    ret = utils_annotation.detect_list(annotation_dir, cif_dir, [args.protein], new_dir, whole_protein=True)
if (ret == 0):
    sys.exit()
last_time = utils_common.time_log("1 detected computed", last_time)



detect_file = f"{annotation_dir}/{args.protein}-detected.sses.json"
#link_json = None
utils_color_protein.switch_color_to_chain(detect_file)

generate_whole_protein_from_scratch(cif_file, args.protein, new_dir, detect_file)
last_time = utils_common.time_log("2 layout ready ", start_time)


# draw single image per chain and whole together
svg_args = argparse.Namespace(
    add_descriptions="1",
    data=False,
    debug=False,
    irotation=None,
    otemplate=None,
    layouts=[f"{new_dir}/layout-{args.protein}.json"],
    opac='0',
    orotation=None,
    itemplate=None,
    output=f"{new_dir}/image-{args.protein}.svg"
)

utils_generate_svg.compute(svg_args)
last_time = utils_common.time_log("4 image ready", last_time)
last_time = utils_common.time_log("5 OK", start_time)

sys.exit()
#for pc in add:
#    svg_args.layouts = [f"{new_dir}/layout-{pc[0]}_{pc[1]}.json"]
#    svg_args.output = f"{new_dir}/image-{pc[0]}_{pc[1]}.svg"
#    utils_generate_svg.compute(svg_args)



# now we have to go through all chains
# found out all domain ranges which are already their diagrams
# and aggregate this intervals
# count layoutn of chain  "A" of imput alpha fold enntity
layouts = []
# set start layout - None, or transformed one
if args.start_layout:
    # here we have to transform start layout to the same annotation
    sl_dir=f"{new_dir}/start_layout"
    os.makedirs(sl_dir)
    sl_filename = args.start_layout[args.start_layout.rfind("/")+1:]
    origsl_filename = f"{sl_dir}/{sl_filename}.orig"
    shutil.copy(args.start_layout, origsl_filename)
    sl = sl_filename.split("-")
    cend = sl[2].find(".")
    pc=(f"{sl[1]}-F1", sl[2][3:cend])
    cif_file = f"{cif_dir}/{pc[0]}.cif"
    generate_cl(cif_file, pc[0], pc[1], new_dir, annotation_dir, sl_dir, "None", link_json)
    start_l = f"{sl_dir}/start_layout.json"
    print(f"utils_transform_layout_sses_names.transform_layout_annotation({origsl_filename}, f{sl_dir}/{sl_filename}, {start_l}")
    utils_transform_layout_sses_names.transform_layout_annotation(origsl_filename, f"{sl_dir}/{sl_filename}", start_l)
else:
    start_l = "None"

#for pc in add:
#    cif_file = f"{cif_dir}/{pc[0]}.cif"
#    layout=f"{new_dir}/layout-{pc[0]}_{pc[1]}.json"
#    generate_cl(cif_file, pc[0], pc[1], new_dir, annotation_dir, new_dir, start_l, link_json)
#    layouts.append(layout)

last_time = utils_common.time_log("3 layout computed", last_time)

full_layout_file = f"{new_dir}/layout-{args.pc}.json"
chain_layout_files = f"{new_dir}/layout-*_*.json"
generate_layout(cif_file, full_layout_file, chain_layout_files, new_dir)


# draw single image per chain and whole together
svg_args = argparse.Namespace(
    add_descriptions="1",
    data=False,
    debug=False,
    irotation=None,
    otemplate=None,
    layouts=[f"{new_dir}/layout-{args.pc}.json"],
    opac='0',
    orotation=None,
    itemplate=None,
    output=f"{new_dir}/image-{args.pc}.svg"
)

utils_generate_svg.compute(svg_args)
for pc in add:
    svg_args.layouts = [f"{new_dir}/layout-{pc[0]}_{pc[1]}.json"]
    svg_args.output = f"{new_dir}/image-{pc[0]}_{pc[1]}.svg"
    utils_generate_svg.compute(svg_args)


svg_args.orotation=f"{new_dir}/rotation.json"
svg_args.otemplate=f"{new_dir}/ctemplate.json"
svg_args.layouts = layouts
svg_args.output = f"{new_dir}/image-multi.svg"
utils_generate_svg.compute(svg_args)

svg_args.orotation=None
svg_args.otemplate=None
svg_args.irotation=f"{new_dir}/rotation.json"
svg_args.itemplate=f"{new_dir}/ctemplate.json"
svg_args.layouts = layouts
svg_args.output = f"{new_dir}/image-multi2.svg"
utils_generate_svg.compute(svg_args)

last_time = utils_common.time_log("4 image created", last_time)


#if (not args.debug):
#    shutil.rmtree(f"{cif_dir}")
#    shutil.rmtree(f"{annotation_dir}")
last_time = utils_common.time_log("5 OK", start_time)
