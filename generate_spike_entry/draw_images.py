#!/usr/bin/env python3
# tool which generate for given protein its directory with
# * protein images
# and image per chain


import logging
import sys
import argparse
import time
import os

import utils_common
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


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a spike protein as a whole')
parser.add_argument('--input_layouts', action="store",
                    help="which will be put the result alphafold directory")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
old_dir = args.input_layouts
new_dir = f"{args.input_layouts}/{timestamp}"
os.makedirs(new_dir)

# start logger
logger = start_logger(args.debug, args.info, new_dir, "aux")
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# draw single image per chain and whole together
svg_args = argparse.Namespace(
    add_descriptions="0",
    data=False,
    debug=False,
    irotation=None, #f"{old_dir}/rotation.json",
    itemplate=None,
    layouts=[f"{old_dir}"],
    opac='1',
    orotation=f"{new_dir}/rotation.json",
    otemplate=f"{new_dir}/ctemplate.json",
    output=f"{new_dir}/image-multi1.svg",
    ligands=None,
)
utils_generate_svg.compute(svg_args)
utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png",{new_dir})

#svg_args.irotation=f"{new_dir}/rotation.json"
svg_args.itemplate = f"{new_dir}/ctemplate.json"
svg_args.orotation = None
svg_args.output = f"{new_dir}/image-multi3.svg"
utils_generate_svg.compute(svg_args)

svg_args.opac="0"
svg_args.output=f"{new_dir}/image-multi2.svg"
utils_generate_svg.compute(svg_args)

for fn in os.listdir(old_dir):
    if fn.startswith("layout-") and fn.endswith(".json"):
        svg_args.layouts = [f"{old_dir}/{fn}"]
        svg_args.output = f"{new_dir}/{fn[:-5]}.svg"
        svg_args.itemplate = None
        svg_args.add_descriptions="1"
        utils_generate_svg.compute(svg_args)

last_time = utils_common.time_log("4 image ready", last_time)
last_time = utils_common.time_log("5 OK", start_time)

