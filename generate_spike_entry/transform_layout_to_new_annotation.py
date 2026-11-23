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
import copy
import json

import utils_stats
import utils_common
import utils_annotation
import utils_transform_layout_sses_names
import utils_color_protein
import utils_convert_c_to_p
import utils_synchronize_rotation_files
import convert_protein_method_chain_to_std

sys.path.append('..')
import utils_generate_svg


# read arguments
parser = argparse.ArgumentParser(
    description='generate single and multiple diagrams for a list of spike protein chains')
parser.add_argument('--layout', required=True, action="store",
                    help="layout which have to renamed to new annotation names")
parser.add_argument('--annotation_l', required=True, action="store",
                    help="new annotation layout")
parser.add_argument('--output', required=True, action="store",
                    help="name of output layout with new annotation names")
args = parser.parse_args()

ret = utils_synchronize_rotation_files.sync_layout_sses_names(args.annotation_l, args.layout, args.output)

