#!/usr/bin/env python3

import argparse
import time
import utils_molespecypic_output

# read arguments
parser = argparse.ArgumentParser(
    description='create json containing protein tunnels data for molespecypic input')
parser.add_argument('--dir', required=True, action="store",
                    help="directory containing 'annotation' dir with annotation files, '2e' dir with tunnel files, 'cifs' dir with cif files")
parser.add_argument('--output', required=False, action="store",
                    help="output directory name")

args = parser.parse_args()

main_dir = args.dir
dir_split = main_dir.split()
if args.output:
    dir_name = f"{args.output}/{dir_split[-1].split('-')[2]}"
else:
    dir_name = dir_split[-1].split("-")[2]

split_dir_name = dir_name.split("_")
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
utils_molespecypic_output.do_molespecypic_output(main_dir, f"{dir_name}_{timestamp}.json", split_dir_name[1])
