#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import json
import argparse
import sys
sys.path.append("..")
import layout_loader_ligands
import layout_loader_channels

# input is list of chain layouts

# Load params - cif file of whole protein domain and list of chains layouts
parser = argparse.ArgumentParser(
    description='generate diagram for a protein structure')
parser.add_argument('--output', required=True, action="store",
                    help="output file name")
parser.add_argument('--protein', required=True, action="store",
                    help="TODO")
parser.add_argument('--channels', action="store",
                    help="channel file name")
parser.add_argument('--sl_layout', action="store",
                    help="layout file name which will be used as SSEs and ligands source")
parser.add_argument('--annotation_dir', action="store",
                    help="annotation directory - channel files should be there")

args = parser.parse_args()

output = args.output
protein = args.protein
channels_fn = args.channels
annot_dir = args.annotation_dir

# read relevant data to per domain structure
with open(args.sl_layout, "r") as fd:
    data = json.load(fd)
    layout_3d = data["3dlayout"]
    input_sizes = data["size"]

print("INFO: Generating 2D layout for %s" % protein)

l = layout_loader_ligands.Layout()
l.read_sl_layout_from_json(args.sl_layout)

if channels_fn != None:
    channels = layout_loader_channels.Channels(channels_fn, "ALL", mol_online=True)
    l.compute_channels(channels, annot_dir, layout_3d=layout_3d, input_sizes=input_sizes, flag="add_channels") # we do not split the whole protein to parts,
                                                                                # all protein is taken as a one bunch (domain)
else:
    channels = None

# save the result layout
l.add_channels_to_protein(args.sl_layout, output, channels)
print("INFO: Protein layout saved into %s" % output)
