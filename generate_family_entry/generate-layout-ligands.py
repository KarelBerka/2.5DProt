#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse
import logging
import json
import sys
import os.path
sys.path.append("..")
import layout_loader_ligands
import layout_loader_channels


# Load params - cif file of whole protein domain and list of chains layouts
parser = argparse.ArgumentParser(
    description='generate diagram for a protein structure')
parser.add_argument('--annotation_file', required=True, action="store",
                    help="annotation file")
parser.add_argument('--stats_file', required=True, action="store",
                    help="annotation file")
parser.add_argument('--domain', required=True, action="store",
                    help="TODO")
parser.add_argument('--chain', required=True, action="store",
                    help="TODO")
parser.add_argument('--output', required=True, action="store",
                    help="output file name")
parser.add_argument('--start_layout', required=False, action="store",
                    help="output file name")
parser.add_argument('--ostats', required=False, action="store",
                    help="TODO")
parser.add_argument('--channels', action="store",
                    help="channel file name")
parser.add_argument('--ligands', action="store",
                    help="ligand file name")
args = parser.parse_args()

# Load params
annot_fn = args.annotation_file
filename_sses_stats = args.stats_file
domain_id = args.domain
chain_id = args.chain
output = args.output

if args.start_layout:
    start_layout_filename = args.start_layout
else:
    start_layout_filename = None

if args.ostats and (args.ostats != "None"):
    ostats = args.ostats
else:
    ostats = None

if args.ligands:
    ligands = layout_loader_ligands.Ligands(args.ligands, chain_id)
else:
    ligands = None

if args.channels:
    channels = layout_loader_channels.Channels(args.channels, "ALL")
else:
    channels = None

# start logger
logger = logging.getLogger()
logger.handlers = []
logger.setLevel(logging.DEBUG)
bh = logging.StreamHandler()
bh.setLevel(logging.DEBUG)
logger.addHandler(bh)

logging.debug(f"INFO: Generating 2D layout for {annot_fn} - chain {chain_id} domain {domain_id}")
logging.debug(f"INFO: Using SSEs stats from {filename_sses_stats}")
logging.debug(f"INFO: Using start layout from {start_layout_filename}")

# Load secondary scructures stats file
sses_stats = layout_loader_ligands.SSESStats(filename_sses_stats)

# Load protein annotation
# TODO: chain_id domain_id are switched
annotation = layout_loader_ligands.Annotation(annot_fn, domain_id, chain_id)
layout_3d = annotation.layout_3d()

# If there is a start layout, load start protein
if start_layout_filename is not None and start_layout_filename != "None":
    try:
        with open(start_layout_filename, "r") as fd:
            data = json.load(fd)
        start_layout = data['layout']
        start_angles = data['angles']
    except IOError:
        start_layout = None
        start_angles = None
else:
    start_layout = None
    start_angles = None

# Generate new layout and save it
layout = layout_loader_ligands.Layout(sses_stats=sses_stats, annotation=annotation)

if start_layout is not None:
    layout.compute(domain_id, chain_id, os.path.dirname(annot_fn), start_layout=start_layout, start_angles=start_angles)
else:
    layout.compute(domain_id, chain_id, os.path.dirname(annot_fn))

layout.save_stats(domain_id, chain_id, statsfn=ostats)
logging.debug("INFO: Stats saved into %s" % ostats)

if ligands is not None and len(ligands.names) > 0:
    layout.compute_ligands(ligands, near=True)

if channels is not None and len(channels.names) > 0:
    layout.compute_channels(channels, os.path.dirname(annot_fn), layout_3d=layout_3d, flag="domain")
layout.save(chain_id, domain_id, output, ligands=ligands, channels=channels)
logging.debug("INFO: Layout saved into %s" % output)
