#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import layout_loader
import logging
import json

# start logger
logger = logging.getLogger()
logger.handlers=[]
logger.setLevel(logging.DEBUG)
bh = logging.StreamHandler()
bh.setLevel(logging.DEBUG)
logger.addHandler(bh)

# Load params
annot_fn = sys.argv[1]
filename_sses_stats = sys.argv[2]
chain_id = sys.argv[3]
domain = sys.argv[4]
output = sys.argv[5]


if len(sys.argv) > 6:
    start_layout_filename = sys.argv[6]
else:
    start_layout_filename = None

if len(sys.argv) > 7:
    ostats = sys.argv[7]
else:
    ostats = None

logging.debug(f"INFO: Generating 2D layout for {annot_fn} - chain {chain_id} domain {domain}")
logging.debug(f"INFO: Using SSEs stats from {filename_sses_stats}")
logging.debug(f"INFO: Using start layout from {start_layout_filename}")

# Load secondary scructures stats file
sses_stats = layout_loader.SSESStats(filename_sses_stats)

# Load protein annotation
annotation = layout_loader.Annotation(annot_fn, chain_id, domain)

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

if ostats == "None":
    ostats = None

# Generate new layout and save it
layout = layout_loader.Layout(sses_stats=sses_stats, annotation=annotation)

if start_layout is not None:
    layout.compute(chain_id, domain, start_layout=start_layout, start_angles=start_angles)
else:
    layout.compute(chain_id, domain)

layout.save_stats(chain_id, domain, statsfn=ostats)
layout.save(chain_id, domain, output)
logging.debug("INFO: Stats saved into %s" % ostats)
logging.debug("INFO: Layout saved into %s" % output)
