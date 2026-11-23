#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import os.path
import json
sys.path.append("..")
import layout_loader_ligands
import layout_loader_channels

# Load params
annot_fn = sys.argv[1]
chain_id = sys.argv[2]
output = sys.argv[3]
protein_name = sys.argv[4]

if (len(sys.argv) > 5) and sys.argv[5] == "alpha_fold":
    alpha_fold = True
else:
    alpha_fold = False

if (len(sys.argv) > 6):
    start_layout_filename = sys.argv[6]
else:
    start_layout_filename = None

if (len(sys.argv) >7):
    channels_fn = sys.argv[7]
else:
    channels_fn = None


# If there is a start layout, load start protein
if start_layout_filename is not None and start_layout_filename != "None":
    try:
        with open(start_layout_filename, "r") as fd:
            data = json.load(fd)
        start_layout = data['layout']
        start_angles = data['angles']
    except IOError:
        start_layout = {}
        start_angles = {}
else:
    start_layout = {}
    start_angles = {}

print("INFO: Generating 2D layout for %s" % (protein_name))

# Load protein
annotation = layout_loader_ligands.Annotation(annot_fn, chain_id, "ALL")
layout_3d = annotation.layout_3d()

if channels_fn:
    channels = layout_loader_channels.Channels(channels_fn, "ALL")
else:
    channels = None

# Generate new layout and save it
layout = layout_loader_ligands.Layout(sses_stats=None, annotation=annotation)
layout.compute(chain_id, "ALL", os.path.dirname(annot_fn), start_layout=start_layout, start_angles=start_angles, alpha_fold=alpha_fold, )

if channels is not None and len(channels.names) > 0:
    layout.compute_channels(channels, "", layout_3d=layout_3d, flag="domain")

layout.save_stats(chain_id, "AlL", statsfn=None)
layout.save(chain_id, "ALL", output, channels=channels)
print("INFO: Layout saved into %s" % output)
