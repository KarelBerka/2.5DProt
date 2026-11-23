#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
sys.path.append('..')
import layout_loader
import logging
import loader

def start_logger(debug, info, directory, protein):
    root = logging.getLogger()
    if debug:
        root.setLevel(logging.DEBUG)
    else:
        root.setLevel(logging.INFO)
#    logfile = f"{directory}/main-{protein}.log"
#    fh = logging.FileHandler(logfile)
#    if debug:
#        fh.setLevel(level=logging.DEBUG)
#    else:
#        fh.setLevel(level=logging.INFO)
#    root.addHandler(fh)

    return root

start_logger(True, False, "./", sys.argv[1])

# Load params
#filename_sses_stats = sys.argv[1]
annot_fn = sys.argv[1]
chain_id = sys.argv[2]
domain = sys.argv[3]
output = sys.argv[4]

if (len(sys.argv) > 5) and sys.argv[5] != "None":
    filename_sses_stats = sys.argv[5]
else:
    filename_sses_stats = None

if len(sys.argv) > 6:
    start_layout_filename = sys.argv[6]
else:
    start_layout_filename = None

# Load secondary scructures stats file
print("filename_sses_stats", filename_sses_stats)
if filename_sses_stats != None:
    sses_stats = layout_loader.SSESStats(filename_sses_stats)

# Load protein
annotation = layout_loader.Annotation(annot_fn, chain_id, domain)

# If there is a start layout, load start protein
if start_layout_filename is not None and start_layout_filename != "None":
    try:
        start_protein = loader.Layout(start_layout_filename)
        start_layout = start_protein.data['layout']
        start_angles = start_protein.data['angles']
    except IOError:
        start_layout = None
        start_angles = None
else:
    start_layout = None
    start_angles = None



# Generate new layout and save it
#layout = loader.Layout(sses_stats=sses_stats, protein=protein)
print(filename_sses_stats)
if filename_sses_stats != None:
    layout = layout_loader.Layout(sses_stats=sses_stats, annotation=annotation)
else:
    layout = layout_loader.Layout(annotation=annotation)

if start_layout is not None:
    layout.compute(chain_id, domain, start_layout=start_layout, start_angles=start_angles)
else:
    layout.compute(chain_id, domain)

#layout.save_stats(chain_id, domain, statsfn=ostats)
layout.save(chain_id, domain, output)
