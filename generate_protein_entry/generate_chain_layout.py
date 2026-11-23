#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import json
import argparse
import os.path
import sys
sys.path.append("..")
import layout_loader_ligands
import layout_loader_channels

# input is list of chain layouts

# Load params - cif file of whole protein domain and list of chains layouts
parser = argparse.ArgumentParser(
    description='generate diagram for a protein structure')
parser.add_argument('--annotation_file', required=True, action="store",
                    help="annotation file")
parser.add_argument('--output', required=True, action="store",
                    help="output file name")
parser.add_argument('--protein', required=True, action="store",
                    help="TODO")
parser.add_argument('--chains', default=[''], type=str, nargs="*", action="store",
                    help="TODO")
parser.add_argument('--channels', action="store",
                    help="channel file name")
parser.add_argument('--ligands', action="store",
                    help="ligand file name")
parser.add_argument('--mol_online', action="store_true",
                    help="input format is mol_online like")


args = parser.parse_args()

annot_fn = args.annotation_file
output = args.output
protein = args.protein

ligands_fn = args.ligands
channels_fn = args.channels

cd_list = []
layout = {}
a_layout = {}
h_type = {}
bcc = {}
layout_3d = {}
hres = {}
chain_id = {}

i = 6
# read protein lists - read layouts and a_layouts per domain, chain
for pcd_full in args.chains:
    i = i + 1
    if "*" not in pcd_full:
        # the set of domains of the given protein is nonempty
        # set protein, domain, chain triple
        pcd = pcd_full[pcd_full.rfind("/")+1:]

        pcd_split = pcd.split("-")
        pcd_split = pcd_split[1].split(".")
        pcd_split = pcd_split[0].split("_")
        protein = pcd_split[0]
        domain = pcd_split[1][-2:]
        chain = pcd_split[1][:-2]
        if (len(chain) == 0):
            chain = domain
            domain = "ALL"
        if (len(domain) == 0):
            domain = "ALL"
        if protein == protein:
            cd_list.append((chain, domain))
        else:
            print("DEBUG: layout", pcd_full, protein, "does not match to protein", pp)

        # read relevant data to per domain structure
        with open(pcd_full, "r") as fd:
            data = json.load(fd)
        if chain not in layout:
            layout[chain] = {}
            a_layout[chain] = {}
            h_type[chain] = {}
            bcc[chain] = {}
            layout_3d[chain] = {}
            hres[chain] = {}
            chain_id[chain] = {}
        layout[chain][domain] = data['layout']
        a_layout[chain][domain] = data['angles']
        h_type[chain][domain] = data['h_type']
        bcc[chain][domain] = data['beta_connectivity']
        layout_3d[chain][domain] = data['3dlayout']
        hres[chain][domain] = data['helices_residues']
        chain_id[chain][domain] = data['chain_id']
    else:
        # there is no domain in the given protein (contains only set of ligands) e.g. 137d protein
        cd_list = []

layouts = {"layout": layout,
           "a_layouts": a_layout,
           "h_type": h_type,
           "bcc": bcc,
           "layout_3d": layout_3d,
           "hres": hres,
           "chain_id": chain_id, }

print("INFO: Domain, chain lists", cd_list)
print("INFO: Generating 2D layout for %s" % protein)

if channels_fn != None:
    channels = layout_loader_channels.Channels(channels_fn, "ALL", mol_online=args.mol_online)
else:
    channels = None

# Load protein
annotation = layout_loader_ligands.Annotation(annot_fn, "ALL", "ALL")

# at first read all domains and chains
chains = []
domains = {}
for i in cd_list:
    print("DEBUG: chain", i[0], " domain ", i[1], " processed")
    if i[0] not in chains:
        chains.append(i[0])
    if i[0] not in domains:
        domains[i[0]] = [i[1]]
    else:
        domains[i[0]].append(i[1])

# generate layout using domain, chains as a clusters and evolution algorithms
l = layout_loader_ligands.Layout(annotation=annotation)
l.compute_whole_protein(layouts, chains, domains, os.path.dirname(args.annotation_file))

if ligands_fn != None:
    ligands = layout_loader_ligands.Ligands(ligands_fn, "ALL", args.mol_online)
    l.compute_ligands(ligands, gcl=True, layout_3d=layout_3d)
else:
    ligands = None

if channels_fn != None:
    if args.mol_online:
        channels = layout_loader_channels.Channels(channels_fn, "ALL", mol_online=True)
    else:
        channels = layout_loader_channels.Channels(channels_fn, "ALL", mol_online=False)

    l.compute_channels(channels, os.path.dirname(annot_fn), layout_3d=layout_3d, flag="protein")
else:
    channels = None

# save the result layout
l.save_protein(chains, domains, h_type, bcc, hres, output, protein, layout_3d, ligands=ligands, channels=channels, mol_online=args.mol_online)
print("INFO: Protein layout saved into %s" % output)
