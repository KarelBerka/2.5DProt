#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import json
import sys
sys.path.append("..")
import layout_loader


# input is list of chain layouts

# Load params - cif file of whole protein domain and list of chains layouts

annot_fn = sys.argv[1]
output = sys.argv[2]
pp = sys.argv[3]
cd_list = []
layout = {}
a_layout = {}
h_type = {}
bcc = {}
layout_3d = {}
hres = {}

i = 4
# read protein lists - read layouts and a_layouts per domain, chain
while (i < len(sys.argv)):
    pcd_full = sys.argv[i]
    print("pcd_full", pcd_full)
    i = i + 1
    pcd = pcd_full[pcd_full.rfind("/")+1:]
    pcd_split = pcd.split("-")
    pcd_split = pcd_split[1].split(".")
    pcd_split = pcd_split[0].split("_")
    print("pcd_split", pcd_split)
    protein = pcd_split[0]
    domain = pcd_split[1][-2:]
    chain = pcd_split[1][:-2]
    if (len(chain) == 0):
        chain = domain
        domain = "ALL"
    if (len(domain) == 0):
        domain = "ALL"
    if protein == pp:
        cd_list.append((chain, domain))
    else:
        print("DEBUG: layout", pcd_full, "does not match to protein", pp)
    with open(pcd_full, "r") as fd:
        data = json.load(fd)
    if chain not in layout:
        layout[chain] = {}
        a_layout[chain] = {}
        h_type[chain] = {}
        bcc[chain] = {}
        layout_3d[chain] = {}
        hres[chain] = {}
    layout[chain][domain] = data['layout']
    a_layout[chain][domain] = data['angles']
    h_type[chain][domain] = data['h_type']
    bcc[chain][domain] = data['beta_connectivity']
    layout_3d[chain][domain] = data['3dlayout']
    hres[chain][domain] = data['helices_residues']

layouts = {"layout":layout, "a_layouts":a_layout, "h_type":h_type, "bcc":bcc, "layout_3d":layout_3d, "hres":hres}

print("INFO: Domain, chain lists", cd_list)
print("INFO: Generating 2D layout for %s" % pp)

# Load protein
annotation = layout_loader.Annotation(annot_fn, "ALL", "ALL")

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
l = layout_loader.Layout(annotation=annotation)
l.compute_whole_protein(layouts, chains, domains)

# save the result layout
l.save_protein(chains, domains, h_type, bcc, hres, output, pp, layout_3d)
print("INFO: Protein layout saved into %s" % output)
