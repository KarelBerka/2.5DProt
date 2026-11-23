#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Find best universal start from "number" random - choosen protein-domain-chains
#
# Usage:
# ./find_best_universal_start.py 5 generated-1.10.630.10/1.10.630.10/uc/*annotated*

import utils_rmsds
import random
import logging


def find_universal_start(pcd_list, cif_directory, number):
    # pymol = utils_rmsds.pymol_init()

    # choose "number" pdcs
    if len(pcd_list) <= number:
        sample_list = pcd_list
    else:
        sample_list = random.sample(pcd_list, number)

    # find its informations
#    data = {}
    logging.getLogger("utils_rmsds").setLevel(logging.CRITICAL)
    for pcd in sample_list:
        protein, chain, domain, dn = pcd
#        protein_filename = f"{cif_directory}/{dn}.cif"
#
#        protein = loader.Protein(protein_filename)
#        query = "%s:%s%s" % (protein, chain, domain)
#
#        try:
#            pymol_var_chained = utils_rmsds.pymol_load(pymol, protein, chain, domain)
#        except:
#            logging.error(f"utils_rmsd._pymol_load {protein} {chain} {domain} failed")
#        else:
#            data[query] = {'layout_filename': pcd,
#                           'protein_filename': protein_filename,
#                           'protein': protein,
#                           'chain': chain,
#                           'domain': domain,
#                           'pymol_var_chained': pymol_var_chained,
#            }
#
#
#    # for each pdc pair count rmsd
#    rmsd_3d_pymol = {}
#    for p1,d1 in data.items():
#        if p1 not in rmsd_3d_pymol:
#            rmsd_3d_pymol[p1] = {}
#        for p2,d2 in data.items():
#            if p2 in rmsd_3d_pymol and p1 in rmsd_3d_pymol[p2]:
#                rmsd_3d_pymol[p1][p2] = rmsd_3d_pymol[p2][p1]
#                continue
#            if p1 == p2:
#                rmsd_3d_pymol[p1][p2] = 0.0
#                continue
#            else:
#                try:
#                    rmsd_3d_pymol[p1][p2], rmsd_3d_pymol_match = \
#                        utils_rmsds.get_rmsd_3d_pymol(pymol, d1['pymol_var_chained'], d2['pymol_var_chained'])
#                except ValueError:
#                    rmsd_3d_pymol[p1][p2] = 99
#
#    # for each pdc count the sum of rmsds with the others
#    rmsd_sum = {}
#    for i in rmsd_3d_pymol:
#       rmsd_sum[i] = 0
#       for j in rmsd_3d_pymol:
#           rmsd_sum[i] = rmsd_sum[i] + rmsd_3d_pymol[i][j]
#
#    # find the minimum - its the winner
#    rmsd_min = None
#    pdc = None
#    for i in rmsd_sum:
#        if (rmsd_min == None ) or (rmsd_min > rmsd_sum[i]):
#           pdc = i
#           rmsd_min = rmsd_sum[i]
#    # we have the winner
#    #TODO: have to be fixed - should be pdc
    return pcd
