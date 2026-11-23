#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import random
import logging
import utils
import math
import tempfile
import time
from collections import OrderedDict
import os
import subprocess
import json


def count_error(for_ligand, ligands, all_s, prime_s,
                       sses_sizes, sses_layout, sses_chain_id, l_s_distances, flags):

    if flags == "ligands":
        return count_ligand_error(for_ligand, ligands, all_s, prime_s,
                       sses_sizes, sses_layout, sses_chain_id, l_s_distances)
    else:
        return count_channel_error(for_ligand, ligands, all_s, prime_s,
                       sses_sizes, sses_layout, sses_chain_id, l_s_distances)


# Scale factor converting 3D Angstrom distances to 2D canvas coordinate units
ANGSTROM_TO_2D_SCALE = 0.38


def count_ligand_error(for_ligand, ligands, all_s, prime_s,
                       sses_sizes, sses_layout, sses_chain_id, l_s_distances):
    err = 0
    l_layout = ligands.layout
    l_chain_id = ligands.auth_chain

    for l1 in for_ligand:
        serr = 0

        for sse in all_s:
           l_dist = math.dist(l_layout[l1], sses_layout[sse])
           raw_td_dist = l_s_distances[l1][sse]
           td_dist = raw_td_dist * ANGSTROM_TO_2D_SCALE
           mindist = min(l_dist / ANGSTROM_TO_2D_SCALE, raw_td_dist)
           lt_diff = abs(l_dist - td_dist)
           if mindist < 15:
               dif = lt_diff * 5
           elif mindist < 25:
               dif = lt_diff * 3
           elif mindist < 35:
               dif = lt_diff / 2
           else:
               dif = 0
               # ligands are residues in the same ligand
           if sse not in prime_s and len(prime_s) >= 2:
               # secondary sses and there is nontrivial number of prime sses
               dif = dif / 2
           if (sses_chain_id[sse] != l_chain_id[l1]) and (l_chain_id[l1] != ""):
               # the sse is in the different chain
               dif = dif / 4
           serr += dif

        lerr = 0
        for l2 in ligands.names:
           l_dist = math.dist(l_layout[l1], l_layout[l2])
           raw_td_dist = l_s_distances[l1][l2]
           td_dist = raw_td_dist * ANGSTROM_TO_2D_SCALE
           mindist = abs(min(l_dist / ANGSTROM_TO_2D_SCALE, raw_td_dist))
           lt_diff = abs(l_dist - td_dist)
           if mindist < 10:
               dif = lt_diff * 3
           elif mindist < 25:
               dif = lt_diff
           else:
               dif = 0
           lerr = lerr + dif
        err += serr + lerr
    return err


def count_channel_error(for_ligand, ligands, all_s, prime_s,
                       sses_sizes, sses_layout, sses_chain_id, l_s_distances):
    err = 0
    l_layout = ligands.layout
    l_chain_id = ligands.auth_chain

    for l1 in for_ligand:
        serr = 0

        for sse in all_s:
           l_dist = math.dist(l_layout[l1], sses_layout[sse])
           raw_td_dist = l_s_distances[l1][sse]
           td_dist = raw_td_dist * ANGSTROM_TO_2D_SCALE
           mindist = min(l_dist / ANGSTROM_TO_2D_SCALE, raw_td_dist)
           lt_diff = (abs(l_dist - td_dist) + 3*max(0, l_dist-td_dist))/4
           if mindist < 15:
               dif = lt_diff * (9-mindist/4)* (9-mindist/4)
           elif mindist < 25:
               dif = lt_diff * (9-mindist/4)*(9-mindist/4)
           elif mindist < 35:
               dif = lt_diff / 2
           else:
               dif = 0
               # ligands are residues in the same ligand
           if sse not in prime_s and len(prime_s) >= 2:
               # secondary sses and there is nontrivial number of prime sses
               dif = dif / 2
           if (sses_chain_id[sse] != l_chain_id[l1]) and (l_chain_id[l1] != ""):
               # the sse is in the different chain
               dif = dif / 4
           serr += dif

        lerr = 0
        for l2 in ligands.names:
           l_dist = math.dist(l_layout[l1], l_layout[l2])
           raw_td_dist = l_s_distances[l1][l2]
           td_dist = raw_td_dist * ANGSTROM_TO_2D_SCALE
           mindist = abs(min(l_dist / ANGSTROM_TO_2D_SCALE, raw_td_dist))
           lt_diff = max(0, l_dist - td_dist)
           if mindist < 10:
               dif = lt_diff * 3
           elif mindist < 25:
               dif = lt_diff
           else:
               dif = 0
           lerr = lerr + dif
        err += serr + lerr
    return err


def draw_curve(bef_start_layout, start_layout, end_layout, af_end_layout, buffer, ligands):
    # between start_layout point and end_layout point distribute buffer points
    # (there are from 1 to 7 points)
    # points are on curves, based on point before start (bef_start_layout) and
    # point after end af_end_layout
    counter = 1
    for j in buffer:

        if len(buffer) == 7:
            # protein 4v4g (22000 ligands)
            if counter == 1:
               r0 = 0.70
               r1 = 0.30
               r2 = 0
            elif counter == 2:
               r0 = 0.85
               r1 = 0.15
               r2 = 0
            elif counter == 3:
               r0 = 0.9
               r1 = 0.1
               r2 = 0
            elif counter == 4:
               r0 = 1
               r1 = 0
               r2 = 0
            elif counter == 5:
               r0 = 0.9
               r1 = 0
               r2 = 0.1
            elif counter == 6:
               r0 = 0.85
               r1 = 0
               r2 = 0.15
            elif counter == 7:
               r0 = 0.70
               r1 = 0
               r2 = 0.30

        if len(buffer) == 6:
            # protein ?4v9p,4v9o (18000 ligands), ?4v4g (22000 ligands)
            if counter == 1:
               r0 = 0.70
               r1 = 0.30
               r2 = 0
            elif counter == 2:
               r0 = 0.85
               r1 = 0.15
               r2 = 0
            elif counter == 3:
               r0 = 0.9
               r1 = 0.1
               r2 = 0
            elif counter == 4:
               r0 = 0.9
               r1 = 0
               r2 = 0.1
            elif counter == 5:
               r0 = 0.85
               r1 = 0
               r2 = 0.15
            elif counter == 6:
               r0 = 0.70
               r1 = 0
               r2 = 0.30

        if len(buffer) == 5:
            # protein 4v9p, 4v9o(18000 ligands)
            if counter == 1:
               r0 = 0.70
               r1 = 0.30
               r2 = 0
            elif counter == 2:
               r0 = 0.85
               r1 = 0.15
               r2 = 0
            elif counter == 3:
               r0 = 1
               r1 = 0
               r2 = 0
            elif counter == 4:
               r0 = 0.85
               r1 = 0
               r2 = 0.15
            elif counter == 5:
               r0 = 0.70
               r1 = 0
               r2 = 0.30

        if len(buffer) == 4:
            if counter == 1:
               r0 = 0.70
               r1 = 0.30
               r2 = 0
            elif counter == 2:
               r0 = 0.85
               r1 = 0.15
               r2 = 0
            elif counter == 3:
               r0 = 0.85
               r1 = 0
               r2 = 0.15
            elif counter == 4:
               r0 = 0.70
               r1 = 0
               r2 = 0.30

        if len(buffer) == 3:
            if counter == 1:
               r0 = 0.80
               r1 = 0.20
               r2 = 0
            elif counter == 2:
               r0 = 1
               r1 = 0
               r2 = 0
            elif counter == 3:
               r0 = 0.8
               r1 = 0
               r2 = 0.2

        if len(buffer) == 2:
            if counter == 1:
               r0 = 0.80
               r1 = 0.20
               r2 = 0
            elif counter == 2:
               r0 = 0.8
               r1 = 0
               r2 = 0.2

        if len(buffer) == 1:
            if counter == 1:
               r0 = 0.70
               r1 = 0.15
               r2 = 0.15

        if ((r1 != 0) and bef_start_layout) and ((r2 != 0) and af_end_layout):
            ligands.layout[j][0] = start_layout[0] + float(counter) / (len(buffer) + 1) * (end_layout[0] - start_layout[0]) * r0 + \
                float(counter) / (len(buffer) + 1) * (start_layout[0] - bef_start_layout[0]) * r1 + \
                float(counter) / (len(buffer) + 1) * (af_end_layout[0] - end_layout[0]) * r2
            ligands.layout[j][1] = start_layout[1] + float(counter) / (len(buffer) + 1) * (end_layout[1] - start_layout[1]) * r0 + \
                float(counter) / (len(buffer) + 1) * (start_layout[1] - bef_start_layout[1]) * r1 + \
                float(counter) / (len(buffer) + 1) * (af_end_layout[1] - end_layout[1]) * r2
        elif ((r1 != 0) and bef_start_layout):
            ligands.layout[j][0] = start_layout[0] + float(counter) / (len(buffer) + 1) * (end_layout[0] - start_layout[0]) * r0 + \
                float(counter) / (len(buffer) + 1) * (start_layout[0] - bef_start_layout[0]) * r1
            ligands.layout[j][1] = start_layout[1] + float(counter) / (len(buffer) + 1) * (end_layout[1] - start_layout[1]) * r0 + \
                float(counter) / (len(buffer) + 1) * (start_layout[1] - bef_start_layout[1]) * r1
        elif ((r2 != 0) and af_end_layout):
            ligands.layout[j][0] = start_layout[0] + float(counter) / (len(buffer) + 1) * (end_layout[0] - start_layout[0]) * r0 + \
                float(counter) / (len(buffer) + 1) * (af_end_layout[0] - end_layout[0]) * r2
            ligands.layout[j][1] = start_layout[1] + float(counter) / (len(buffer) + 1) * (end_layout[1] - start_layout[1]) * r0 + \
                float(counter) / (len(buffer) + 1) * (af_end_layout[1] - end_layout[1]) * r2
        else:
            ligands.layout[j][0] = start_layout[0] + float(counter) / (len(buffer) + 1) * (end_layout[0] - start_layout[0]) * r0
            ligands.layout[j][1] = start_layout[1] + float(counter) / (len(buffer) + 1) * (end_layout[1] - start_layout[1]) * r0

        counter += 1


def tmp_write(tmp, ligands_names, l_s_distances):
    tmp.write("lsdistances2 ")
    for sse1 in ligands_names:
        tmp.write(f"{sse1}:")
        for sse2 in ligands_names:
            tmp.write(f"{l_s_distances[sse1][sse2]},")
        tmp.write(";")
    tmp.write("\n")



def modify_ligand_layout_via_evolution(ligands, all_s, prime_s, sses_sizes,
                                       sses_layout, sses_chain_id,
                                       l_s_distances, layout_3d, flags={},
                                       flag="ligands", tmp_dir="",
                                       channels_positions={}, sses_positions={}):

    if len(ligands.names) > 4000:
       # maximum 4v4g (22000 lig), 4v9o, 4v9p (18000 lig)
       max_ligands = 3000
    elif len(ligands.names) > 2000:
       max_ligands = 2000
    else:
       max_ligands = len(ligands.names)
    max_generation = 100000

    ligands_names_orig = ligands.names
    if len(ligands.names) > max_ligands:
        # if there is too much ligands e.g. 22000 (4v4g)
        # this process is too slow and the diagram too messy,
        # thus we choose only limited number of ligands

       score = max_ligands / len(ligands.names)
       names = []
       for c in range(len(ligands.names)):
           if (c+1) * score > len(names):
               names.append((ligands.names)[c])
       ligands.names = names
    l_number = len(ligands.names)

    # set the maximal step size value
    if "step_size" in flags:
        STEP_SIZE = flags["step_size"]
    else:
        STEP_SIZE = 10
        for i in ligands.names:
            for j in ligands.names:
                dist = math.dist(ligands.center[i], ligands.center[j])
                if dist > STEP_SIZE:
                     STEP_SIZE = dist

            if layout_3d:
                for c in layout_3d:
                    for d in layout_3d[c]:
                        for j in layout_3d[c][d]:
                            line_diff = [layout_3d[c][d][j][1][0]- layout_3d[c][d][j][0][0],
                                 layout_3d[c][d][j][1][1]- layout_3d[c][d][j][0][1],
                                 layout_3d[c][d][j][1][2]- layout_3d[c][d][j][0][2]]
                            dist = utils.point_to_line_dist(layout_3d[c][d][j][0], line_diff, ligands.center[i],)
                            if dist > STEP_SIZE:
                                STEP_SIZE = dist


    # set the number of iterations
    max_generations = 300000
    if l_number < 25:
        CHECK_EVERY = l_number * 100   # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
    elif l_number < 80:
        CHECK_EVERY = l_number * 75
    elif l_number < 2000:
        CHECK_EVERY = l_number * 50 
    elif l_number == 2001:
        CHECK_EVERY = 60000
    else: 
        CHECK_EVERY = 100000
        max_generations = 500000

    if flag == "channels" and os.path.isfile("../gcc/ligands_main_loop"): 
        # do main loop in language C -
        # prepare language c input:
        aux_sses_sizes = {}
        for i in sses_sizes:
            aux_sses_sizes[i] = int(sses_sizes[i])

        tmpname_in = f"{tmp_dir}/c_input_channels"
        all_s = list(all_s)
        prime_s = list(prime_s)
        data = {
             "ligandsnames": ligands.names,                 #
             "ligands_layout": ligands.layout,              #
             "all_s": all_s,                                #
             "prime_s": prime_s,                            #
             "sses_sizes": aux_sses_sizes,                  #
             "sses_layout": sses_layout,                    #
             "lsdistances": l_s_distances,                  #
             "flag": flag,                                  # -
             "STEP_SIZE": int(STEP_SIZE),                   #
             "CHECK_EVERY": CHECK_EVERY}                    #

        with open(tmpname_in, "w") as tmp:
             json.dump(data, tmp)

        # all input prepared
        # do do main loops
        
        command = f"../gcc/ligands_main_loop {tmpname_in} {tmp_dir}/c_output_channels"
        print("command", command)
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()

        # read output results
        layout = {}
        with open(f"{tmp_dir}/c_output_channels", "r") as fd:
            line = fd.readline()
            while line:
                sline=line.split(",")
                if sline[0] == "layout":
                    ligands.layout[sline[1]] = [float(sline[2]),float(sline[3])]
                line = fd.readline()
        return ligands.layout

    else:
      logging.debug("We modify ligands positions")
      logging.debug(f" ligands number {l_number} CHECK_EVERY {CHECK_EVERY} ligands {ligands.names}")
      MINIMAL_IMPROVEMENT = 0.02   # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout

      generation = 1
      error_checkpoint = None
      coef = 1
##      d_chan = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


      while True:
        generation += 1
        # choose the moved sse
        for_ligand = random.choice(ligands.names)
        # count original error
        error_old = count_error([for_ligand], ligands, all_s, prime_s,
                                       sses_sizes, sses_layout, sses_chain_id,
                                       l_s_distances, flag)
        # set the new coordinates
        dx = random.random() * (STEP_SIZE/coef) - STEP_SIZE/(2*coef)
        dy = random.random() * (STEP_SIZE/coef) - STEP_SIZE/(2*coef)
        ligands.layout[for_ligand] = [ligands.layout[for_ligand][0] + dx, ligands.layout[for_ligand][1] + dy]

        # count new error
        error_new = count_error([for_ligand], ligands, all_s, prime_s,
                                       sses_sizes, sses_layout, sses_chain_id,
                                       l_s_distances, flag)
        if error_new < error_old:
            # there was nontrivial change -> leave the layout as it is
            percentage = 100 - error_new / error_old * 100
            logging.debug("[%s/%s] By move of %s(x:%.3f,y:%.3f - z %s)(x,y:%.3f,%.3f) error changed from %.3f to %.3f, i.e. %.2f%% "
                          % (generation, max_generations, for_ligand, ligands.layout[for_ligand][0], ligands.layout[for_ligand][1],
                             STEP_SIZE/(2*coef), ligands.layout[for_ligand][0]-dx, ligands.layout[for_ligand][1]-dy, error_old, error_new, percentage))
            error_old = error_new
        else:
            # we have to return layout back
            ligands.layout[for_ligand] = [ligands.layout[for_ligand][0]-dx, ligands.layout[for_ligand][1]-dy]

        if generation % CHECK_EVERY == 0:
            coef = coef * 2
            # If there was no significant improvement in last CHECK_EVERY iterations, there is no reason to continue
            errors_summary = count_error(ligands.names, ligands, all_s,
                                                prime_s, sses_sizes,
                                                sses_layout, sses_chain_id,
                                                l_s_distances, flag)
            if generation == CHECK_EVERY:
                # If we are in the first round, there might be "lost" ligands which are in absolutely different position
                # if there is a lot of ligands, we have to put them to the best positions
                e = {}
                esum = 0
                lxmin = None
                lxmax = None
                lymin = None
                lymax = None
                for l_name in ligands.names:
                    # we found the average errors of ligands
                    e[l_name] = count_error([l_name], ligands, all_s, prime_s,
                                                   sses_sizes, sses_layout,
                                                   sses_chain_id, l_s_distances, flag)
                    esum += e[l_name]
                    lx = ligands.layout[l_name][0]
                    ly = ligands.layout[l_name][1]
                    if lxmin is None or lx < lxmin:
                        lxmin = lx
                    if lxmax is None or lx > lxmax:
                        lxmax = lx
                    if lymin is None or ly < lymin:
                        lymin = ly
                    if lymax is None or ly > lymax:
                        lymax = ly
                eavg = esum / len(e)
                auxlist = []
                for lig in e:
                    if e[lig] > eavg * 1.8:
                         auxlist.append(lig)
                for lig in auxlist:
                    # if the error is bigger then twice the standart - we the ligand could have wrong position
                    buffer = ligands.layout[lig]
                    for i in range(100):
                        # we try to find better position for the ligand (we go through the net 5 x 5 points through the whole
                        # square where are ligands now
                        ligands.layout[lig] = (lxmin+(lxmax-lxmin)*(i % 10)/9, lymin+(lymax-lymin)*(i // 10)/9)
                        error_new = count_error([lig], ligands, all_s, prime_s,
                                                       sses_sizes, sses_layout,
                                                       sses_chain_id, l_s_distances, flag)
                        if error_new >= e[lig]:
                            ligands.layout[lig] = buffer
                        else:
                            e[lig] = error_new
                            buffer = ligands.layout[lig]

            if error_checkpoint is None:
                if l_number < 100:
                    coef = coef * 2
            # if we have small improvenment or ligands are too distant ad error is zero stop
            if (error_checkpoint is not None) and ((errors_summary < 0.1) or (error_checkpoint - errors_summary) < (errors_summary * MINIMAL_IMPROVEMENT)):
                logging.debug("aaaa In last %s generations (step %s), total error decreased from %.2f to %.2f, which is below minimal required %.2f%%, so returning what we have"
                              % (CHECK_EVERY, generation, error_checkpoint, errors_summary, MINIMAL_IMPROVEMENT * 100))
                break

            error_checkpoint = errors_summary

        if generation >= max_generations:
            logging.debug("aaaa Finished after %s generations, maybe with more generations we can get better results? {max_generations}", max_generations)
            break

    if ligands_names_orig != ligands.names:
        # we have to return all missing ligands to their best positions
        ligands_set = ligands.names
        ligands.names = ligands_names_orig
        buffer = []
        for i in ligands.names:
            # we go through intervals, where ligands have their layouts (from start to end)
            # the ligands in the middle (buffer) is set based on them
            if i in ligands_set:
                c = ligands_set.index(i)
                if c > 0:
                    start_layout = ligands.layout[ligands_set[c-1]]
                    end_layout = ligands.layout[i]
                    if c > 1:
                        bef_start_layout = ligands.layout[ligands_set[c-2]]
                    else:
                        bef_start_layout = None
                    if len(ligands_set) > c+1:
                        af_end_layout = ligands.layout[ligands_set[c+1]]
                    else:
                        af_end_layout = None
                    draw_curve(bef_start_layout, start_layout, end_layout, af_end_layout, buffer, ligands)

                    buffer = []
            else:
                buffer.append(i)

    return ligands.layout
