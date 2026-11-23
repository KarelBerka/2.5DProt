#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import math
import random
import logging
from collections import OrderedDict
import utils
import utils_angleslayout
import tempfile
import subprocess
import json
import sys

def get_fake_layout(clusters, start_layout, relevant_sses, layout):
    """
    Just returns layout where SSEs are placed with coords [0, 0], [1, 1]...
    """
    for c in clusters:
        sses = c.split(":")

        is_in = False
        for i in sses:
           if i in relevant_sses:
              is_in = True
        if is_in == False:
            # no relevant sses, not interested
            continue

        # this have to be count
        ccoords = (0, 0)
        counter = 0
        for j in sses:
            if j in start_layout:
               ccoords = (ccoords[0]+start_layout[j][0],
                           ccoords[1]+start_layout[j][1])
               counter += 1
        if counter > 0:
            layout[c] = (ccoords[0]/counter, ccoords[1]/counter)
        else:
            layout[c] = (0, 0)


    return layout

def find_best_rotation_reflection(movable_clusters, layout, a_layout, \
                c_layout, ca_layout, sses_sizes, relevant_sses, movable_sses, \
                sses_distance, start_layout, flags):
    for for_helic in movable_clusters:
        sses = for_helic.split(":")

        if (len(sses) > 1):
            # cluster is nontrivial
            minerr = None
            angle = None
            reflection = False

            # at first try angles 0, 10, 20 ... 350
            i = 0
            while i < 36:
                a_layout[for_helic] = i * 10.0
                e1, e2 = count_error(layout, a_layout, c_layout, ca_layout, sses_sizes, relevant_sses, \
                    movable_sses, for_helic, sses_distance, start_layout, flags)
                err = e1 + e2
                if (minerr == None) or (minerr > err):
                    angle = i * 10.0
                    minerr = err
                i = i + 1

            # now we have to test reflection and try angles 0, 10, 20, ... 350
            i = 0
            # reflection itself
            for sse in c_layout[for_helic]:
                c_layout[for_helic][sse] = (c_layout[for_helic][sse][0], -c_layout[for_helic][sse][1])
                ca_layout[for_helic][sse] = 360 - ca_layout[for_helic][sse]

            # and the rotation afterwards
            while i < 36:
                a_layout[for_helic] = i * 10.0
                e1,e2 = count_error(layout, a_layout, c_layout, ca_layout, sses_sizes, relevant_sses, \
                    movable_sses, for_helic, sses_distance, start_layout, flags)
                err = e1 + e2
                if (minerr > err):
                    angle = i * 10.0
                    reflection = True
                    minerr = err
                i = i + 1

            # choose the best rotation / reflection and use it
            if (reflection == False):
                # we have to reflect cluster back
                for sse in c_layout[for_helic]:
                    c_layout[for_helic][sse] = (c_layout[for_helic][sse][0], -c_layout[for_helic][sse][1])
                    ca_layout[for_helic][sse] = 360 - ca_layout[for_helic][sse]
            a_layout[for_helic] = angle
    return a_layout, c_layout, ca_layout


def get_distance_matrix_2d(layout, only_for, selected_cluster):
    """
    Compute distance matrix of given layout
      layout ... {'A': (x_a, y_a), 'B': (x_b, y_b)...}
      only_for ... if set to some helic letter, only compute distances for
                   that helic,
                   if set to list only compute distances for given list of helices,
                   if set to None, compute for each
    """
    out = {}
    for key1 in only_for:
        out[key1] = {}
        for k2 in selected_cluster:
            for key2 in k2.split(":"):
                out[key1][key2] = math.dist(layout[key1], layout[key2])
    return out

def plain_error_counting(dist1, dist2, size1, size2, const):
   if min(dist1,dist2) < const:
      res = abs(dist1 - dist2) * (size1+size2)
   else:
      res = 0
   return res


def get_measured_error_2d(a, b, sizes, relevant_sses, movable_sses, selected_sses, flags):
    """
    Compute sum of absolute values of differences in two provided dictionaries
      a ... {'C': {'A': 5, 'B': 13}, 'A': {...}, 'B': {...}}
      b ... {'C': {'A': 4, 'B': 12}, 'A': {...}, 'B': {...}}
      keys ... ['A', 'B', 'C'], defaults to keys from a
    """
    summary = 0
    for key1 in relevant_sses:
        if key1 not in a:
            continue
        for k2 in selected_sses:
            for key2 in k2.split(":"):
                if key1 != key2 and key2 in a[key1]:
                    dif = plain_error_counting(a[key1][key2], b[key1][key2], sizes[key1], sizes[key2], 40)
                    if ("dbs" in flags) and flags["dbs"] == "spike-protein":
                        # for spike protein dbs, chains 1,2,3 layout is generated together. There is necessary to have more wieght in
                        # keeping the same distances in each protein chain. The same situation is in spike proteins chains 3-6 (or more if they are present)
                        if (key1[1] == "1" and key2[1] == "1") or (key1[1] == "2" and key2[1] == "2") or (key1[1] == "3" and key2[1] == "3") or \
                           (key1[1] == "4" and key2[1] == "4") or (key1[1] == "5" and key2[1] == "5") or (key1[1] == "6" and key2[1] == "6") or \
                           (key1[1] == "6" and key2[1] == "6"):
                            summary += dif * 4
                        elif (key1[1] == "4" and key2[1] == "1") or (key1[1] == "4" and key2[1] == "2") or (key1[1] == "4" and key2[1] == "3") or \
                           (key1[1] == "5" and key2[1] == "1") or (key1[1] == "5" and key2[1] == "2") or (key1[1] == "5" and key2[1] == "3") or \
                           (key1[1] == "6" and key2[1] == "1") or (key1[1] == "6" and key2[1] == "2") or (key1[1] == "6" and key2[1] == "3") or \
                           (key2[1] == "4" and key1[1] == "1") or (key2[1] == "4" and key1[1] == "2") or (key2[1] == "4" and key1[1] == "3") or \
                           (key2[1] == "5" and key1[1] == "1") or (key2[1] == "5" and key1[1] == "2") or (key2[1] == "5" and key1[1] == "3") or \
                           (key2[1] == "6" and key1[1] == "1") or (key2[1] == "6" and key1[1] == "2") or (key2[1] == "6" and key1[1] == "3")or \
                           (key2[1] == "7" and key1[1] == "1") or (key2[1] == "7" and key1[1] == "2") or (key2[1] == "7" and key1[1] == "3"):
                            summary += dif / 6
                        else:
                            summary += dif
                    else:
                        # for other dbs layouts are computed per chains
                        summary += dif
    return summary


def move_error(start_layout, layout, sses_sizes, relevant_sses, selected_sses, flags):
    m_error = 0

    for s in selected_sses:
        sses_in = 0
        sses_out = 0
        md = 0
        for i in s.split(":"):
            if (i in layout) and (i in sses_sizes):
                if (i in start_layout):
                    dif = math.dist(start_layout[i],layout[i])
                    if ("dbs" in flags) and flags["dbs"] == "spike-protein":
                        # for spike protein dbs, chains 1,2,3 are computed together, their cumulativ shape is
                        # cone, which means there are several minimal positions. We have to have bigger stress not
                        # to reach different minimals for single bunch of layouts, thus we have to increase 
                        # start layout weight
                        if (i[1] == "1") or (i[1] == "2") or (i[1] == "3"):
                            md += dif * len(relevant_sses) * sses_sizes[i] / 20
                        else:
                            md += dif * len(relevant_sses) * sses_sizes[i] / 20
                    else:
                        # for other dbs there are not so big and problematic structures
                        m_error += dif * len(relevant_sses) * sses_sizes[i] / 20
                    sses_in += sses_sizes[i]
                else:
                    sses_out += sses_sizes[i]
        if sses_out != 0 and sses_in != 0:
            # if only part of sheet is in start layout, we have to multiply
            # move error by the size of whole strands in sheet and divide by
            # the size of strands which are in start layout
            m_error += md * (sses_in+sses_out) / sses_in
        else:
            m_error += md

    return m_error



def count_error(layout, a_layout, c_layout, ca_layout, sses_sizes, relevant_sses, \
    movable_sses, for_helic, sses_distance, start_layout, flags):

    sses_layout, sses_layout_angles = count_sse_coords_based_on_cluster_coords(\
         layout, a_layout, c_layout, ca_layout)

    if (type(for_helic) != list):
       fh = {for_helic}
    else:
       fh = for_helic

    # we have to compute only the positions
    sd_layout = get_distance_matrix_2d(sses_layout, relevant_sses, fh)
    error = get_measured_error_2d(sses_distance, sd_layout, \
        sses_sizes, relevant_sses, movable_sses, fh, flags)
    m_error = move_error(start_layout, sses_layout, sses_sizes, \
        relevant_sses, fh, flags)
    return error,m_error


def count_sse_coords_based_on_cluster_coords(layout, a_layout, c_layout, ca_layout):

    sses_layout = OrderedDict()
    sses_angles_layout = OrderedDict()

    for cl in layout:
        angle = a_layout[cl]/180 * math.pi
        for sse in cl.split(":"):
            # count down sse1 coords based on clusters coords, cluster layout and cluster_layout_angles
            x = c_layout[cl][sse][0]
            y = c_layout[cl][sse][1]
            ssedif = math.sqrt(x*x+y*y)
            if y == 0:
               beta = 1
            else:
               beta = math.atan(x/y)
            if y<0:
                coord_sse1_x = layout[cl][0] - (ssedif*math.sin(angle+beta))
                coord_sse1_y = layout[cl][1] - (ssedif*math.cos(angle+beta))
            else:
                coord_sse1_x = layout[cl][0] + (ssedif*math.sin(angle+beta))
                coord_sse1_y = layout[cl][1] + (ssedif*math.cos(angle+beta))

            sses_layout[sse] = (coord_sse1_x, coord_sse1_y)
            sses_angles_layout[sse] = (a_layout[cl] + ca_layout[cl][sse]) % 360
    return sses_layout, sses_angles_layout


def modify_layout_via_evolution(layout, a_layout, c_layout, ca_layout, \
        sses_distance, sses_sizes, movable_clusters, fixed_clusters, start_layout, \
        tmp_dir, flags = {}):

    dmax = 0
    # TODO: should not be all sse_distance members only movable_clusters + fixed clusters together
    for i in sses_distance:
        for j in sses_distance[i]:
            if sses_distance[i][j] > dmax:
                dmax = sses_distance[i][j]

    if "step_size" in flags:
        STEP_SIZE = flags["step_size"]
    else:
        STEP_SIZE = dmax/2

    if len(movable_clusters) < 6:
        CHECK_EVERY = 1000
    elif len(movable_clusters) < 12:
        CHECK_EVERY = 2000
    elif len(movable_clusters) < 25:
        CHECK_EVERY = 3000   # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
    elif len(movable_clusters) < 80:
        CHECK_EVERY = 4000   # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
    else:
        CHECK_EVERY = 5000

    logging.debug("We modify layout and nontrivial clusters orientation")
    logging.debug(f"movabe_clusters number {len(movable_clusters)} CHECK_EVERY {CHECK_EVERY} movable clusters {movable_clusters}")

    movable_sses=[]
    for cl in movable_clusters:
       for i in cl.split(":"):
          movable_sses.append(i)

    relevant_sses = movable_sses
    for cl in fixed_clusters:
       for i in cl.split(":"):
          relevant_sses.append(i)

    # try to find the best rotation/reflection of all nontrivial clusters
    a_layout, c_layout, ca_layout = find_best_rotation_reflection(\
                movable_clusters, layout, a_layout, c_layout, ca_layout, \
                sses_sizes, relevant_sses, movable_sses,  sses_distance, \
                start_layout, flags)

    if True:
        # use language c to be faster
        if STEP_SIZE >0.01:
            tmp_input_file = f"{tmp_dir}/c_input_sses"
            tmp_output_file = f"{tmp_dir}/c_output_sses"
            aux = OrderedDict(sses_sizes)
            data = {
               "movable_clusters": list(movable_clusters), #
               "layout": layout,                           #
               "a_layout": a_layout,                       #
               "c_layout": c_layout,
               "ca_layout": ca_layout,
               "sses_sizes": aux,                          #
               "relevant_sses": relevant_sses,             #
               "movable_sses": movable_sses,               #
               "sses_distance": sses_distance,             #
               "start_layout": start_layout,               #
               "flags": flags,                             #  do not red now
               "STEP_SIZE": STEP_SIZE,                     #
               "CHECK_EVERY": CHECK_EVERY}                 #
            with open(tmp_input_file, "w") as tmp:
               json.dump(data, tmp)

            command = f"../gcc/sses_main_loop {tmp_input_file} {tmp_output_file}"
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            r = p.wait()

            layout = {}
            with open(f"{tmp_output_file}", "r") as fd:
                  line = fd.readline();
                  while line:
                      sline=line.split(",")
                      if sline[0] == "full_layout":
                          layout[sline[1]] = [float(sline[2]),float(sline[3])]
                      if sline[0] == "full_cluster_layout":
                           for cl in c_layout:
                               for sse in c_layout[cl]:
                                   if sline[2] == sse:
                                       c_layout[cl][sse] = [float(sline[3]),float(sline[4])]
                                       ca_layout[cl][sse] = float(sline[5])
                      line = fd.readline()
        return layout, a_layout, c_layout, ca_layout
    else:
      # use python language - slower version
      ANGLES_SIZE = 45
      MINIMAL_IMPROVEMENT = 0.01   # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout
      max_generations = 60000
      generation = 0
      error_checkpoint = None
      coef = 1

      while True:
        generation += 1
        # choose the moved sse
        for_helic = random.choice(list(movable_clusters))
        esn, esnm = count_error(layout, a_layout, c_layout, ca_layout, \
             sses_sizes, relevant_sses, movable_sses, for_helic, sses_distance, start_layout, flags)
        errors_summary_old = esn+esnm

        # set the new coordinates
        dx = random.random() * (STEP_SIZE/coef) - STEP_SIZE/(2*coef)
        dy = random.random() * (STEP_SIZE/coef) - STEP_SIZE/(2*coef)
        layout[for_helic] = (layout[for_helic][0]+dx, layout[for_helic][1]+dy)

        # count new error
        esn, esnm = count_error(layout, a_layout, c_layout, ca_layout, \
             sses_sizes, relevant_sses, movable_sses, for_helic, sses_distance, start_layout, flags)
        errors_summary_new = esn+esnm
        if errors_summary_new < errors_summary_old:
            percentage = 100 - errors_summary_new / errors_summary_old * 100
            logging.debug("[%s/%s] By move of %s(x:%.3f,y:%.3f - z %s)(x,y:%.3f,%.3f) error changed from %.3f to %.3f %3f %3f, i.e. %.2f%% "
                % (generation, max_generations, for_helic, layout[for_helic][0], layout[for_helic][1], STEP_SIZE/(2*coef), layout[for_helic][0]-dx, layout[for_helic][1]-dy, 
                    errors_summary_old, errors_summary_new, esn, esnm, percentage))
            errors_summary_old = errors_summary_new
        else:
            layout[for_helic] = (layout[for_helic][0]-dx, layout[for_helic][1]-dy)

        # now we try to rotate nontrivial clusters to improve their orientation
        if (len(for_helic.split(":")) > 1):
            counter = 1
            while counter < 3:
               counter = counter + 1
               difangle = random.random() * (ANGLES_SIZE/coef) - ANGLES_SIZE/(2*coef)
               a_layout[for_helic] = (a_layout[for_helic] + difangle ) % 360

               esn, esnm = count_error(layout, a_layout, c_layout, ca_layout, \
                  sses_sizes, relevant_sses, movable_sses, for_helic, sses_distance, start_layout, flags)
               errors_summary_new = esn+esnm

               if errors_summary_new < errors_summary_old:
                   percentage = 100 - errors_summary_new / errors_summary_old * 100
                   logging.debug("[%s/%s] By rotation of %s ( %.3f from %.3f  dif %.3f (max %s) ) error changed from %.3f to %.3f, i.e. %.2f%%"
                       % (generation, max_generations, for_helic, a_layout[for_helic], (a_layout[for_helic] - difangle ) % 360, difangle, ANGLES_SIZE/(2*coef),
                          errors_summary_old, errors_summary_new, percentage))
                   errors_summary_old = errors_summary_new
               else:
                   a_layout[for_helic] = (a_layout[for_helic] - difangle ) % 360

        # we have to try reflection (if the cluster is nontrivial and <500 generations
        if (len(for_helic.split(":")) > 1) and (generation < 500):

            for i in c_layout[for_helic]:
                c_layout[for_helic][i] = (c_layout[for_helic][i][0], -c_layout[for_helic][i][1])
                ca_layout[for_helic][i] = 360 - ca_layout[for_helic][i]
            difangle = random.random() * 360 - 180
            a_layout[for_helic] = a_layout[for_helic]+difangle

            esn, esnm = count_error(layout, a_layout, c_layout, ca_layout, \
                 sses_sizes, relevant_sses, movable_sses, for_helic, sses_distance, start_layout, flags)
            errors_summary_new = esn+esnm


            if errors_summary_new < errors_summary_old:
                percentage = 100 - errors_summary_new / errors_summary_old * 100
                logging.debug("[%s/%s] By reflection of %s error changed from %.3f to %.3f, i.e. %.2f%%"
                    % (generation, max_generations, for_helic, errors_summary_old, errors_summary_new, percentage))
                errors_summary_old = errors_summary_new
            else:
                a_layout[for_helic] = a_layout[for_helic]-difangle
                for i in c_layout[for_helic]:
                    c_layout[for_helic][i] = (c_layout[for_helic][i][0], -c_layout[for_helic][i][1])
                    ca_layout[for_helic][i] = 360 - ca_layout[for_helic][i]

        if generation % CHECK_EVERY == 0:
            # If there was no significant improvement in last 10k iterations, there is no reason to continue
            esn, esnm = count_error(layout, a_layout, c_layout, ca_layout, \
                 sses_sizes, relevant_sses, relevant_sses, movable_sses,  sses_distance, start_layout, flags)
            errors_summary = esn+esnm

            if error_checkpoint != None and \
                (error_checkpoint - errors_summary) < (errors_summary * MINIMAL_IMPROVEMENT):
                logging.debug(f"start angles layout is {a_layout}")
                logging.debug(f"start ca layout is {ca_layout}")
                logging.debug(f"start layout {c_layout}")

                logging.debug("In last %s generations (step %s), total error decreased from %.2f to %.2f, which is below minimal required %.2f%%, so returning what we have" \
                    % (CHECK_EVERY, generation, error_checkpoint, errors_summary, MINIMAL_IMPROVEMENT*100))
                return layout, a_layout, c_layout, ca_layout
            error_checkpoint = errors_summary

        if (STEP_SIZE / coef) > 2:
            if (len(movable_clusters)<6) and (generation % 1000 == 0) \
            or (len(movable_clusters)<12) and (generation % 2000 == 0) \
            or (len(movable_clusters)<25) and (generation % 4000 == 0) \
            or (len(movable_clusters)<100) and (generation % 6000 == 0) \
            or ( generation % 8000 == 0):
                coef = coef * 2

        if generation >= max_generations:
            logging.debug("Finished after %s generations, maybe with more generations we can get better results? {max_generations}")
            return layout, a_layout, c_layout, ca_layout


def count_sse_protein_coords_based_on_cluster_coords(layout, a_layout, c_layout, ca_layout, angles):

    sses_layout = {}
    sses_angles_layout = {}
    for d in layout:
       sses_layout[d] = {}
       sses_angles_layout[d] = {}
       for c in layout[d]:
          sses_layout[d][c] = {}
          sses_angles_layout[d][c] = {}
          angle = a_layout[d][c]/180 * math.pi
          for sse in c_layout[d][c].keys():
            # count down sse1 coords based on clusters coords, cluster layout and cluster_layout_angles
            x = c_layout[d][c][sse][0]
            y = c_layout[d][c][sse][1]
            ssedif = math.sqrt(x*x+y*y)
            if y == 0:
               beta = 1
            else:
               beta = math.atan(x/y)
            if y<0:
                coord_sse1_x = layout[d][c][0] - (ssedif*math.sin(angle+beta))
                coord_sse1_y = layout[d][c][1] - (ssedif*math.cos(angle+beta))
            else:
                coord_sse1_x = layout[d][c][0] + (ssedif*math.sin(angle+beta))
                coord_sse1_y = layout[d][c][1] + (ssedif*math.cos(angle+beta))

            sses_layout[d][c][sse] = (coord_sse1_x, coord_sse1_y)
            sses_angles_layout[d][c][sse] = (a_layout[d][c] + ca_layout[d][c][sse]) % 360
    if angles == 0:
       return sses_layout
    else:
       return sses_layout, sses_angles_layout


def count_protein_error(layout, a_layout, c_layout, ca_layout, \
             for_domain, for_chain, sses_distance, sses_sizes):
   err_summ = 0
   sses_layout = count_sse_protein_coords_based_on_cluster_coords(\
         layout, a_layout, c_layout, ca_layout, 0)
   for d in layout:
       for c in layout[d]:
          if (d != for_domain) or (c != for_chain):
             for sse1 in sses_layout[d][c]:
                if sse1 in sses_distance[d][c][for_domain][for_chain]:
                   for sse2 in sses_layout[for_domain][for_chain]:
                       if sse2 in sses_distance[d][c][for_domain][for_chain][sse1]:
                           dif2d = math.dist(sses_layout[d][c][sse1],sses_layout[for_domain][for_chain][sse2])
                           dif3d = sses_distance[d][c][for_domain][for_chain][sse1][sse2]
                           err_summ = err_summ + abs(dif3d-dif2d) * ( sses_sizes[d][c][sse1] +sses_sizes[for_domain][for_chain][sse2])
   return err_summ

def count_all_protein_error(layout, a_layout, c_layout, ca_layout, \
             sses_distance, sses_sizes):
   err_all_sum = 0
   for d in layout:
       for c in layout[d]:
          err_all_sum = err_all_sum + count_protein_error(layout, a_layout, c_layout, ca_layout, \
             d, c, sses_distance, sses_sizes)
   return err_all_sum


def modify_protein_layout_via_evolution(layout, a_layout, c_layout, ca_layout, \
        sses_distance, sses_sizes, tmp_dir):

    nchain = 0
    for i in layout:
       nchain = nchain + len(i)

    if (nchain == 0):
        return layout, a_layout, c_layout, ca_layout

    if (nchain<6):
        STEP_SIZE = 25
    else:
        STEP_SIZE = 50
    ANGLES_SIZE = 45#180#45
    if nchain<6:
        CHECK_EVERY = 500
    elif nchain<12:
        CHECK_EVERY = 1000
    else:
        CHECK_EVERY = 2000   # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
    logging.debug(f"domain number {len(layout)} chain number {nchain} CHECK_EVERY {CHECK_EVERY}")
    MINIMAL_IMPROVEMENT = 0.01   # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout
    max_generations = 60000

    generation = 0
    error_checkpoint = None
    coef = 1

    ne_domain = []
    for i in list(layout):
       if layout[i] != {}:
          ne_domain.append(i)

    if ((len(ne_domain) == 1)) and (len(layout[ne_domain[0]])== 1):
        return layout, a_layout, c_layout, ca_layout

    while True:
        generation += 1
        # choose the moved sse
        for_domain = random.choice(list(ne_domain))
        for_chain = random.choice(list(layout[for_domain]))
        errors_summary_old = count_protein_error(layout, a_layout, c_layout, ca_layout, \
             for_domain, for_chain, sses_distance, sses_sizes)

        # set the new coordinates
        dx = random.random() * (STEP_SIZE/coef) - STEP_SIZE/(2*coef)
        dy = random.random() * (STEP_SIZE/coef) - STEP_SIZE/(2*coef)
        layout[for_domain][for_chain] = (layout[for_domain][for_chain][0]+dx, layout[for_domain][for_chain][1]+dy)
        errors_summary_new = count_protein_error(layout, a_layout, c_layout, ca_layout, \
             for_domain, for_chain, sses_distance, sses_sizes)
        if errors_summary_new < errors_summary_old:
            percentage = 100 - errors_summary_new / errors_summary_old * 100
            print("DEBUG: [%s/%s] By move of %s %s(x:%.3f,y:%.3f - z %s)(x,y:%.3f,%.3f) error changed from %.3f to %.3f, i.e. %.2f%% "
                % (generation, max_generations, for_domain, for_chain, layout[for_domain][for_chain][0], \
                layout[for_domain][for_chain][1], STEP_SIZE/(2*coef), layout[for_domain][for_chain][0]-dx, \
                layout[for_domain][for_chain][1]-dy, errors_summary_old, errors_summary_new, percentage))
            errors_summary_old = errors_summary_new
        else:
            layout[for_domain][for_chain] = (layout[for_domain][for_chain][0]-dx, layout[for_domain][for_chain][1]-dy)

        # now we try to rotate nontrivial clusters to improve their orientation
        if (len(layout[for_domain][for_chain]) > 1):
            counter = 1
            while counter < 3:
               counter = counter + 1
               difangle = random.random() * (ANGLES_SIZE/coef) - ANGLES_SIZE/(2*coef)
               a_layout[for_domain][for_chain] = (a_layout[for_domain][for_chain] + difangle ) % 360
               errors_summary_new = count_protein_error(layout, a_layout, c_layout, ca_layout, \
                  for_domain, for_chain, sses_distance, sses_sizes)

               if errors_summary_new < errors_summary_old:
                   percentage = 100 - errors_summary_new / errors_summary_old * 100
                   print("DEBUG: [%s/%s] By rotation of %s %s ( %.3f from %.3f  dif %.3f (max %s) ) error changed from %.3f to %.3f, i.e. %.2f%%"
                       % (generation, max_generations, for_domain, for_chain, a_layout[for_domain][for_chain], \
                       (a_layout[for_domain][for_chain] - difangle ) % 360, difangle, ANGLES_SIZE/(2*coef),
                       errors_summary_old, errors_summary_new, percentage))
                   errors_summary_old = errors_summary_new
               else:
                   a_layout[for_domain][for_chain] = (a_layout[for_domain][for_chain] - difangle ) % 360

        # we have to try reflection (if the cluster is nontrivial and <500 generations
        if (len(layout[for_domain][for_chain]) > 1) and (generation < 500):
            for i in c_layout[for_domain][for_chain]:
                c_layout[for_domain][for_chain][i] = (c_layout[for_domain][for_chain][i][0], -c_layout[for_domain][for_chain][i][1])
                ca_layout[for_domain][for_chain][i] = 360 - ca_layout[for_domain][for_chain][i]
            difangle = random.random() * 360 - 180
            a_layout[for_domain][for_chain] = a_layout[for_domain][for_chain]+difangle
            errors_summary_new = count_protein_error(layout, a_layout, c_layout, ca_layout, \
                 for_domain, for_chain, sses_distance, sses_sizes)

            if errors_summary_new < errors_summary_old:
                percentage = 100 - errors_summary_new / errors_summary_old * 100
                print("DEBUG: [%s/%s] By reflection of %s %s error changed from %.3f to %.3f, i.e. %.2f%%"
                    % (generation, max_generations, for_domain, for_chain, errors_summary_old, errors_summary_new, percentage))
                errors_summary_old = errors_summary_new
            else:
                a_layout[for_domain][for_chain] = a_layout[for_domain][for_chain]-difangle
                for i in c_layout[for_domain][for_chain]:
                    c_layout[for_domain][for_chain][i] = (c_layout[for_domain][for_chain][i][0], \
                           -c_layout[for_domain][for_chain][i][1])
                    ca_layout[for_domain][for_chain][i] = 360 - ca_layout[for_domain][for_chain][i]


        if generation % CHECK_EVERY == 0:
            # If there was no significant improvement in last 10k iterations, there is no reason to continue
            errors_summary = count_all_protein_error(layout, a_layout, c_layout, ca_layout, \
                 sses_distance, sses_sizes)
            if error_checkpoint != None and \
                ((errors_summary == 0) or (error_checkpoint - errors_summary) < (errors_summary * MINIMAL_IMPROVEMENT)):
                return layout, a_layout, c_layout, ca_layout
            error_checkpoint = errors_summary

        if (STEP_SIZE / coef) > 2:
            if (nchain<6) and (generation % 1000 == 0) \
            or (nchain<12) and (generation % 2000 == 0) \
            or ( generation % 4000 == 0):
                coef = coef * 2

        if generation >= max_generations:
            print ("DEBUG: Finished after %s generations, maybe with more generations we can get better results?" % max_generations)
            return layout, a_layout, c_layout, ca_layout




