#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math
import random
import time
import tempfile
import logging
import utils
from collections import OrderedDict
import utils_matrix2layout
import subprocess


A_MAX = 30

def count_2d_angle(v1, v2):
    angle = abs(v1 - v2)
    # FIXME: Is this safe transformation?
    if (angle > 180):
        angle = 360 - angle
    return angle


def count_measured_angles_error(angle, dangle, sses_sizes, relevant_sses, selected_sses):
    error = 0
    for g in selected_sses:
        if len(g.split(":"))==1:
            for h in relevant_sses:
                if (h != g) and (len(g.split(":")) == 1) and (h,g) in angle:
                    error = error + abs(angle[h,g] - \
                        count_2d_angle(dangle[h], dangle[g]))*(sses_sizes[h]+sses_sizes[g])
    return error


def get_fake_angles_layout(clusters, start_layout, relevant_sses):
    """
    Just returns layout where SSEs are placed with angles [0] for all
    """
    layout = OrderedDict()
    for c in clusters:
        sses = c.split(":")
        is_in = False
        for i in sses:
            if i in relevant_sses:
                is_in = True
            # is_in ..  cluster is in relevant set
        if is_in == False:
            continue

        if len(sses) == 1:
            # if it is trivial cluster - look to start layout
            if sses[0] in start_layout:
                layout[c] = start_layout[sses[0]]
            else:
                layout[c] = 0
        else:
            layout[c] = 0
    return layout


def count_sse_angles_based_on_cluster_coords(a_layout, ca_layout):
    sses_angles_layout = OrderedDict()

    for cl in a_layout:
        angle = a_layout[cl] / 180 * math.pi
        for sse in cl.split(":"):
            sses_angles_layout[sse] = (a_layout[cl] + ca_layout[cl][sse])
    return sses_angles_layout


def move_angles_error(dangle, start_angles, sses_sizes, relevant_sses, selected_sses):
    m_error = 0
    for i in selected_sses:
        if (i in start_angles) and (i in dangle):
            dif = abs(start_angles[i]-dangle[i])
            if dif > 180:
               dif = 360 - 180
            m_error = m_error + dif * len(relevant_sses) * sses_sizes[i] / 2
    return m_error


def count_angles_error(a_layout, ca_layout, sses_sizes, relevant_sses, \
    for_helic, sses_angles, start_angles):

    sses_angles_layout = count_sse_angles_based_on_cluster_coords(\
         a_layout, ca_layout)

    if (type(for_helic) != list):
       fh = {for_helic}
    else:
       fh = for_helic

    error = count_measured_angles_error(sses_angles, \
               sses_angles_layout, sses_sizes, relevant_sses, fh)
    m_error = move_angles_error(sses_angles_layout, start_angles, sses_sizes,  relevant_sses, fh)
    return error + m_error


def final_angles_rotation(a_layout, relevant_sses, layout, \
       c_layout, ca_layout, sses_sizes, sses_positions):
    # now we have to rotate all angles to some base rotation
    # if there are clusters (fixed == 1) we can rotate only to opposit side TODO - is it necessary? van we do that better way?
    # if there is no cluster, we can rotate to arbitrarry angles

    if (len(relevant_sses)) < 2:
        # too small to do any rotation
        return a_layout, ca_layout

    sses = []
    for i in relevant_sses:
       for j in (i.split(":")):
          if j in sses_positions:
              sses.append(j)
    number = len(sses)

    sses_layout, sses_a_layout = utils_matrix2layout.count_sse_coords_based_on_cluster_coords(\
          layout, a_layout, c_layout, ca_layout)

    # find centers
    centre_2d = [0, 0]
    centre_3d = [0, 0, 0]
    for i in sses:
       centre_2d = [centre_2d[0] + sses_layout[i][0], centre_2d[1] + sses_layout[i][1]]
       centre_3d = [centre_3d[0] + ((sses_positions[i][0][0]+sses_positions[i][1][0]) / 2),
                centre_3d[1] + ((sses_positions[i][0][1]+sses_positions[i][1][1]) / 2),
                centre_3d[2] + ((sses_positions[i][0][2]+sses_positions[i][1][2]) / 2)]
    centre_2d = [centre_2d[0] / number, centre_2d[1] / number]
    centre_3d = [centre_3d[0] / number, centre_3d[1] / number, centre_3d[2] / number]

    # find best sses (the biggest one)
    size = -1
    best_sse = ""
    for i in sses:
        if sses_sizes[i] > size:
           best_sse = i
           size = sses_sizes[i]

    # find the positions of the best_sse
    best_sses_3d_start = sses_positions[best_sse][0]
    best_sses_3d_center = [sses_positions[best_sse][0][0]+sses_positions[best_sse][1][0]/2, \
        sses_positions[best_sse][0][1]+sses_positions[best_sse][1][1]/2, \
        sses_positions[best_sse][0][2]+sses_positions[best_sse][1][2]/2]
    best_sses_3d_end = sses_positions[best_sse][1]
    best_sses_2d_start = (sses_layout[best_sse][0]-math.sin(sses_a_layout[best_sse]/180*math.pi), sses_layout[best_sse][1]-1*math.cos(sses_a_layout[best_sse]/180*math.pi))
    best_sses_2d_center = sses_layout[best_sse]
    best_sses_2d_end = (sses_layout[best_sse][0]+math.sin(sses_a_layout[best_sse]/180*math.pi), sses_layout[best_sse][1]+math.cos(sses_a_layout[best_sse]/180*math.pi))

    if len(sses) != len(relevant_sses):
       if ((math.dist(centre_3d, best_sses_3d_start) > math.dist(centre_3d, best_sses_3d_end)) and \
          ((math.dist(centre_2d, best_sses_2d_start) < math.dist(centre_2d, best_sses_2d_end)))) or \
          ((math.dist(centre_3d, best_sses_3d_start) <  math.dist(centre_3d, best_sses_3d_end)) and \
          ((math.dist(centre_2d, best_sses_2d_start) > math.dist(centre_2d, best_sses_2d_end)))):
            for i in a_layout:
                a_layout[i] = (a_layout[i]+180) % 360
    else:
        best_sses_3d_center = [(best_sses_3d_start[0] + best_sses_3d_end[0])/2, (best_sses_3d_start[1] + best_sses_3d_end[1])/2, (best_sses_3d_start[2] + best_sses_3d_end[2])/2]
        alpha_3d = utils.count_angle_in_triangle(best_sses_3d_center, best_sses_3d_start, centre_3d)
        alpha_2d = utils.count_angle_in_triangle_2d(best_sses_2d_center, best_sses_2d_start, centre_2d)
        for i in a_layout:
            a_layout[i] = (a_layout[i]-alpha_2d+alpha_3d) % 360
    return a_layout, ca_layout

def modify_angles_layout_via_evolution(a_layout, ca_layout, \
        sses_angles, sses_sizes, movable_clusters, fixed_clusters, \
        start_angles, layout, c_layout, sses_positions): #, layout3d):

    ANGLES_SIZE = 45#180
    if (len(movable_clusters)<6):
        CHECK_EVERY = 1000
    elif (len(movable_clusters)<12):
        CHECK_EVERY = 2000
    else:
        CHECK_EVERY = 3000

    # now we have to rotate trivial clusters to inmploce their angles 
    # nontrivial have the firm angle
    movable_cl = []
    relevant_sses = []
    for cl in movable_clusters:
        if (len(cl.split(":")) == 1):
            movable_cl.append(cl)
        for i in (cl.split(":")):
            relevant_sses.append(i)
    logging.debug(f"movable_cluster {movable_cl} relevant_sses {relevant_sses}")

    if (movable_cl == []):
        return a_layout, ca_layout

    if False:
      # this is not speedup - do not used now
      #if STEP_SIZE >0.01:
      # movable_cl,  sses_angles, start_Angles, sses_positions
      # nacteno: a_layout, ca_layout, layout, c_layout, movable_clusters,sses_sizes, relevant_sses,
      # nenacteno:
      aux1 = {}
      for sse_pair in sses_angles.keys():
          if sse_pair[0] not in aux1:
              aux1[sse_pair[0]] = {}
          if sse_pair[1] not in aux1[sse_pair[0]]:
             aux1[sse_pair[0]][sse_pair[1]] = sses_angles[sse_pair]

      with tempfile.TemporaryDirectory() as tmpdir:
            tmpname = f"{tmpdir}/aux"
            with open(tmpname, "w") as tmp:
               aux = OrderedDict(sses_sizes)
               tmp.write((f"movable_clusters {movable_clusters}\n"))
               tmp.write((f"a_layout {a_layout}\n"))
               tmp.write((f"ca_layout {ca_layout}\n"))
               tmp.write((f"sses_sizes {aux}\n"))
               tmp.write((f"relevant_sses {relevant_sses}\n"))
               tmp.write((f"movable_cl {movable_cl}\n"))
               tmp.write((f"CHECK_EVERY {CHECK_EVERY}\n"))

               tmp.write((f"sses_angles {aux1}\n")) 
               # aux transformed to  // read line: sses_distance {'H3': {'H3': 0.0, 'H4': 10.449864918385309, 'H5': 10.327968099413358, 'H
               # now it is sses_angles {('E2', 'E2'): 0, ('E2', 'E3'): 73, ('E2', 'E4'): 115, ('E2', 'H5'): 79, ('E2', 'H7'): 129, ('E2', 'E8'): 40, ('E2', 'E9'): 120, ('E2', 'H10'): 103, ('E2', 'E11'): 3, ('E2', 'H12'): 94, ('E2', 'H16'): 94, ('E2', 'H17'): 39, ('E2', 'H19'): 119, ('E2', 'H21'): 102, ('E2', 'H22'): 122, ('E2', 'H23'): 27, ('E2', 'H25'): 116, ('E2', 'E27'): 172, ('E2', 'E28'): 70, ('E2', 'H29'): 109, ('E2', 'H31'): 69, ('E2', 'H33'): 88, ('E2', 'H35'):

               tmp.write((f"start_angles {OrderedDict(start_angles)}\n"))
               tmp.flush()
               tmp.close()

            command = f"../gcc/angles {tmpdir}/aux {tmpdir}/in"
            print("command\n", command)
            p = subprocess.Popen(command, shell=True) #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            r = p.wait()

            with open(f"{tmpdir}/in", "r") as fd:
                line = fd.readline();
                while line:
                    sline=line.split(",")
                    if sline[0] == "angle_layout":
                        a_layout[sline[1]] = float(sline[2])
                    if sline[0] == "angle_cluster_layout":
                        for cl in c_layout:
                            for sse in c_layout[cl]:
                                if sline[2] == sse:
                                    ca_layout[cl][sse] = float(sline[3])
                    line = fd.readline()

            if (fixed_clusters == set()) and (start_angles == {}) :
                return final_angles_rotation(a_layout, movable_clusters, \
                   layout, c_layout, ca_layout, sses_sizes, sses_positions)
            else:
                return a_layout, ca_layout

    else:
      # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
      MINIMAL_IMPROVEMENT = 0.01  # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout
      max_generations = 60000

      generation = 0
      error_checkpoint = None
      coef = 1
      while True:
        generation += 1
        # choose the moved sse
        for_helic = random.choice(list(movable_cl))
        errors_summary_old = count_angles_error(a_layout, ca_layout, sses_sizes, relevant_sses, \
             for_helic, sses_angles, start_angles)

        difangle = random.random() * (ANGLES_SIZE/coef) - ANGLES_SIZE/(2*coef)
        a_layout[for_helic] = (a_layout[for_helic] + difangle ) % 360
        errors_summary_new = count_angles_error(a_layout, ca_layout, \
             sses_sizes, relevant_sses, for_helic, sses_angles, start_angles)

        if errors_summary_new < errors_summary_old:
            percentage = 100 - errors_summary_new / errors_summary_old * 100
            logging.debug("[%s/%s] By rotation of %s (%s ->%s %s ) error changed from %.3f to %.3f, i.e. %.2f%%"
                   % (generation, max_generations, for_helic, a_layout[for_helic]-difangle, difangle, ANGLES_SIZE/(2*coef),
                   errors_summary_old, errors_summary_new, percentage))
            errors_summary_old = errors_summary_new
        else:
            a_layout[for_helic] = (a_layout[for_helic] - difangle ) % 360

        if generation % CHECK_EVERY == 0:
            # If there was no significant improvement in last 10k iterations, there is no reason to continue
            errors_summary = count_angles_error(a_layout, ca_layout, \
                 sses_sizes, relevant_sses, relevant_sses, sses_angles, start_angles)

            if error_checkpoint != None and \
                (error_checkpoint - errors_summary) < (errors_summary * MINIMAL_IMPROVEMENT):
                logging.debug("In last %s generations (step %s) , total error decreased from %.2f to %.2f, which is below minimal required %.2f%%, so returning what we have" \
                    % (CHECK_EVERY, generation, error_checkpoint, errors_summary, MINIMAL_IMPROVEMENT*100))
                if (fixed_clusters == set()) and (start_angles == {}) :
                    return final_angles_rotation(a_layout, movable_clusters, \
                         layout, c_layout, ca_layout, sses_sizes, sses_positions)
                else:
                    return a_layout, ca_layout
            error_checkpoint = errors_summary

        if (ANGLES_SIZE / coef) > 2:
            if (len(movable_clusters)<6) and (generation % 500 == 0) \
            or (len(movable_clusters)<12) and (generation % 1000 == 0) \
            or ( generation % 2000 == 0):
                coef = coef * 2

        if generation >= max_generations:
            logging.debug("Finished after %s generations, maybe with more generations we can get better results?" % max_generations)
            if (fixed_clusters == set()) and (start_angles == {}):
               return final_angles_rotation(a_layout, movable_clusters, \
                      layout, c_layout, ca_layout, sses_sizes, sses_positions)
            else:
               return a_layout, ca_layout
