#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math
import random
import logging
from collections import OrderedDict
import utils
COEF = 0.2  # TODO:


def get_distance_matrix_2d(layout, only_for=None, selected_cluster=None):
    """
    Compute distance matrix of given layout
      layout ... {'A': (x_a, y_a), 'B': (x_b, y_b)...}
      only_for ... if set to some helic letter, only compute distances for
                   that helic,
                   if set to list only compute distances for given list of helices,
                   if set to None, compute for each
    """
    if only_for:
        if not isinstance(only_for, list):
           only_for = [only_for]
    else:
       only_for = layout.keys()

    if selected_cluster:
       selected = selected_cluster
    else:
       selected = layout.keys()

    out = {}
    for key1 in only_for:
        out[key1] = {}
        for key2 in selected:
            out[key1][key2] = math.sqrt(
                (layout[key1][0]-layout[key2][0])**2
                + (layout[key1][1]-layout[key2][1])**2)
    return out


def get_measured_error_2d(a, b, sizes, keys=None, start_clusters=None):
    """
    Compute sum of absolute values of differences in two provided dictionaries
      a ... {'C': {'A': 5, 'B': 13}, 'A': {...}, 'B': {...}}
      b ... {'C': {'A': 4, 'B': 12}, 'A': {...}, 'B': {...}}
      keys ... ['A', 'B', 'C'], defaults to keys from a
    """
    if keys is None:
        keys = a.keys()
    if start_clusters is None:
        start_clusters = keys
    summary = 0
    for key1 in keys:
        if key1 in a.keys():
            for key2 in start_clusters:
                if key1 != key2 and key2 in a[key1]:
                    min_d = min([a[key1][key2], b[key1][key2]])
                    if (min_d <= 40):
                       summary += abs(a[key1][key2] - b[key1][key2]) * (sizes[key1] + sizes[key2]) * (1-(min_d/40))
    return float(summary)/2


def circle(centerx, centery, d, number, counter):
   x = centerx + d * math.cos(2 * math.pi / number * counter)
   y = centery + d * math.sin(2 * math.pi / number * counter)
   return (x, y)


def iterate_cluster_evolution(layout, relevant_helices, distances, sizes,
         limit_max, nearby, minimal_improvement, no_start, chain, domain):

    CHECK_EVERY = 1000   # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
    MAX_GENERATIONS = 50000

    generation = 0
    coef = 1
    error_checkpoint = None
    if len(relevant_helices) == 0:
        return layout
    logging.debug(f"Evolution with this set of relevant hel {limit_max} {relevant_helices}")
    distances_layout = get_distance_matrix_2d(layout)
    errors_summary = get_measured_error_2d(distances, distances_layout, sizes, relevant_helices)
    no_start_num = len(no_start)
    while True:
        generation += 1
        if generation % CHECK_EVERY == 0:
            # If there was no significant improvement in last 10k iterations, there is no reason to continue
            if (error_checkpoint != None and
                (error_checkpoint - errors_summary) < (errors_summary * minimal_improvement)) \
                or ((error_checkpoint != None) and (error_checkpoint - errors_summary) < 1):
                logging.debug(f"cluster {generation} In the last {CHECK_EVERY} generations, total error decreased from {error_checkpoint}, {errors_summary} which is less then  {minimal_improvement} {limit_max}")
                return layout
            error_checkpoint = errors_summary

        if generation >= MAX_GENERATIONS:
            logging.debug (f"Finished after {MAX_GENERATIONS} generations, maybe with more generations we can get better results ({limit_max})?")
            return layout

        layout_new = layout.copy()

        if (generation < no_start_num * 20):
            for_helic = random.choice(no_start)
        else:
            # choose the moved sse
            for_helic = random.choice(relevant_helices)
        # set the new coordinates
        difx = random.random() * limit_max/coef * 2 - limit_max/coef
        dify = random.random() * limit_max/coef * 2 - limit_max/coef
        x, y = layout_new[for_helic]
        layout_new[for_helic] = (x+difx, y+dify)
        dmax = find_max_distance_in_cluster_sses(layout_new, for_helic, nearby, chain, domain)
        if (dmax < limit_max):

            distances_layout = get_distance_matrix_2d(layout_new)
            errors_summary_new = get_measured_error_2d(distances, distances_layout, sizes, relevant_helices)
            if errors_summary_new < errors_summary:
                percentage = 100 - errors_summary_new / errors_summary * 100
                logging.debug("DEBUG cl%s: [%s/%s] By move of %s error changed from %.3f to %.3f, i.e. %.2f%% (%.2f, %.2f, ((%.2f, %.2f)) limit, max %s)"
                    % (limit_max, generation, MAX_GENERATIONS, for_helic,
                    errors_summary, errors_summary_new, percentage, x, y,
                    difx, dify,  limit_max/coef))
                layout = layout_new
                errors_summary = errors_summary_new
            if ((len(relevant_helices) > 5) and (generation % 2000 == 0)) or \
               ((len(relevant_helices) <= 5) and (generation % 500 == 0)):
               coef = coef * 2


def get_1d_error(layout, for_helic, bcc_near, chain, domain):
   err = 0
   for sse in bcc_near(for_helic):
       # this if is necessary for per protein stuff, it is possible that the cluster contains stuff from more chains
        if sse[0] in layout:
            dest = abs(layout[for_helic][1] - layout[sse[0]][1])
            err = err + (2.5-dest)*(2.5-dest)
   return err


def get_1d_error_all(layout, bcc_near, chain, domain):
   err = 0
   for sse1 in layout.keys():
      for sse2 in bcc_near(sse1):
       # this if is necessary for per protein stuff, it is possible that the cluster contains stuff from more chains
         if sse2[0] in layout:
             dest = abs(layout[sse1][1] - layout[sse2[0]][1])
             err = err + (2.5-dest)*(2.5-dest)
   return err


def find_max_distance_in_cluster_sses(cluster_layout, sses, nearby, chain, domain):
    max_dist = 0
    for j in nearby(sses):
        # this if is necessary for per protein stuff, it is possible that the cluster contains stuff from more chains
        if j[0] in cluster_layout:
            dist = math.dist(cluster_layout[sses], cluster_layout[j[0]])
            if dist > max_dist:
                max_dist = dist

    return max_dist


def find_max_distance_in_cluster(cluster_layout, nearby, chain, domain):
    max_dist = 0
    for i in cluster_layout:
        for j in nearby(i):
            # this if is necessary for per protein stuff, it is possible that the cluster contains stuff from more chains
            if j[0] in cluster_layout:
                dist = math.dist(cluster_layout[i], cluster_layout[j[0]])
                if dist > max_dist:
                    max_dist = dist
    return max_dist


def put_clusters_layout_to_given_interval(cluster_layout, ratio):
    centerx = 0
    centery = 0
    number = 0
    for i in cluster_layout:
        number = number + 1
        centerx = centerx + cluster_layout[i][0]
        centery = centery + cluster_layout[i][1]
    # we have to move all sses 
    centerx = centerx/number
    centery = centery-number

    for i in cluster_layout:
        cluster_layout[i]=((cluster_layout[i][0]-centerx)*ratio, (cluster_layout[i][1]-centery)*ratio)

    return cluster_layout

def find_pair_layout(distance_matrix, relevant_clusters):
    cluster_layout = {}
    if relevant_clusters[0] in distance_matrix and relevant_clusters[1] in distance_matrix:
        d = distance_matrix[relevant_clusters[0]][relevant_clusters[1]]
        if (d<2):
           d = 2
        if (d>4):
           d = 4
    else:
        d = 3
    cluster_layout[relevant_clusters[0]] = (0, -d/2)
    cluster_layout[relevant_clusters[1]] = (0, d/2)
    return cluster_layout

def modify_cluster_layout_via_evolution(cluster_layout, \
        distances, sizes, nearby, relevant_clusters, no_start, chain, domain):
    """
    Try random move of random SSE in the 2D cluster layout and only keep these if it
    lead to a lower total error the layout. Repeat untill there is significant
    improvement in last 10k attempts.
    """
    relevant_clusters = list(relevant_clusters)

    logging.debug(f"relevant_clusters {relevant_clusters}")
    # Tune first level sec structures
    # if there is no start layout
    if (no_start != []):
        logging.debug(f"there is no start layout try all parts {len(relevant_clusters)} {cluster_layout}")
        # and cluster have > 15 sses go to 30 and 25
        if (len(relevant_clusters)>16):
            cluster_layout = iterate_cluster_evolution(cluster_layout, \
               relevant_clusters, distances, sizes, 30, nearby, 0.02, no_start, chain, domain)

            max_dist = find_max_distance_in_cluster(cluster_layout, nearby, chain, domain)
            # put sses to positions which meets the limitations 0-7 distances
            if (max_dist > 25):
                cluster_layout = put_clusters_layout_to_given_interval(cluster_layout, 25/max_dist)
                cluster_layout = iterate_cluster_evolution(cluster_layout, \
                    relevant_clusters, distances, sizes, 25, nearby, 0.02, no_start, chain, domain)

        # if there is 15 ... 6 clusters with no start point start iteration here
        max_dist = find_max_distance_in_cluster(cluster_layout, nearby, chain, domain)
        # put sses to positions which meets the limitations 0-18 distances
        if  ((len(relevant_clusters) <= 16) and (len(relevant_clusters) > 6)) or (max_dist > 18):
            if (max_dist > 18):
                cluster_layout = put_clusters_layout_to_given_interval(cluster_layout, 18/max_dist)
            cluster_layout = iterate_cluster_evolution(cluster_layout, \
                relevant_clusters, distances, sizes, 18, nearby, 0.02, no_start, chain, domain)

        max_dist = find_max_distance_in_cluster(cluster_layout, nearby, chain, domain)
        # put sses to positions which meets the limitations 0-13 distances
        if (max_dist > 13):
            cluster_layout = put_clusters_layout_to_given_interval(cluster_layout, 13/max_dist)
            cluster_layout = iterate_cluster_evolution(cluster_layout, \
                relevant_clusters, distances, sizes, 13, nearby, 0.02, no_start, chain, domain)

    max_dist = find_max_distance_in_cluster(cluster_layout, nearby, chain, domain)
    # if the cluster is bigger then 6 and there is complete
    # or between 6 and 3  --| both cases have 9A as a start point
    # or the distances are bigger than 9A
    if ((len(relevant_clusters) > 6) and (no_start == [])) or \
       ((len(relevant_clusters) <= 6) and (len(relevant_clusters) > 3)) or \
       (max_dist > 9):
       if (max_dist > 9):
           cluster_layout = put_clusters_layout_to_given_interval(cluster_layout, 9/max_dist)
           # put sses to positions which meets the limitations 0-9 distances
       cluster_layout = iterate_cluster_evolution(cluster_layout, \
               relevant_clusters, distances, sizes, 9, nearby, 0.02, no_start, chain, domain)

    max_dist = find_max_distance_in_cluster(cluster_layout, nearby, chain, domain)
    if ((len(relevant_clusters) <= 3) or (max_dist > 6)):
        if (max_dist > 6):
            cluster_layout = put_clusters_layout_to_given_interval(cluster_layout, 6/max_dist)
            # put sses to positions which meets the limitations 0-6 distances
        cluster_layout = iterate_cluster_evolution(cluster_layout, \
                    relevant_clusters, distances, sizes, 6, nearby, 0.01, no_start, chain, domain)

    # put sses to positions which meets the limitations 2-4 distances
    max_dist = find_max_distance_in_cluster(cluster_layout, nearby, chain, domain)
    if (max_dist > 4):
        cluster_layout = put_clusters_layout_to_given_interval(cluster_layout, 4/max_dist)
        cluster_layout = iterate_cluster_evolution(cluster_layout, \
                relevant_clusters, distances, sizes, 4, nearby, 0.01, no_start, chain, domain)
    return cluster_layout


def tune_cluster_layout_to_2_3_near_distances(layout, near, chain, domain):
    " finaly tune cluster layout to keep distances between 2-3, in one direction"

    CHECK_EVERY = 1000   # check every N generations if there was improvement at least MINIMAL_IMPROVEMENT
    MINIMAL_IMPROVEMENT = 0.02   # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout
    max_generations = 20000

    error_checkpoint = None
    generation = 0
    movable_helices = list(layout.keys())
    coef = 1
    logging.debug(f"2-3 evolution with this set of relevant, fixed and movable helices {movable_helices}")
    while True:
        generation += 1
        if generation % CHECK_EVERY == 0:
            # If there was no significant improvement in last 10k iterations, there is no reason to continue
            errors_summary = get_1d_error_all(layout, near, chain, domain)
            if (error_checkpoint != None and \
                (error_checkpoint - errors_summary) < (errors_summary * MINIMAL_IMPROVEMENT)) or \
                (error_checkpoint != None and (error_checkpoint - errors_summary) < 1) or \
                (errors_summary == 0):
                logging.debug ( f"{generation}  In the last {CHECK_EVERY} generations({generation}), total error decreased from {error_checkpoint} {errors_summary} which is less then  {MINIMAL_IMPROVEMENT}")

                return layout
            error_checkpoint = errors_summary

        if generation >= max_generations:
            logging.debug ("Finished after %s generations, maybe with more generations we can get better results?" % max_generations)
            return layout

        layout_new = layout.copy()
        for_helic = random.choice(list(movable_helices))
        # set the new coordinates
        dif = random.random() * 2/coef - 1/coef
        layout_new[for_helic] = (layout_new[for_helic][0], layout_new[for_helic][1]+dif)
        errors_old = get_1d_error(layout, for_helic, near, chain, domain)
        errors_new = get_1d_error(layout_new, for_helic, near, chain, domain)
        if errors_new < errors_old:

            percentage = 100 - errors_new / errors_old * 100
            logging.debug("cl(2-3): [%s/%s] move %s (%.3f, %.3f) ec from %.3f to %.3f, i.e. %.2f%%"
                % (generation, max_generations, for_helic, layout[for_helic][1], layout_new[for_helic][1], errors_old, errors_new, percentage))
            layout = layout_new

        if ((len(movable_helices) > 5) and (generation % 2000 == 0)) or \
           ((len(movable_helices) <= 5) and (generation % 500 == 0)):
           coef = coef * 2

    return layout


def set_start_layout(cluster_sses, start_layout):
   # set layout for all cluster_sses members - 
   # if theer is start layout use if, if it is not use (0,0)
   cluster_layout = {}
   no_start = []
   start = []
   center_layout = (0, 0)
   for i in cluster_sses:
      if i in start_layout:
          start.append(i)
          center_layout = (center_layout[0] + start_layout[i][0],
                           center_layout[1] + start_layout[i][1])
      else:
          # there is no start layout, thus we use (0,0)
          no_start.append(i)
          cluster_layout[i] = (0,0)
   if (start != []):
       center_layout = (center_layout[0]/len(start),
                        center_layout[1]/len(start))
       for i in start:
           # there is start layout, we use the distance from center of
           # all start_layout ones
           cluster_layout[i] = (start_layout[i][0]-center_layout[0], \
                   start_layout[i][1]-center_layout[1])
   return cluster_layout, no_start


def rotate_cluster_layout(cluster_layout, start, end):
    # rotate cluster layout to have start-end on (0,0, (0,dist(start,end))
    preplayout=OrderedDict()
    dx=cluster_layout[end][0]-cluster_layout[start][0]
    dy=cluster_layout[end][1]-cluster_layout[start][1]
    angle = math.pi/2- math.atan2(dy,dx)
    if angle > math.pi/2:
        angle = angle - math.pi
    for sse in cluster_layout:
        dsax = cluster_layout[sse][0]-cluster_layout[start][0]
        dsay = cluster_layout[sse][1]-cluster_layout[start][1]
        prepx=math.cos(angle)*dsax-math.sin(angle)*dsay
        prepy=math.sin(angle)*dsax+math.cos(angle)*dsay
        preplayout[sse]=(prepx,prepy)
    return preplayout


def center_cluster_layout(cluster_layout):
    sumx = 0
    sumy = 0
    sumc = 0
    for i in cluster_layout:
        sumc = sumc + 1
        sumx = sumx + cluster_layout[i][0]
        sumy = sumy + cluster_layout[i][1]
    sumx = sumx / sumc
    sumy = sumy / sumc
    for i in cluster_layout:
          cluster_layout[i] = (cluster_layout[i][0] - sumx, cluster_layout[i][1] - sumy)
    return cluster_layout


# depth first search - go through all path through the graph
def dfs(ret, start, bcc_near, buffer, level, chain, domain):
    for i in bcc_near(start):
        if i[0] not in buffer:
            bl = len(buffer)
            buffer.append(i[0])
            level[i[0]] = level[start]+1
            ret = dfs(ret, i[0], bcc_near, buffer, level, chain, domain)
            del(buffer[bl])
            del(level[i[0]])
        else:
            if (level[start]-level[i[0]]) > 3:
                return True
        if (ret == True):
           return ret
    return ret


# find whether there is a circle bigger then 3 sses
def cluster_is_circle(clusters_sses, bcc_near, chain, domain):
    start = clusters_sses[0]
    buffer = [start]
    level = {}
    level[start]=1
    ret = dfs(False, start, bcc_near, buffer, level, chain, domain)
    return ret


def compute_clusters_layout(clusters, distance_matrix, sizes, bcc_near, \
       relevant_sses, bcc_orientation, start_layout, start_angles, prime_sses, \
       chain, domain):
    """
    Compute clusters layout - in clusters with >1 sses (sheets),
    there have to be generated layout. The distances between sheets have to be <2 and >4
    based on beta connections
    """
    all_clusters_layout = OrderedDict()
    all_clusters_angles_layout = OrderedDict()

    for cluster in clusters:
        # add the first sse of the given layout
        clusters_sses = cluster.split(":")

        is_in = False
        for i in clusters_sses:
            if i in relevant_sses:
                is_in = True
        if not is_in:
            continue  # there is no need to continue cluster does not include any visible sse

        cluster_layout, no_start = set_start_layout(clusters_sses, start_layout)
        if (len(clusters_sses) == 2) and (no_start != []):
            # there are only two sses in cluster, there is no start layout this pair
            cluster_layout = find_pair_layout(distance_matrix,
                                              relevant_clusters=clusters_sses)

        if len(clusters_sses) > 2:
            cluster_layout = modify_cluster_layout_via_evolution(cluster_layout,
                                                                 distance_matrix,
                                                                 sizes,
                                                                 bcc_near,
                                                                 clusters_sses,
                                                                 no_start,
                                                                 chain,
                                                                 domain)

        if len(clusters_sses) > 1:
            if cluster_is_circle(clusters_sses, bcc_near, chain, domain) is True:
                clusters_sses.sort()
                n = 1
                while (n > len(clusters_sses)) or (clusters_sses[len(clusters_sses)-n] not in prime_sses):
                   n = n + 1
                if (n > len(clusters_sses)):
                   n = 1
                start = clusters_sses[len(clusters_sses)-n]   # TODO have to be the first (lex order)
                end = utils.find_most_distant_point(cluster_layout, clusters_sses, start)
            else:
                relevant_clusters_sses = []
                for i in clusters_sses:
                   if i in relevant_sses:
                      relevant_clusters_sses.append(i)
                start, end = utils.find_most_distant_pair(cluster_layout, relevant_clusters_sses)
            if start != end:
                cluster_layout = rotate_cluster_layout(cluster_layout, start, end)
                cluster_layout = tune_cluster_layout_to_2_3_near_distances(cluster_layout, bcc_near, chain, domain)

        all_clusters_layout[cluster] = center_cluster_layout(cluster_layout)
        start_sse = clusters_sses[0]

        dangle = {}
        if len(clusters_sses) > 1:
           if (start in cluster_layout) and (end in cluster_layout) and (cluster_layout[start][1] > cluster_layout[end][1]):
              aux = start
              start = end
              end = aux
           # it is not trivial cluster
           if (start in start_layout) and (end in start_layout) and (start != end) and (start_sse in start_layout):

               difx = (start_layout[start])[0]-(start_layout[end])[0]
               dify = (start_layout[start])[1]-(start_layout[end])[1]
               new_cl_end_layout = [(cluster_layout[start])[0]-difx,
                                    (cluster_layout[start])[1]-dify]
               angle = utils.count_angle_in_triangle_2d(new_cl_end_layout, cluster_layout[start], cluster_layout[end])
               if (difx < 0) and (dify < 0):
                   angle = - angle
               if (abs((start_angles[start_sse]+angle) % 360 - 270) < abs((start_angles[start_sse]+angle) % 360 - 90)):
                  value = 270
               else:
                  value = 90
           else:
               value = 90
           cl, sameset, opositeset = bcc_orientation(start_sse)
           dangle[start_sse] = value
           for h in sameset:
                  dangle[h] = value
           for h in opositeset:
                  dangle[h] = (value+180) % 360
        else:
           dangle = {}
           dangle[clusters_sses[0]] = 0
        all_clusters_angles_layout[cluster] = dangle

    logging.debug(f"all_clusters_layout {all_clusters_layout}")
    logging.debug(f"all_clustars_angles_layout {all_clusters_angles_layout}")
    return all_clusters_layout, all_clusters_angles_layout

