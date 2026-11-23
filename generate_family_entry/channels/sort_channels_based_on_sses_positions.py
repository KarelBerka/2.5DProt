#!/usr/bin/env python3
# go through unannotated channels and write which channel is the nearest one

import argparse
import json
import sys
import os
import numpy as np
import time
import math
sys.path.append("../..")
import utils_vector


MIN_RATIO = 0.2

MAX_RATIO = 3

def count_file_lines(fn):
    l_count = 0
    with open(fn, "r") as fr:
        line = fr.readline()
        while line:
            sline =line.split(" ")
            l_count += 1
            line = fr.readline()
    return l_count


def short_ok_long_by_lines(gdir):
    # filter too short < MIN_RATIO * avg
    # and long > MAX_RATIO * avg channels
    # ength defined by number of lines in channel files
    lines_sum = 0
    files_count = 0
    # at first set the average number of lines of channel files
    for file_name in os.listdir(gdir):
        if file_name.startswith("gnuplot-"):
            lines_sum += count_file_lines(f"{gdir}/{file_name}")
            files_count += 1
    # avg channel length
    avg = lines_sum/files_count
    short_channels = []
    long_channels = []
    ok_channels = []
    for file_name in os.listdir(gdir):
        if file_name.startswith("gnuplot-"):
             counter = count_file_lines(f"{gdir}/{file_name}")
             if counter < avg*MIN_RATIO:
                 short_channels.append(file_name)
             elif counter < avg*MAX_RATIO:
                 ok_channels.append(file_name)
             else:
                 long_channels.append(file_name)

    print(f"limits {avg*MIN_RATIO}, {avg*MAX_RATIO}")
    print(f"short {len(short_channels)}, ok {len(ok_channels)}, long {len(long_channels)}")
    return short_channels, ok_channels, long_channels


def short_ok_long_by_name(gdir):
    # filter too short < MIN_RATIO * avg
    # and long > MAX_RATIO * avg channels
    # length defined by channel name
    lines_sum = 0
    files_count = 0
    # at first set the average number of lines of channel files
    for file_name in os.listdir(gdir):
        if file_name.startswith("gnuplot-"):
            fsplit = file_name.split("_")
            print("fsplit is", fsplit)
            lines_sum += float(fsplit[-1])
            files_count += 1

    # avg channel length
    avg = lines_sum/files_count
    print(f"avg: {avg}, files_count {files_count}")
    short_channels = []
    long_channels = []
    ok_channels = []
    for file_name in os.listdir(gdir):
        if file_name.startswith("gnuplot-"):
             fsplit = file_name.split("_")
             counter = float(fsplit[-1])
             if counter < avg*MIN_RATIO:
                 short_channels.append(file_name)
             elif counter < avg*MAX_RATIO:
                 ok_channels.append(file_name)
             else:
                 long_channels.append(file_name)

    print(f"limits {avg*MIN_RATIO}, {avg*MAX_RATIO}")
    print(f"short {len(short_channels)}, ok {len(ok_channels)}, long {len(long_channels)}")
    return short_channels, ok_channels, long_channels


def get_sses_set(annot_fn, pcd):
    sses_list = []
    # TODO - filtr only main sses
    with open(f"{annot_fn}", "r") as fr:
         data = json.load(fr)
         for sse in data[pcd]["secondary_structure_elements"]:
             sses_list.append({"start": sse["start_vector"], "end": sse["end_vector"]})
#    print("--", sses_list)
    return sses_list


def transp(p):
    return np.array([p[0]], [p[1]], [p[2]])


def line_cross_triangle(triangle, line):
    p1 = np.array([triangle[0][0], triangle[0][1], triangle[0][2]]) # point 1
    p2 = np.array([triangle[1][0], triangle[1][1], triangle[1][2]]) # point 2
    p3 = np.array([triangle[2][0], triangle[2][1], triangle[2][2]]) # point 3

    s1 = np.array([line["start"][0], line["start"][1], line["start"][2]])
    s2 = np.array([line["end"][0], line["end"][1], line["end"][2]])

    p21 = np.subtract(p2,p1) # vector from p1 to p2
    p31 = np.subtract(p3,p1) # vector from p1 to p3

    triangle_norm = np.cross(p21, p31) # triangle norm

    N12 = np.cross((np.subtract(p2, p1)), triangle_norm) # normal vector to side p21
    N23 = np.cross((np.subtract(p3, p2)), triangle_norm) # normal vector to side p32
    N31 = np.cross((np.subtract(p1, p3)), triangle_norm) # normal vector to side p13

    D = -np.dot(triangle_norm,p1) # distance of plane N . x to point p1 (to define plane of triangle p1,p2,p3)

    t = -(D + np.dot(triangle_norm,s1))/ (np.dot(triangle_norm,(s2-s1)))# intersection of triangle and line

    if t < 0 or t > 1: # it is segment not line
        return False
    else:
        R = s1 + t * (s2 - s1) # intersection of tringle and segmet
        dist_r_n12 = np.dot((R-p1), N12)/np.linalg.norm(N12)
        dist_r_n23 = np.dot((R-p2), N23)/np.linalg.norm(N23)
        dist_r_n31 = np.dot((R-p3), N31)/np.linalg.norm(N31)
        if dist_r_n12 < 0 and dist_r_n23 < 0 and dist_r_n31 < 0:
            return True
        else:
            return False


def is_sse_between_channels(c1, c2, sse):
    # return True = c2 is not in cluster with c1
    ratio_c1 = count_file_lines(c1) /9
    ratio_c2 = count_file_lines(c2) /9

    c1_points = []
    with open(c1, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            if line_counter > point_counter * ratio_c1:

                sline = line.split(" ")
                c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                point_counter += 1
            line_counter += 1
            sline = line.split(" ")
            line = fr.readline()
        c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    while len(c1_points) < 10:
        c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])

    c2_points = []
    with open(c2, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            if line_counter > point_counter * ratio_c2:
                sline = line.split(" ")
                c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                point_counter += 1
            line_counter += 1
            sline = line.split(" ")
            line = fr.readline()
        c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    while len(c2_points) < 10:
        c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])

    ret = False
    dist = (c1_points[9][0]-c2_points[9][0])*(c1_points[9][0]-c2_points[9][0])+(c1_points[9][1]-c2_points[9][1])*(c1_points[9][1]-c2_points[9][1])+(c1_points[9][2]-c2_points[9][2])*(c1_points[9][2]-c2_points[9][2])

    # find out whether sse is not between c1 and c2 modelled by triangles
    for n in range(0,9):
        intersection = line_cross_triangle([c1_points[n], c1_points[n+1], c2_points[n]], sse)
        if intersection:
                return True

    for n in range(0,9):
        intersection = line_cross_triangle([c1_points[n+1], c2_points[n+1], c2_points[n]], sse)
        if intersection:
                return True

    return False


def is_near(c1, c2, sse):
    # return True = c2 is not in cluster with c1
    ratio_c1 = count_file_lines(c1) /9
    ratio_c2 = count_file_lines(c2) /9

    c1_points = []
    with open(c1, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            if line_counter > point_counter * ratio_c1:

                sline = line.split(" ")
                c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                point_counter += 1
            line_counter += 1
            sline = line.split(" ")
            line = fr.readline()
        c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    while len(c1_points) < 10:
        c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])

    c2_points = []
    with open(c2, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            if line_counter > point_counter * ratio_c2:
                sline = line.split(" ")
                c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                point_counter += 1
            line_counter += 1
            sline = line.split(" ")
            line = fr.readline()
        c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    while len(c2_points) < 10:
        c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    ret = False
    dist = (c1_points[9][0]-c2_points[9][0])*(c1_points[9][0]-c2_points[9][0])+(c1_points[9][1]-c2_points[9][1])*(c1_points[9][1]-c2_points[9][1])+(c1_points[9][2]-c2_points[9][2])*(c1_points[9][2]-c2_points[9][2])
    if dist > 169:
         #if the ends od channel is more than 15 A distant from the original
#         print("dale nez 15:", end = " ")
         return False
    else:
         return True


def is_start_near(c1, c2, sse):
    # return True = c2 is not in cluster with c1
    ratio_c1 = count_file_lines(c1) /9
    ratio_c2 = count_file_lines(c2) /9

    c1_points = []
    with open(c1, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            if line_counter > point_counter * ratio_c1:

                sline = line.split(" ")
                c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                point_counter += 1
            line_counter += 1
            sline = line.split(" ")
            line = fr.readline()
        c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    while len(c1_points) < 10:
        c1_points.append([float(sline[0]), float(sline[1]), float(sline[2])])

    c2_points = []
    with open(c2, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            if line_counter > point_counter * ratio_c2:
                sline = line.split(" ")
                c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                point_counter += 1
            line_counter += 1
            sline = line.split(" ")
            line = fr.readline()
        c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    while len(c2_points) < 10:
        c2_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
    ret = False
    dist = (c1_points[0][0]-c2_points[0][0])*(c1_points[0][0]-c2_points[0][0])+(c1_points[0][1]-c2_points[0][1])*(c1_points[0][1]-c2_points[0][1])+(c1_points[0][2]-c2_points[0][2])*(c1_points[0][2]-c2_points[0][2])
    if dist > 100:
         #if the ends od channel is more than 15 A distant from the original
#         print("dale nez 15:", end = " ")
         return False
    else:
         return True





def find_channels_from_the_same_cluster(start_channel, channels_set, channel_dir, sses_set, simple, ommit_list, limit_max, start):
    new_cluster = [start_channel]
    rest = []
    print(f"find_channels_from_the_same_cluster: start cluster generating with start {start_channel}", end=" ")
    counter = 0
    for c in channels_set:
        counter += 1
        print(f"{counter} len(new_cluster){len(new_cluster)} . {c} ", end = " ")
        if c in ommit_list:
            print(" - ommit")
            continue
        the_same_type = True
        if simple:
            for sse in sses_set:
                # thus there is no sses between this two channels
#                    print("tady", f"{channel_dir}/{cc}", f"{channel_dir}/{c}")
                if not is_near(f"{channel_dir}/{start_channel}", f"{channel_dir}/{c}", sse) or \
                   is_sse_between_channels(f"{channel_dir}/{start_channel}", f"{channel_dir}/{c}", sse) or \
                   (start and not  is_start_near(f"{channel_dir}/{start_channel}", f"{channel_dir}/{c}", sse)):
                    the_same_type = False
                    rest.append(c)
                    break
        else:
            counter2 = 0
            for cc in new_cluster:
                counter2 += 1
                # we want to fint the channel which is similar
                if c == cc:
                    # skip the start one
                    continue

                for sse in sses_set:
                    # thus there is no sses between this two channels
#                    print("tady", f"{channel_dir}/{cc}", f"{channel_dir}/{c}")
                    if not is_near(f"{channel_dir}/{cc}", f"{channel_dir}/{c}", sse) or is_sse_between_channels(f"{channel_dir}/{cc}", f"{channel_dir}/{c}", sse):
                        the_same_type = False
                        rest.append(c)
                        break
                if not the_same_type:
                    break
                if counter2 >100:
                    break

        if the_same_type:
            csplit = c.split("_")
            if float(csplit[-1])>float(limit_max):
                print(" prilis dlouhe", c)
            else:
            #    print("----", csplit[-1], limit_max,)
#            sys.exit()
                new_cluster.append(c)
                print("*  dobre")
        else:
            print("")
        if counter > 20000:
            print("")
            return new_cluster, rest

    return new_cluster, rest


def find_channels_from_the_same_cluster2(start_channel, channels_set, channel_dir, sses_set, simple, ommit_list, limit_max):
    # near and ommitd by the vprevious version
    new_cluster = [start_channel]
    rest = []
    print(f"find_channels_from_the_same_cluster: start cluster generating with start {start_channel}", end=" ")
    counter = 0
    for c in channels_set:
        counter += 1
        print(f"{counter} len(new_cluster){len(new_cluster)} . {c} ", end = " ")
        the_same_type = False
        if simple:
            for sse in sses_set:
                # thus there is no sses between this two channels
#                    print("tady", f"{channel_dir}/{cc}", f"{channel_dir}/{c}")
                a =  is_near(f"{channel_dir}/{start_channel}", f"{channel_dir}/{c}", sse)
                b = is_sse_between_channels(f"{channel_dir}/{start_channel}", f"{channel_dir}/{c}", sse)
                print("----", c, limit_max)
                sys.exit()
#                print(f"{a},{b}", end="  ")
                if not a:
                    the_same_type = False
                    rest.append(c)
                    break
                if b:
                    the_same_type = True
                #    rest.append(c)
                #    break
        else:
            counter2 = 0
            for cc in new_cluster:
                counter2 += 1
                # we want to fint the channel which is similar
                if c == cc:
                    # skip the start one
                    continue

                for sse in sses_set:
                    # thus there is no sses between this two channels
#                    print("tady", f"{channel_dir}/{cc}", f"{channel_dir}/{c}")
                    if is_sse_between_channels(f"{channel_dir}/{cc}", f"{channel_dir}/{c}", sse):
                        the_same_type = False
                        rest.append(c)
                        break
                if not the_same_type:
                    break
                if counter2 >100:
                    break

        if the_same_type:
            new_cluster.append(c)
            print("*")
        else:
            rest.append(c)
            print("")

        if counter > 20000:
            return new_cluster, rest

    return new_cluster, rest





# read arguments
parser = argparse.ArgumentParser(
    description='TODO')
parser.add_argument('--dir', required=True, action="store",
                    help="channel files directory (gnuplot-{family_nme}-{timestamp})")
parser.add_argument('--name', required=True, action="store",
                    help="name")
parser.add_argument('--annotation_dir', required=True, action="store",
                    help="annotation files directory generated-{family_nme}-{timestamp}/annotation")
parser.add_argument('--start_channel', action="store",
                    help="generate only cluster contaiting start channel (e.g. gnuplot-1spu03B.txt_29.14)")
parser.add_argument('--simple', action="store_true",
                    help="compare clusters only with the selected start representative")
parser.add_argument('--ommit', action="store",
                    help="list of ommited domains")
parser.add_argument('--limit_max', action="store",
                    help="maximal channel length")
parser.add_argument('--start', action="store_true",
                    help="maximal channel length")


args = parser.parse_args()

ch_list = {}
coef = 5
channels_clusters = {}

#short_ch, ok_ch, long_ch = short_ok_long_by_lines(args.dir)
short_ch, ok_ch, long_ch = short_ok_long_by_name(args.dir)

rest = ok_ch
counter = 1

ommit_list = []
if args.ommit:
    with open(args.ommit, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            sline = line.split("'")
            ssline = sline[1].split("/")
            ommit_list.append(ssline[2])
#            print("++++",ssline[2])
            line = fr.readline()
#print("lll", len(ommit_list))
#sys.exit()


while counter < 35 and len(rest) > 0:
    channels_set = rest
    if (args.start_channel):
        start_channel = args.start_channel
    else:
        start_channel = channels_set[0]
    if (args.start):
        start = True
    else:
        start = False
    if (args.simple):
        simple = True
    else:
        simple = False
    if (args.limit_max):
        limit_max = args.limit_max
    else:
        limit_max = 1000
    print("start_channel", start_channel, len(channels_set))
    aux = start_channel.split(".")[0][8:]
    pcd = f"{aux[0:4]}{aux[6:]}{aux[4:6]}"
    annot_fn = f"{args.annotation_dir}/{pcd}-annotated.sses.json"
    print(f"get_sses_set({annot_fn}, {pcd})")
    sses_set = get_sses_set(annot_fn, pcd)     # return sses set - its start and end coordinates

    c_name = f"C{counter}"
    channels_clusters[c_name], rest = find_channels_from_the_same_cluster(start_channel, channels_set, args.dir, sses_set, simple, ommit_list, limit_max, start)
    o_fn=f"lists/list_{args.name}_cl_{counter}"
    print("cluster", counter, "number", len(channels_clusters[c_name]))
    print("pisu asi sem:", o_fn)
    with open(o_fn, "w") as fw:
        for i in channels_clusters[c_name]:
            fw.write(f"'{args.dir}/{i}' lt rgb 'white'\n")
    if args.start_channel:
        counter = 36
    else:
        counter += 1

print("rest is ", len(rest))
