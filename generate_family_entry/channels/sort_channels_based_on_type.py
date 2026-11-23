#!/usr/bin/env python3
# go through unannotated channels and write which channel is the nearest one

import argparse
import json
import sys
import os
import numpy
import time
import math
sys.path.append("../..")
import utils_vector


def d(p1, p2):
#    print("pp", p1, p2)
    return math.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]))


def diff(n1, n2, ch_list):
    max_diff = 0
    for snap1 in ch_list[n1]:
        min_diff = 30
        for snap2 in ch_list[n2]:
            diff_p_to_p = d(snap1, snap2)
            if diff_p_to_p < min_diff:
                min_diff = diff_p_to_p
        if min_diff > max_diff:
            max_diff = min_diff
    return(max_diff)


def set_cluster(start, ch_rest, t):
    cl = [start]
#    ch_rest.remove(cl)
#    print("cl", cl, ch_rest)
    counter = 0
    while counter < 100:
        cc_min = 100
        counter += 1
        cc_c = None
        for i in ch_rest:
            c_min = 0
            for j in cl:
                if t[j][i] > c_min:
                    c_min = t[j][i]
            if c_min < cc_min:
                cc_min = c_min
                cc_c = i
#            print("vybral jsem i", i , "z", len(ch_rest))
#            print("vyslo", c_min, cc_min)

#        print("ii", i, cc_min)
        if cc_min < 20:
            cl.append(cc_c)
            ch_rest.remove(cc_c)
        else:
            return(cl)

# read arguments
parser = argparse.ArgumentParser(
    description='generate gnuplot input data for input channel annotation files from input directory')
parser.add_argument('--dir_gnuplot', required=True, action="store",
                    help="channel files directory")
parser.add_argument('--dir_channels', required=True, action="store",
                    help="channel files directory")

parser.add_argument('--name', required=True, action="store",
                    help="jmeno")

args = parser.parse_args()

c_types = {}
counter = 0
for file_name in os.listdir(args.dir_gnuplot):
    if file_name.startswith("gnuplot-"):
         dom=file_name[8:].split(".")[0]
         chan_id= file_name.split("_")[-1]

         # go through all channels from the given directory
         with open(f"{args.dir_channels}/channel-{dom[0:4]}.json", "r") as fr:
              data = json.load(fr)
              i = False
              for d in data["Channels"]:
                 for c in data["Channels"][d]:
                     if float(c["Id"]) == float(chan_id):
                         # if the channel is named, put it to relevant group
                         if c["Type"] not in c_types:
                             c_types[c["Type"]] = []
                         c_types[c["Type"]].append((dom, chan_id, file_name))
         counter += 1


print(c_types.keys())
for ct in c_types.keys():
    print(ct, len(c_types[ct]))
    ct2=""
    for i in ct:
        if i != " ":
            ct2 = ct2 + i
        else:
            ct2 = ct2 + "_"
    with open(f"list_cl_{args.name}_{ct2}.txt", "w") as fw:
        for c in c_types[ct]:
            fw.write(f"'{args.dir_gnuplot}/{c[2]}' lt rgb 'orange'\n")

sys.exit()

t = {}
 
for n2 in ch_list.keys():
    print(n2, end=" ")
print()

for n1 in ch_list.keys():
    t[n1]={}
    for n2 in ch_list.keys():
        dd = diff(n1, n2, ch_list)
        t[n1][n2] = dd
        print(int(dd), end=" ")
    print()

print("chlist", len(ch_list.keys()))

ch_all = list(ch_list.keys())
ch_rest = ch_all
clusters = []
while len(ch_rest) != 0:
    chosen = ch_rest[0]
    cluster = set_cluster(chosen, ch_rest, t)
    print("cluster", cluster)
    clusters.append(cluster)
print("len c", len(clusters))
print("clusters", clusters)

counter = 1
color_set = ["1,0,0", "0,1,0", "0,0,1", "0.5,0.5,0", "0.5,0,0.5", "0,0.5,0.5", "0.3,0.3,0.3", "0.5,0,0", "0,0.5,0", "0,0,0.5", "0.75,0.75,0", "0,0.75,0.75", "0.75,0,0.75", "0.5,0.75,0", "0,0.5,0.75", "0.5,0,0.75","0.75,0.5,0", "0,0.75,0.5", "0.75,0,0.5", "0.75,0.5,1", "1,0.75,0.5", "0.75,1,0.5","0.5,0.75,1", "1,0.5,0.75", "0.5,1,0.75","0.75,0.5,0"]
for cluster in clusters:
    counter += 1
    with open(f"list_cl_{args.name}_{counter}.txt", "w") as fw:
        for c in cluster:
            fw.write(f"'{args.dir}/{c}' lt rgb 'orange'\n")
    if len(color_set) > counter:
        color = color_set[counter]
    else:
        color = "1,1,1"
    print(f"lists/list_cl_{args.name}_{counter}.txt None {color}")


