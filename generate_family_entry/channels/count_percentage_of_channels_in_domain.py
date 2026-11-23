#!/usr/bin/env python3
# input - layout list in main dir
#       - chanels list (usually in annotation dir if the famili is generated with --debug
# go through layout list, count channels sizes and compare with the hole size of channels from orig files

import argparse
import json
import os
import math
import sys


def measure_channel_length(ch_name, p_name, ch_dir):
    ch_file_name = f"./{ch_dir}/channel-{p_name}.json"
    with open(ch_file_name, "r") as fc:
        data = json.load(fc)
        for ct in data['Channels']:
            for c in data['Channels'][ct]:
                if float(c['Id']) == float(ch_name):
                    start = None
                    sum_len = 0
                    for ch in c['Profile']:
                        if start is None:
                            start = [ch['X'], ch['Y'], ch['Z']]
                        else:
                            end = [ch['X'], ch['Y'], ch['Z']]
                            ch_len = math.sqrt((start[0]-end[0])*(start[0]-end[0])+(start[1]-end[1])*(start[1]-end[1])+(start[2]-end[2])*(start[2]-end[2]))
                            sum_len += ch_len
                            start = end

                    return sum_len


def measure_layout_channel(ch_list, data, p_name, ch_dir):
    # measure channel part length from layout file
    channel_old_name = ""
    channel_old_coords = [0, 0, 0]
    channel_len = 0
    len_sum = 0
    for i in data:
        c = 0
        mid = 0
        while c < len(i):
            if i[c].isalpha():
                mid = c
                break
            c += 1
        channel_name = i[:mid]
        channel_coords = data[i]['3d_layout']
        if (channel_name == channel_old_name):
            channel_len = math.sqrt(((channel_coords[0]-channel_old_coords[0])*(channel_coords[0]-channel_old_coords[0]))+((channel_coords[1]-channel_old_coords[1])*(channel_coords[1]-channel_old_coords[1]))+((channel_coords[2]-channel_old_coords[2])*(channel_coords[2]-channel_old_coords[2])))
            channel_old_coords = channel_coords
            len_sum += channel_len
        elif channel_old_name != "":
            ch_len_sum = measure_channel_length(channel_old_name, p_name, ch_dir)
            if ch_len_sum != 0 and ch_len_sum >= len_sum:
                ch_list.append(int(len_sum*100/ch_len_sum))
            channel_old_name = channel_name
            channel_old_coords = channel_coords
            len_sum = 0
        else:
            channel_old_name = channel_name
            channel_old_coords = channel_coords
            len_sum = 0
    return ch_list


# read arguments
parser = argparse.ArgumentParser(
    description='TODO')
parser.add_argument('--layout_dir', required=True, action="store",
                    help="directory containing layout files")
parser.add_argument('--channels_dir', required=True, action="store",
                    help="directory containing channels files")
parser.add_argument('--verbose', required=False, action="store",
                    help="do more verbose output")

args = parser.parse_args()

ch_list = []

for fn in os.listdir(args.layout_dir):
    if fn.startswith("layout") and fn.endswith(".json"):
        with open(f"{args.layout_dir}/{fn}", "r") as fdl:
            data = json.load(fdl)
            if len(data["channels"]) > 0:
                ch_list = measure_layout_channel(ch_list, data["channels"], fn[7:11], args.channels_dir)

c = 0
s = 0
aux = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}
for i in ch_list:
    c += 1
    s += i
    aux[int(i / 10)] += 1
if c == 0:
    sys.exit()

print(f"{args.layout_dir}", end=": ")
for i in aux:
    print(f"{int(aux[i]* 100/c)}", end=" ")
print(" || ", c, s/c)
