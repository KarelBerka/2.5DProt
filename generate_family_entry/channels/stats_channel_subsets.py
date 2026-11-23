#!/usr/bin/env python3
import argparse
import os


def generate_axis_set(ch_filename):
    axis = {'x': [],
            'y': [],
            'z': [],
            }
    with open(ch_filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.split(" ")
            axis['x'].append(float(sline[0]))
            axis['y'].append(float(sline[1]))
            axis['z'].append(float(sline[2]))
            line = fd.readline()
    return axis


def avg(data_list):
    return sum(data_list)/len(data_list)


def generate_stats(filename):
    domain_set = []
    radius_set = []
    channel_counter = 0
    print("gmi", filename)
    with open(filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.strip().split("'")
            ssline = sline[1].split("/")
            sssline = ssline[2].split("_")
#           # only MOLE channels
            if "." not in sssline[1]:
                # only CAVER channels
                fname = ssline[-1][8:]
                d_name = fname.split(".")[0]
                if d_name not in domain_set:
                    domain_set.append(d_name)
            line = fd.readline()
    print()
    print("-------------")
    print(f"file name: {filename}, channels: {channel_counter}, in domains: {len(domain_set)}, avg_radius: {avg(radius_set)}, max_radius: {max(radius_set)}")
    print(f"<tr><td>{filename}</td><td>{channel_counter}</td><td>{len(domain_set)}</td><td>{avg(radius_set)}</td><td>{max(radius_set)}</td></tr>")
    return domain_set


def intersection(list1, list2):
    c = 0
    for i in list1:
        if i in list2:
            c += 1
    return c


def intersection2(list1, list2):
    c = []
    for i in list1:
        if i in list2:
            c.append(i)
    return c


# read arguments
parser = argparse.ArgumentParser(
    description='show input channels and sses files using mayavi tool')
parser.add_argument('--directory', required=True, action="store",
                    help="channel data directory - usually generated-<family>-<time>/annotation dir")
parser.add_argument('--channel_list_dir', required=True, action="store",
                    help="directory containing set of channels in the given choice (in format TODO)")
args = parser.parse_args()

channels = {}

for fn in os.listdir(args.channel_list_dir):
    channels[fn] = generate_stats(f"{args.channel_list_dir}/{fn}")

domains = {}
for chan in channels:
    for d_name in channels[chan]:
        if d_name in domains:
            domains[d_name].append(chan)
        else:
            domains[d_name] = [chan]

print("                  ", end=" || ")
for ch1 in channels:
    print(f"{ch1:12}", end="|")
print()

for ch1 in channels:
    print(f"{ch1:18}", end=" || ")
#    print("ch1", ch1, channels[ch1])
    for ch2 in channels:
        print(f"{len(channels[ch2]):3d}({intersection(channels[ch2], channels[ch1]):3d}){len(channels[ch1]):3d}", end=" | ")
    print()


print("procenta")

print("                  ", end=" || ")
for ch1 in channels:
    print(f"{ch1}", end="|")
print()

for ch1 in channels:
    print(f"{ch1}({len(channels[ch1])})", end=" || ")
#    print("ch1", ch1, channels[ch1])
    for ch2 in channels:
        a = len(channels[ch2])
        b = intersection(channels[ch2], channels[ch1])
        c = len(channels[ch1])
        r = a*100/len(domains)
#        print(f"{int(b*100/c):3d} ({int(r)})%", end= " | ")
#        print(f"{int(b*100/c/r*100)-100:3d}", end= " | ")
        print(f"{int(b*100/c/r*100)-100}", end=" ")
    print()

print("ld", len(domains))
