#!/usr/bin/env python3
import argparse
import json


def read_file_list(filename):
    fl = []
    with open(filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.split(" ")
            fl.append([sline[0], sline[1], sline[2][:-1]])
            line = fd.readline()
    return fl


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


def generate_stats(filename, directory):
    domain_set = []
    radius_set = []
    channel_counter = 0
    print("gmi", filename)
    with open(filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.strip().split("'")
            ssline = sline[1].split("/")
            fname = ssline[-1][8:]
            d_name = fname.split(".")[0]
            c_name = fname.split("_")[-1]
            channel_file_name = f"{directory}/{d_name}-channel.json"
            channel_counter += 1
#            print("sline", channel_file_name, c_name)
            with open(channel_file_name, "r") as chfn:
                data = json.load(chfn)
#               print("chfn", data["Channels"])
                for ch_type in data["Channels"]:
                    for channel in data["Channels"][ch_type]:
                        if str(c_name) == str(channel["Id"]):
                            print("*", end="")
                            if d_name not in domain_set:
                                domain_set.append(d_name)
                            r_set = []
                            for p in channel["Profile"]:
                                r_set.append(p["Radius"])
                            radius_set.append(min(r_set))
#                             print("JOOO")
            line = fd.readline()
    print()
    print("-------------")
    print(f"file name: {filename}, channels: {channel_counter}, in domains: {len(domain_set)}, avg_radius: {avg(radius_set)}, max_radius: {max(radius_set)}")
    print(f"<tr><td>{filename}</td><td>{channel_counter}</td><td>{len(domain_set)}</td><td>{avg(radius_set)}</td><td>{max(radius_set)}</td></tr>")
    domain_set.sort()
#    print("domain_set", domain_set)


# read arguments
parser = argparse.ArgumentParser(
    description='show input channels and sses files using mayavi tool')
parser.add_argument('--directory', required=True, action="store",
                    help="channel data directory - usually generated-<family>-<time>/annotation dir")
parser.add_argument('--channel_list', required=True, action="store",
                    help="set of channels in the given choice (in format TODO)")
args = parser.parse_args()


generate_stats(args.channel_list, args.directory)
