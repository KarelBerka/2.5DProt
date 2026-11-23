#!/usr/bin/env python3

import argparse
import os
import json
import sys
import math


def test_unknown_channel_and_known_channel(unknownch, ch):
    mindiff = 0
    for sses in unknownch["Profile"]:
        ssmin = 200
        for sses2 in ch["Profile"]:
            diff = math.sqrt((sses["X"]-sses2["X"])*(sses["X"]-sses2["X"])+(sses["Y"]-sses2["Y"])*(sses["Y"]-sses2["Y"])+(sses["X"]-sses2["X"])*(sses["X"]-sses2["X"]))
            if diff < ssmin:
                ssmin = diff
        if mindiff < ssmin:
            mindiff = ssmin
    return mindiff


# test unknown channel against all annotated chasnnels
def test_unknown_channel(unknownch, data):
    mind = 1000
    ind = None
    if "Annotations" in data:
        for an in data["Annotations"]:
            if "Channels" in data:
                for cc1 in data["Channels"]:
                    for ch in data["Channels"][cc1]:
                        if ch["Id"] == an["Id"]:
                            diff = test_unknown_channel_and_known_channel(unknownch, ch)
                            if diff < mind:
                                mind = diff
                                ind = {an["Name"], an["Id"]}
    return mind, ind


# test unknown channel against all unknown channels
def test_unknown_channel2(unknownch, data):
    mind = 1000
    ind = None
    if "Channels" in data:
        for cc1 in data["Channels"]:
            for ch in data["Channels"][cc1]:
                if ch["Id"] != unknownch["Id"]:
                    diff = test_unknown_channel_and_known_channel(unknownch, ch)
                    if diff < mind:
                        mind = diff
                        ind = {ch["Id"], "---"}
    return mind, ind


# start
def test_channel_data_file(channel_dir, channel_file):
    with open(f"{channel_dir}/{channel_file}", "r") as fd:
        print(f" fn {channel_dir}/{channel_file}")
        data = json.load(fd)
        if "Channels" in data:
            for cc1 in data["Channels"]:
                for ch in data["Channels"][cc1]:
                    annotated = False
                    if "Annotations" in data:
                        for an in data["Annotations"]:
                            if ch["Id"] == an["Id"]:
                                print(f"1Annotations {ch['Id']}", an["Id"], an["Name"])
                                annotated = True
                    if not annotated:
                        # have to be tested
                        mind, ind = test_unknown_channel(ch, data)
                        if mind < 1:
                            print(f"{channel_file} diff {ch['Id']} {ind} {mind} JO")
                        else:
                            print(f"{channel_file} diff {ch['Id']} {ind} {mind} NE")
                        mind, ind = test_unknown_channel2(ch, data)
                        if mind < 1:
                            print(f"{channel_file} diff {ch['Id']} {ind} {mind} JO")
                        else:
                            print(f"{channel_file} diff {ch['Id']} {ind} {mind} NE")


# read arguments
parser = argparse.ArgumentParser(
    description='generate statistic regarding channel with channel names from the input directory')
parser.add_argument('--directory', required=True, action="store",
                    help="channel files directory")

args = parser.parse_args()

channels_sses_min_dist = {}
unknown_sses_min_dist = {}
print("00 - nacitam tunely a jmenovane tunely")
for fn in os.listdir(args.directory):
    if fn.endswith(".json") and fn.startswith("channel"):
        # it is a channel file
        # zjisti jestli v unknown neni nejakej known
        test_channel_data_file(args.directory, fn)
sys.exit()
