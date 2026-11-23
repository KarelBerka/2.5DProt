#!/usr/bin/env python3
# tool which generate for given alphafold family
# its diagrams


import sys
import argparse
import time
import os
import subprocess
import shutil
import json

import utils_common
import svgwrite


def average(alist):
   asum = 0
   for i in alist:
       asum += i
   return asum / len(alist)


# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument('--all_directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
if args.directory[-1] =="/":
    name = args.directory[:-1]
else:
    name = args.directory

#new_dir = f"charges-{name}-{timestamp}"
#os.makedirs(new_dir, exist_ok=True)
#print("new_dir", new_dir)

all_data = {} # list of all sses which contains list of all charges for all res of this sses
fn_order = []
for filename in os.listdir(args.all_directory):
    print("filename", filename)
    with open(f"{args.all_directory}/{filename}", "r") as fd:
        print ("filename", filename)
        jdata = json.load(fd)
        fn_order.append(filename)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in all_data:
                    all_data[sse] = {}
                rmin = 1000#int(jdata["charges"][sse][0])
                for res in jdata["charges"][sse]:
                     if rmin > int(res):
                          rmin = int(res)
                for res in jdata["charges"][sse]:
                     n = int(res)-rmin
                     if n not in  all_data[sse]:
                         all_data[sse][n] = []
                     all_data[sse][n].append(jdata["charges"][sse][res])
#print("all_data", all_data['H2'])

data = {}
fn_order = []
for filename in os.listdir(args.directory):
    with open(f"{args.directory}/{filename}", "r") as fd:
        jdata = json.load(fd)
        fn_order.append(filename)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in data:
                    data[sse] = {}
                rmin = 1000
                for res in jdata["charges"][sse]:
                     if rmin > int(res):
                          rmin = int(res)
                for res in jdata["charges"][sse]:
                     n = int(res)-rmin
                     if n not in data[sse]:
                         data[sse][n] = []
                     data[sse][n].append(jdata["charges"][sse][res])

ssum = 0
for sse in data:
    for res in data[sse]:
        if res in all_data[sse]:
            ssum += abs(int(average(data[sse][res])- average(all_data[sse][res])))

text = ["1", "2", "3", "4a", "4b", "5", "6a", "6b", "7", "8", "9a", "9b", "10"]


for c in text:
    print(c, end = " ")
print()

print(f"{args.directory} {ssum}", end=" ")
#sys.exit()
sses_order = list(data.keys())
c = 0
for sse in sses_order:
    ssum = 0
    for res in data[sse]:
        if res in all_data[sse]:
            ssum += abs(int(average(data[sse][res])- average(all_data[sse][res])))
    print(f"{ssum}", end=" ")
    c += 1
print()

