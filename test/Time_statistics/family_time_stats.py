#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import os
import os.path
import logging

# Load params - input family generation dir
if len(sys.argv) < 2:
    print("ERROR: argument missing")
    print("Usage:")
    print("    ./family_time_stats.py family_generated_dbs_dir")
    sys.exit()
else:
    input_dir = sys.argv[1]

stat = {}
label = {}
single_images = 0
counter = 0

for dirname in os.listdir(input_dir):
    if dirname.startswith("generated-"):
       for filename in os.listdir(f"{input_dir}/{dirname}"):
          fn = f"{input_dir}/{dirname}/{filename}"
          if os.path.isfile(fn) and filename.endswith("l.svg"):
#          if os.path.isfile(fn) and filename.startswith("image-") and filename.endswith(".svg"):
             single_images += 1
          if os.path.isfile(fn) and filename.startswith("main"):
             with open(fn, "r") as fd:
                l={}
                line = fd.readline()
                counter += 1
                while line:
                    if "Timestamp" in line:
                        val = line[:-5].rfind(" ")
                        c = int(line[15])
                        if c not in stat:
                            stat[c] = float(line[val:-5])
                        else:
                            stat[c] += float(line[val:-5])
                        l[c] = float(line[val:-5])
                        label[c] = line[17:val - 20]
                    line = fd.readline()
                if (2 in l) and (6 in l):
                    if (l[2] * 2 > l[6]) and l[2] > 10:
                    # annotation takes more then 1/2 of all time - probably annotation problem
                        print(f"{fn} probably annotation problem: {int(l[2])} - {int(l[2]*100/l[6])}procent z {int(l[6])}")

if stat == {}:
    print("There is no family generated directory")
    sys.exit()

print(f"read records about {counter} families")
for i in stat:
    print(f"part {i} takes {int(stat[i])} minutes/ {int(stat[i])/60/24} days  (in procent {int(stat[i]*100/stat[6])} ) - activity description '{label[i]}'")

ph = single_images * 60 / stat[6]
pd = ph * 24
ndomain = 471275
days_dbs = ndomain / pd

print(f"generates single images {single_images} that means {ph} per hour {pd} per day, {int(days_dbs)} days for whole dbs - {ndomain} domains")

