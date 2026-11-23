#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import os
import os.path
import logging

# Load params - input family generation dir
if len(sys.argv) < 2:
    print("ERROR: Argument missing")
    print("Usage:")
    print("    ./protein_time_stats protein enerated dbs dir")
    sys.exit()

input_dir = sys.argv[1]

stat = {}
label = {}
single_images = 0

for dirname1 in os.listdir(input_dir):
 for dirname2 in os.listdir(f"{input_dir}/{dirname1}"):
  for dirname3 in os.listdir(f"{input_dir}/{dirname1}/{dirname2}"):
#    print("dirname3", dirname3)
    if dirname3.startswith("generated-"):
       for filename in os.listdir(f"{input_dir}/{dirname1}/{dirname2}/{dirname3}"):
          fn = f"{input_dir}/{dirname1}/{dirname2}/{dirname3}/{filename}"
#          if (os.path.isfile(fn) and filename.startswith("image-") and filename.endswith(".svg")):

             
          if os.path.isfile(fn) and filename.startswith("main"):
             single_images += 1
             with open(fn, "r") as fd:
                line = fd.readline()
#                print("fn", fn)
                while line:
                    if "Timestamp" in line:
                       val = line[:-5].rfind(" ")
                       c=int(line[15])
                       if c != 5:
#                           print(float(line[val:-5]))
                           if c not in stat:
                               stat[c] = float(line[val:-5])
                           else:
                               stat[c] += float(line[val:-5])
                           label[c] = line[15:val]
                       else:
                           if (float(line[val:-5]) > 50):
#                                print(float(line[val:-5]), "vetsi filename:", filename, "c:", c, "line[15:val]", line[15:val], "float", float(line[val:-5]), "line", line)
                               if c not in stat:
                                   stat[c] = float(line[val:-5])
                               else:
                                   stat[c] += float(line[val:-5])
                               label[c] = line[15:val]
#                                k = 1
                           else:
#                               print("mensi", filename, c, line[15:val], float(line[val:-5]))
                               if c not in stat:
                                   stat[c] = float(line[val:-5])
                               else:
                                   stat[c] += float(line[val:-5])
                               label[c] = line[15:val]

                    line = fd.readline()
#       print("hh", stat)
#sys.exit()
if stat == {}:
    print("There is no protein generated directory")
    sys.exit()

for i in stat:
    print(f"part: {i}, time perc: {int(stat[i]*100/stat[8])},  label: {label[i]} time(min): {int(stat[i])}, time (hour): {int(stat[i])/60}")


j=190000
ph=single_images*60/stat[8]
pd=ph*24
days_dbs=j/pd

print(f"generates single images {single_images} that means {ph} per hour {pd} per day, {int(days_dbs)} days for whole dbs - {j} domains")

