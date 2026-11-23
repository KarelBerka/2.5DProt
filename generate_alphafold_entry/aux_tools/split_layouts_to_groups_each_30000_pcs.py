#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import json
import os
import shutil
from pathlib import Path


def is_in_tree(tdir, name):
    sname = name.split("-")
    tmp_file_name = f"{tdir}/{sname[1][-1]}/{sname[1][-2]}/{sname[1][-3]}/{name}"
    print("tmp_file_name", tmp_file_name, os.path.isfile(tmp_file_name))
    if os.path.isfile(tmp_file_name):
        return True
    else:
        return False


def add_file_to_dirs(res_dir, name, res2_dir, layout_dir):
    sname = name.split("-")
    tmp_file_name = f"{res_dir}/{sname[1][-1]}/{sname[1][-2]}/{sname[1][-3]}/{name}"
    if not os.path.isdir(f"{res_dir}/{sname[1][-1]}/{sname[1][-2]}/{sname[1][-3]}/"):
        os.makedirs(f"{res_dir}/{sname[1][-1]}/{sname[1][-2]}/{sname[1][-3]}/")
    print("---", f"{layout_dir}/{name}", f"{res2_dir}/{name}")
    shutil.copyfile(f"{layout_dir}/{name}", f"{res2_dir}/{name}")
#    print("pridano", tmp_file_name)
    Path(tmp_file_name).touch()

# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--layout_dir', action="store", required="True", help="input json")
parser.add_argument('--result_dir', action="store", required="True", help="input json")
parser.add_argument('--test_dir', action="store", required="True", help="input json")
parser.add_argument('--start', action="store", required="True", help="input json")
args = parser.parse_args()

layouts = 0
moved = 0
f_split = []
number=int(args.start)
res_dir = f"{args.result_dir}/res{number}"

if not os.path.isdir(f"{args.result_dir}"):
    os.mkdir(f"{args.result_dir}") 
if not os.path.isdir(res_dir):
    os.mkdir(res_dir)
f_split.append(f"res{number}")

counter = 0
for f in os.listdir(args.layout_dir):
     if f.endswith("json") and f.startswith("layout-"):
         layouts += 1
         if not is_in_tree(args.test_dir, f):
             print("is not in tree", f)
             add_file_to_dirs(args.test_dir, f, res_dir, args.layout_dir)
             counter += 1
             moved += 1
             if counter % 30000 == 0:
                  number = number + 1
                  res_dir = f"{args.result_dir}/res{number}"
                  os.mkdir(res_dir)
                  f_split.append(f"res{number}")
                  print("pridava se", number, counter)
              

print("counter", counter)
print("pridano", layouts)
print("do skupin", f_split)