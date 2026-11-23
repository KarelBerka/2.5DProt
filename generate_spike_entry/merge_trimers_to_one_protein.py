#!/usr/bin/env python3
import json
import argparse
import sys
import os

def rename_annotation_chain(chain_n, pref , filename, outputdir, pname):
    with open(filename, "r") as fd:
        data = json.load(fd)
        dom = list(data.keys())[0]
        for i in data[dom]["secondary_structure_elements"]:
            i["chain"] = chain_n
            i["chain_id"] = chain_n
            i["label"] = i["label"][0]+pref+i["label"][1:]

    outputfn = f"{outputdir}/{pname}annotation.sses.json"
    if os.path.isfile(outputfn):
        # add chain annotation file to existing one
        with open(outputfn, "r") as fd:
            data2 = json.load(fd)
            dom2 = list(data2.keys())[0]
        for i in data[dom]["secondary_structure_elements"]:
            data2[dom2]["secondary_structure_elements"].append(i)
        with (open(outputfn, "w")) as fd:
            json.dump(data2, fd, indent=4, separators=(",", ": "))
    else:
        # copy chain annotation file to protein one
        with (open(outputfn, "w")) as fd:
            json.dump(data, fd, indent=4, separators=(",", ": ") )

def rename_diagram_stuff(pref, filename, outputdir, pname):
    with open(filename, "r") as fd:
        data = json.load(fd)
        for i in data["nodes"]:
            i["label"] = i["label"][0]+pref+i["label"][1:]

    outputfn = f"{outputdir}/{pname}diagram.json"
    if os.path.isfile(outputfn):
        # add chain annotation file to existing one
        with open(outputfn, "r") as fd:
            data2 = json.load(fd)
        for i in data["nodes"]:
            data2["nodes"].append(i)
        with (open(outputfn, "w")) as fd:
            json.dump(data2, fd, indent=4, separators=(",", ": "))
    else:
        # copy chain annotation file to protein one
        with (open(outputfn, "w")) as fd:
            json.dump(data, fd, indent=4, separators=(",", ": ") )
   


# read arguments
parser = argparse.ArgumentParser(
    description='convert input list to overprot domain input list')
parser.add_argument('--t1', required=True, action="store", help="trimer 1")
parser.add_argument('--t2', required=True, action="store", help="trimer 2")
parser.add_argument('--t3', required=True, action="store", help="trimer 3")
parser.add_argument('--d1', required=True, action="store", help="trimer 1")
parser.add_argument('--d2', required=True, action="store", help="trimer 2")
parser.add_argument('--d3', required=True, action="store", help="trimer 3")
parser.add_argument('--output_dir', required=True, action="store",
                    help="output dir name")
args = parser.parse_args()

# find protein name
s = args.t1.split('/')
fn = s[-1]
c = 0
ind = -1
while c < len(fn):
    if fn[c].isupper():
        ind=c
        break
    c += 1 

if ind != -1:
    name=fn[:ind]+"-"
else:
    name = ""

#create output_dir
if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)

# aggregate t1 t2 t3 to protein annotation file
rename_annotation_chain("A", "1", args.t1, args.output_dir, name)
rename_diagram_stuff("1", args.d1, args.output_dir, name)
rename_annotation_chain("B", "2", args.t2, args.output_dir, name)
rename_diagram_stuff("2", args.d2, args.output_dir, name)
rename_annotation_chain("C", "3", args.t3, args.output_dir, name)
rename_diagram_stuff("3", args.d3, args.output_dir, name)

