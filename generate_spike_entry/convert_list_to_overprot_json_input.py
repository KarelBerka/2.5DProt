#!/usr/bin/env python3
import json
import argparse

# read arguments
parser = argparse.ArgumentParser(
    description='convert input list to overprot domain input list')
parser.add_argument('--input_list', required=True, action="store",
                    help="input list of p/c row should look like <1tqn A>")
parser.add_argument('--output_list', required=True, action="store",
                    help="output list name")
args = parser.parse_args()

data = []
with open(args.input_list,"r") as fd:
    line = fd.readline()
    while line:
        rec = {}
        s=line.split(" ")
        rec["domain"] = str(s[0])+str(s[1][:-1])+"00"
        rec["pdb"] = str(s[0])
        rec["chain_id"] = s[1][:-1]
        rec["ranges"] = ":"
        rec["auth_chain_id"] = rec["chain_id"]
        rec["auth_ranges"] = rec["ranges"]
        data.append(rec)
        line = fd.readline()


with (open(args.output_list, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ": "))




