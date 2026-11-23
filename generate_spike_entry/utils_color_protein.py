#!/usr/bin/env python3
import json

def find_all_chains(detect_fn):
    with open(detect_fn,"r") as fd:
        data = json.load(fd)
    i = data.keys()
    chains = set()
    for j in data[list(i)[0]]["secondary_structure_elements"]:
        chains.add(j['chain_id'])
    return chains


def switch_color_to_chain(number, detect_fn):
    chains=find_all_chains(detect_fn)
    n = 1
    cnum = {}
    for i in chains:
        cnum[i] = 900.0/len(chains)*n
        n = n+1
    with (open(detect_fn, "r")) as fd:
        data2 = json.load(fd)
    i = data2.keys()
    color = {}
    for sses_rec in data2[list(i)[0]]["secondary_structure_elements"]:
            chain = sses_rec["chain_id"]
            if number == 3:
                if chain == "A": 
                    sses_rec["rainbow"] = "s900"
                elif chain == "B":
                    sses_rec["rainbow"] = "s600"
                elif chain == "C":
                    sses_rec["rainbow"] = "s300"
                else:
                    sses_rec["rainbow"] = "s750"
            if number == 6:
                if chain == "A": 
                    sses_rec["rainbow"] = "s900"
                elif chain == "B":
                    sses_rec["rainbow"] = "s600"
                elif chain == "C":
                    sses_rec["rainbow"] = "s300"
                elif chain == "D": 
                    sses_rec["rainbow"] = "s750"
                elif chain == "E":
                    sses_rec["rainbow"] = "s450"
                elif chain == "F":
                    sses_rec["rainbow"] = "s150"
                else:
                    sses_rec["rainbow"] = "s000"
#           sses_rec["rainbow"] = str("s"+str(int(cnum[chain])))
    with (open(detect_fn, "w")) as fd:
        json.dump(data2, fd, indent=4, separators=(",", ":"))
