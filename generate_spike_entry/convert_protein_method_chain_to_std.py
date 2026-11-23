#!/usr/bin/env python3
import json
import sys
import os
import copy

def rename_layout_sses(filename):
    with open(filename, "r") as fd:
        data = json.load(fd)
        dom = list(data.keys())[0]
        new_data = {}
        for label in data.keys():
            if (label == "layout") or (label == "size") or (label == "helices_residues") or (label == "angles") or (label == "h_type") or (label == "chain_id") or (label == "color") or (label == "domain") or (label == "3dlayout"):
                new_data[label] = {}
                for sses in data[label]:
                    if label == "color":
                        if sses[1]=="0":
                           # this is the first trimer, we have to remove the start "0" symbol
                           new_sses = sses[0]+sses[2:]
                           if sses[2] == "1":
                               new_data[label][new_sses] = "s900"
                           if sses[2] == "2":
                               new_data[label][new_sses] = "s600"
                           if sses[2] == "3":
                               new_data[label][new_sses] = "s300"
                        else:
                            # this is the rest of trimer - we have to increase the start symbol +3 - to have space for trimer
                            new_sses = sses[0]+str(int(sses[1]) + 3) + sses[2:]
                            new_data[label][new_sses] = "grey"
                    elif label == "chain_id":
                        if sses[1]=="0":
                           # this is the first trimer, we have to remove the start "0" symbol
                           new_sses = sses[0]+sses[2:]
                           if sses[2] == "1":
                               new_data[label][new_sses] = "A"
                           if sses[2] == "2":
                               new_data[label][new_sses] = "B"
                           if sses[2] == "3":
                               new_data[label][new_sses] = "C"
                        else:
                            # this is the rest of trimer - we have to increase the start symbol +3 - to have space for trimer
                            new_sses = sses[0]+str(int(sses[1]) + 3) + sses[2:]
                            new_data[label][new_sses] = chr(ord(data[label][sses])+4)
                    else:
                        if sses[1]=="0":
                            # this is the first trimer, we have to remove the start "0" symbol
                            new_sses = sses[0]+sses[2:]
                            new_data[label][sses[0]+sses[2:]] = data[label][sses]
                        else:
                            # this is the rest of trimer - we have to increase the start symbol +3 - to have space for trimer
                            new_sses = sses[0]+str(int(sses[1]) + 3) + sses[2:]
                            new_data[label][new_sses] = data[label][sses]
            if (label == "generated") or (label == "protein"):
                new_data[label] = data[label]
            if (label == "beta_connectivity"):
                new_data[label] = []
                for cl in data[label]:
                    new = ""
                    for sses in cl.split(":"):
                       if sses[0]=="0":
                           # this is the first trimer, we have to remove the start "0" symbol
                           new_sses = sses[0]+sses[2:]
                       else:
                           # this is the rest of trimer - we have to increase the start symbol +3 - to have space for trimer
                           new_sses = sses[0]+str(int(sses[1]) + 3) + sses[2:]
                       if new == "":
                           new = f"{new_sses}"
                       else:
                           new +=  f":{new_sses}"
                    new_data[label].append(new)
    with open(filename, "w") as fd:
        json.dump(new_data, fd, indent=4, separators=(",", ":"))


#rename_layout_sses("layout-7a98.json")
