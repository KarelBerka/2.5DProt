#!/usr/bin/env python3

import json

def transform_layout_annotation(orig_sl, tr_l, start_layout):
    # transform sse names of orig_sl to names present in tr_l
    # and save orig_sl layout with transformed names to start_layout

    print("orig_sl", orig_sl)

    # open orig_sl and read its data
    with open(orig_sl, "r") as fo:
       data_orig = json.load(fo)

    print(data_orig.keys())

    print("tr_l", tr_l)
    # open tr_l and read its data
    with open(tr_l, "r") as ft:
       data_trans = json.load(ft)


    print(data_trans.keys())
    # find transformation function based on res numbers
    data = {}
    data["helices_residues"]={}
    data["layout"]={}
    data["angles"]={}
    data["color"]={}
    for sseo in data_orig["helices_residues"]:
#        print(sseo)
        found = False
        for sset in data_trans["helices_residues"]:
            if data_orig["helices_residues"][sseo] == data_trans["helices_residues"][sset]:
#                print("stejne", sseo, sset, data_orig["helices_residues"][sseo])
                data["helices_residues"][sset]=data_orig["helices_residues"][sseo]
                data["layout"][sset]=data_orig["layout"][sseo]
                data["angles"][sset]=data_orig["angles"][sseo]
                data["color"][sset]=data_orig["color"][sseo]


#        if not found:
#           print("nenaslo se")


    # save the result start layout
    with open(start_layout, "w") as fsl:
        json.dump(data, fsl, indent=4, separators=(",", ":"))
