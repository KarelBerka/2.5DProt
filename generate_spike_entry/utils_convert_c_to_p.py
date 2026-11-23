#!/usr/bin/env python3
import json
import sys
import os

def rename_annotation_chain(chain_n, pref , filename, dfilename, outputfn):
    with open(filename, "r") as fd:
        data = json.load(fd)
        dom = list(data.keys())[0]
        with open(dfilename, "r") as fdd:
            ddata = json.load(fdd)
            ddom = list(ddata.keys())[0]
            for i in data[dom]["secondary_structure_elements"]:
                for j in ddata[ddom]["secondary_structure_elements"]:
                      if i["start"] == j["start"] and i["end"] == j["end"] and i["chain_id"] == j["chain_id"]:
                          i["start_vector"] = j["start_vector"]
                          i["end_vector"] = j["end_vector"]
                i["chain"] = chain_n
                i["chain_id"] = chain_n
                i["label"] = i["label"][0]+pref+i["label"][1:]
            for i in data[dom]["beta_connectivity"]:
                i[0] = i[0][0]+pref+i[0][1:]
                i[1] = i[1][0]+pref+i[1][1:]

    if os.path.isfile(outputfn):
        # add chain annotation file to existing one
        with open(outputfn, "r") as fd:
            data2 = json.load(fd)
            dom2 = list(data2.keys())[0]
        for i in data[dom]["secondary_structure_elements"]:
            data2[dom2]["secondary_structure_elements"].append(i)
        for i in data[dom]["beta_connectivity"]:
            data2[dom2]["beta_connectivity"].append(i)
        with (open(outputfn, "w")) as fd:
            json.dump(data2, fd, indent=4, separators=(",", ": "))
    else:
        # copy chain annotation file to protein one
        with (open(outputfn, "w")) as fd:
            json.dump(data, fd, indent=4, separators=(",", ": "))

def rename_diagram_stuff(pref, filename, outputfn, pname):
    with open(filename, "r") as fd:
        data = json.load(fd)
        for i in data["nodes"]:
            i["label"] = i["label"][0]+pref+i["label"][1:]

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
            json.dump(data, fd, indent=4, separators=(",", ": "))


def convert_annotation_c_to_p(d, number, dirzb, pname, outputdir):

    #create output_dir
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    # aggregate t1 t2 t3 to protein annotation file
    output_name = f"{outputdir}/{pname}-annotation.sses.json"
    for fn in os.listdir(outputdir):
        if fn.startswith(pname) and fn.endswith("-detected.sses.json"):
            d1 = f"{outputdir}/{fn}"

    i = 1
    while i <= number:
        t = ""
        for fn in os.listdir(d[i-1]):
            if fn.startswith(pname) and fn.endswith("annotated.sses.json"):
                t = f"{d[i-1]}/{fn}"
        if t != "":
            rename_annotation_chain(chr(ord("A")+i-1), f"{i}", t, d1, output_name)
        else:
            print(f"{pname} has no annotation in chain {t}")

        i += 1

    t4 = ""
    if dirzb != None:
        for fn in os.listdir(dirzb):
            if fn.startswith(pname) and fn.endswith("annotated.sses.json"):
                t4 = f"{dirzb}/{fn}"
                num = 70 + ord(t4[-23]) - ord("A")
                nch = chr(ord(t4[-23]) - ord("A")+7)
                rename_annotation_chain(nch, str(num), t4, d1, output_name)


def convert_diagram_c_to_p(d, number, pname, outputdir):

    #create output_dir
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    # aggregate t1 t2 t3 to protein annotation file
    output_name = f"{outputdir}/diagram.json"
    i = 1
    while i <= number:
        rename_diagram_stuff(f"{i}", f"{d[1-1]}/diagram.json", output_name, pname)
        i += 1
