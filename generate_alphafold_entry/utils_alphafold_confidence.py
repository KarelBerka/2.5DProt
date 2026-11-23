#!/usr/bin/env python3
import json


def read_res_confidence(file_name):
    res_conf = {}
    with (open(file_name, "r")) as fd:
        line = fd.readline()
        while line:
            if line.startswith("ATOM"):
                res_conf[int(line[22:27])] = float(line[61:66])
            line = fd.readline()
    return res_conf


def switch_htype(layout_fn, consensus_fn):
    # read confidence tag to set sses colors
    with (open(layout_fn, "r")) as fd:
        data = json.load(fd)
    # read consensus overprot file to set sses h_type
    with (open(consensus_fn, "r")) as fd:
        data_con = json.load(fd)
    sse_type = {}
    for hr in data["helices_residues"]:

        # set sses_type based on the sse occurrence - if it is wanted
        found = False
        for sse_rec in data_con["consensus"]["secondary_structure_elements"]:
            if sse_rec["label"] == hr:
                found = True
                if sse_rec["occurrence"] > 0.8:
                    sse_type[hr] = 1
                else:
                    sse_type[hr] = 2
        if not found:
            sse_type[hr] = 2

    data["h_type"] = sse_type

    with (open(layout_fn, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ":"))


def switch_htype_to_one(layout_fn):
    # read confidence tag to set sses colors
    with (open(layout_fn, "r")) as fd:
        data = json.load(fd)

    for hr in data["h_type"]:
        data["h_type"][hr] = 1

    with (open(layout_fn, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ":"))


def switch_color_to_confidence(layout_fn, confidence_fn, consensus_fn=None):
    # read confidence tag to set sses colors
    res = read_res_confidence(confidence_fn)
    with (open(layout_fn, "r")) as fd:
        data = json.load(fd)
    # read consensus overprot file to set sses h_type
    if consensus_fn is not None:
        with (open(consensus_fn, "r")) as fd:
            data_con = json.load(fd)
    color = {}
    sse_type = {}
    for hr in data["helices_residues"]:

        # we calculate the average sses res confidence
        minr, maxr = data["helices_residues"][hr]
        if minr > maxr:
            aux = maxr
            maxr = minr
            minr = aux
        i = minr
        conf = 0

        while i <= maxr:
            conf += res[i]
            i += 1

        if maxr-minr+1 == 0:
            minconf = minr
        else:
            minconf = conf / (maxr-minr+1)

        # based of average sses confidence we use the same color as alfa fold
        if minconf > 90:
            color[hr] = "#0053d6"
        elif minconf > 70:
            color[hr] = "#65cbf3"
        elif minconf > 50:
            color[hr] = "#FFdb13"
        else:
            color[hr] = "#FF7d45"

        # set sses_type based on the sse occurrence - if it is wanted
        if consensus_fn is not None:
            found = False
            for sse_rec in data_con["consensus"]["secondary_structure_elements"]:
                if sse_rec["label"] == hr:
                    found = True
                    if sse_rec["occurrence"] > 0.8:
                        sse_type[hr] = 1
                    else:
                        sse_type[hr] = 2
            if not found:
                sse_type[hr] = 2
        else:
            sse_type[hr] = 1

    data["color"] = color
    data["h_type"] = sse_type

    with (open(layout_fn, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ":"))

