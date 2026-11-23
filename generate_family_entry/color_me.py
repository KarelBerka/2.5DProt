#!/usr/bin/env python3
import json
import math


def switch_color_to_stats(layout_fn, stats_fn):
    with (open(layout_fn, "r")) as fd:
        data = json.load(fd)
    with (open(stats_fn, "r")) as fd:
        data_stats = json.load(fd)
    color = {}
    for hr in data["helices_residues"]:
       # we calculate the sses confidence as the minimum of confidence of all sses res
       # without the two on the ends
        if hr in data_stats:
            minconf = data_stats[hr]["probability"]*100
        else:
            minconf = 10
        print("hr", hr, minconf)

        # based of average sses confidence we use the same color as alfa fold
        if minconf > 90:
            color[hr] = "#0053d6"
        elif minconf > 70:
            color[hr] = "#65cbf3"
        elif minconf > 50:
            color[hr] = "#FFdb13"
        else:
            color[hr] = "#FF7d45"
        print("color[hr]", color[hr])

    data["color"] = color

    with (open(layout_fn, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ":"))


def switch_color_to_variance(layout_fn, output_fn, stats_fn):
    with (open(layout_fn, "r")) as fd:
        data = json.load(fd)
    with (open(stats_fn, "r")) as fd:
        data_stats = json.load(fd)
    color = {}
    for hr in data["helices_residues"]:
        # we calculate the sses confidence as the minimum of confidence of all sses res
        # without the two on the ends
#        print("dd", data_stats["consensus"]["secondary_structure_elements"])
        minconf = 10
        for sses in data_stats["consensus"]["secondary_structure_elements"]:
            if sses['label'] == hr:
                minconf = (15 - math.sqrt(sses["variance"]))*7
                break

        # based of average sses confidence we use the same color as alfa fold
        if minconf > 90:
            color[hr] = "#0053d6"
        elif minconf > 70:
            color[hr] = "#65cbf3"
        elif minconf > 50:
            color[hr] = "#FFdb13"
        else:
            color[hr] = "#FF7d45"

    data["color"] = color

    with (open(output_fn, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ":"))


