#!/usr/bin/env python3

import argparse
import json
import os
import svgwrite


blue = "#3274A1"
orange = "#E1812C"


def add_records(input_data, organism_name):
    # add a records to data_list structure
    # from input data structure
    # column is column name and added data
    # organism_name is the name of organism which will be displayed on diagram
    # if the s[lit plot is wanted use vg_split and vg_split_column
    if len(organism_name) > 5:
        org_name = f"{organism_name[:5]}."
    else:
        if len(organism_name) == 4:
            org_name = f"{organism_name[0][:2]}_{organism_name[1][:2]}"
        else:
            org_name = organism_name[0][:4]
    spliton = organism_name[0].split(" ")
    if len(spliton) == 1:
        org_name = f"{spliton[0][:4]}. "
    elif len(spliton[1]) < 8:
        org_name = f"{spliton[0][:4]}. {spliton[1]}"
    else:
        org_name = f"{spliton[0][:4]}. {spliton[1][:7]}."
    # now the organism name is set
    org_data = {"helices_cumulative": {},
                "strands_cumulative": {},
                "org_name_shortcut": org_name,
                "org_name": organism_name[0],
                "helices_sum": 0,
                "strands_sum": 0,
                "helices_fivetyplus": 0,
                "strands_fivetyplus": 0,
                }
    for sses_num in data["strands"]:
        org_data["strands_cumulative"][int(sses_num)] = data["strands"][sses_num]
        if int(sses_num) >= 50:
            org_data["strands_fivetyplus"] += data["strands"][sses_num]
        org_data["strands_sum"] += data["strands"][sses_num]
    for sses_num in data["helices"]:
        org_data["helices_cumulative"][int(sses_num)] = data["helices"][sses_num]
        if int(sses_num) >= 50:
            org_data["helices_fivetyplus"] += data["helices"][sses_num]
        org_data["helices_sum"] += data["helices"][sses_num]
    return org_data


def violin2(tr, k, org, img, mm):
    m = 50
#    print("org", org)
    c = 0
    while c < m:
        if c in org["strands_cumulative"]:
            dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr-int(org["strands_cumulative"][c]/org["strands_sum"]*100/2*k)-1, mm-c*3+2), stroke_width=4, stroke="black"))
        else:
            dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr-1, mm-c*3+2), stroke_width=4, stroke="black"))
        c += 1
    dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr-int(org["strands_fivetyplus"]/org["strands_sum"]*100/2*k)-1-1, mm-c*3+2), stroke_width=4, stroke="black"))
#    print("ooo", tr+1, tr-int(org["strands_fivetyplus"]/org["strands_sum"]*100/2*k)-1-1, org["strands_fivetyplus"], org["strands_sum"])
    c = 0
    while c < m:
        if c in org["helices_cumulative"]:
            dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr+1+int(org["helices_cumulative"][c]/org["helices_sum"]*100/2*k)*2+1, mm-c*3+2), stroke_width=4, stroke="black"))
        else:
            dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr+2, mm-c*3+2), stroke_width=4, stroke="black"))
        c += 1
    dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr+1+int(org["helices_fivetyplus"]/org["helices_sum"]*100/2*k)*2+1, mm-c*3+2), stroke_width=4, stroke="black"))
    c = 0
    while c < m:
        if c in org["strands_cumulative"]:
            dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr-int(org["strands_cumulative"][c]/org["strands_sum"]*100/2*k), mm-c*3+2), stroke_width=3, stroke=blue))
        c += 1
    dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr-int(org["strands_fivetyplus"]/org["strands_sum"]*100/2*k)-1, mm-c*3+2), stroke_width=3, stroke=blue))
    c = 0
    while c < m:
        if c in org["helices_cumulative"]:
            dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr+1+int(org["helices_cumulative"][c]/org["helices_sum"]*100/2*k)*2, mm-c*3+2), stroke_width=3, stroke=orange))
        c += 1
    dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr+1+int(org["helices_fivetyplus"]/org["helices_sum"]*100/2*k)*2, mm-c*3+2), stroke_width=3, stroke=orange))
    l = len(org["org_name"])
    if org["org_name"].endswith(".json"):
        org["org_name"] = org["org_name"][:-5]
    if l > 15:
        orig_name = org["org_name"]
        sname = orig_name.split(" ")
        if len(sname) == 1:
            output_name = sname[0][:5]
        else:
            output_name = sname[0][:5] + ". " + sname[1]
        l = len(output_name)
    else:
        output_name = org["org_name"]

    dwg.add(dwg.text(output_name, insert=(tr+4-int(l*3), mm+10*3+2), color="red", font_size=12, style="text-anchor:center"))


def violin_all(tr, data_list, k, img, mm):
    # draw cumulative violin graph
    m = 50
    ldata_list = len(data_list)

    c = 0
    # draw black borders of strands graph
    while c < m:
        ossum = 0
        for org in data_list:
            if c in org["strands_cumulative"]:
                ossum += org["strands_cumulative"][c]/org["strands_sum"]
        if c in org["strands_cumulative"]:
            dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr-int(ossum/ldata_list*100/2*k)-1, mm-c*3+2), stroke_width=4, stroke="black"))
        else:
            dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr-1, mm-c*3+2), stroke_width=4, stroke="black"))
        c += 1
    # draw 50+ column
    ossum = 0
    for org in data_list:
        ossum += org["strands_fivetyplus"]/org["strands_sum"]
    dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr-int(ossum/ldata_list*100/2*k)-1-1, mm-c*3+2), stroke_width=4, stroke="black"))


    c = 0
    # draw black borders of helices graph
    while c < m:
        ohsum = 0
        for org in data_list:
            if c in org["helices_cumulative"]:
                ohsum += org["helices_cumulative"][c]/org["helices_sum"]
        if c in org["helices_cumulative"]:
            dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr+1+int(ohsum/ldata_list*100/2*k)*2+1, mm-c*3+2), stroke_width=4, stroke="black"))
        else:
            dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr+2, mm-c*3+2), stroke_width=4, stroke="black"))
        c += 1
    # draw 50+ column
    ohsum = 0
    for org in data_list:
        ohsum += org["helices_fivetyplus"]/org["helices_sum"]
    dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr+1+int(ohsum/ldata_list *100/2*k)*2+1, mm-c*3+2), stroke_width=4, stroke="black"))

    c = 0
    # draw blue strands graph
    while c < m:
        ossum = 0
        for org in data_list:
            if c in org["strands_cumulative"]:
                ossum += org["strands_cumulative"][c]/org["strands_sum"]

        if c in org["strands_cumulative"]:
            dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr-int(ossum/ldata_list*100/2*k), mm-c*3+2), stroke_width=3, stroke=blue))
        c += 1
    # draw 50+ column
    ossum = 0
    for org in data_list:
        ossum += org["strands_fivetyplus"]/org["strands_sum"]
    dwg.add(dwg.line(start=(tr, mm-c*3+2), end=(tr-int(ossum/ldata_list*100/2*k)-1, mm-c*3+2), stroke_width=3, stroke=blue))

    c = 0
    # draw orange of helices graph
    while c < m:
        ohsum = 0
        for org in data_list:
            if c in org["helices_cumulative"]:
                ohsum += org["helices_cumulative"][c]/org["helices_sum"]
        if c in org["helices_cumulative"]:
            dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr+1+int(ohsum/ldata_list*100/2*k)*2, mm-c*3+2), stroke_width=3, stroke=orange))
        c += 1
    # draw 50+ column
    ohsum = 0
    for org in data_list:
        ohsum += org["helices_fivetyplus"]/org["helices_sum"]
    dwg.add(dwg.line(start=(tr+1, mm-c*3+2), end=(tr+1+int(ohsum/ldata_list*100/2*k)*2, mm-c*3+2), stroke_width=3, stroke=orange))


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--dirn', required=True, action="store", help="dir containing input jsons")
parser.add_argument('--all', action="store_true", help="dir containing input jsons")
parser.add_argument('--name', action="store", help="dir containing input jsons")
args = parser.parse_args()

# read jsons from the whole dir
data_list = []
rmsd_data_list = {"helices": {}, "strands": {}, "sses": {}, "len_helices": {}, "len_strands": {}, "avg_helices_number": {}, "proteins": {}}
diff_data_list = {}
prot_number = 0
name_list = set()
for filename in os.listdir(args.dirn):
    fnsplit = filename.split("_")
    if not filename.endswith("json"):
        continue
    with open(f"{args.dirn}/{filename}", 'r') as infile:
        data = json.load(infile)
        prot_number += 1
        data_list.append(add_records(data, fnsplit))

width = 1900
dwg = svgwrite.Drawing(args.name, size=(width, 300))
dwg.add(dwg.rect(insert=(1, 1), size=(width, 300), fill="white"))

tr = 100
mm = 200
if args.all:
    tr += 110
    violin_all(tr, data_list, 3, dwg, mm)
else:
    for org in data_list:
        tr += 100
        violin2(tr, 3, org, dwg, mm)

# frame around
dwg.add(dwg.line(start=(150, 45), end=(150, 215), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(tr+50, 45), end=(tr+50, 215), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(150, 45), end=(tr+50, 45), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(150, 215), end=(tr+50, 215), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(150, mm-0*3+2), end=(145, mm-0*3+2), stroke_width=1, stroke="black"))
dwg.add(dwg.text("0", insert=(135, mm-0*3+2), font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(150, mm-10*3+2), end=(145, mm-10*3+2), stroke_width=1, stroke="black"))
dwg.add(dwg.text("10", insert=(131, mm-10*3+2), font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(150, mm-20*3+2), end=(145, mm-20*3+2), stroke_width=1, stroke="black"))
dwg.add(dwg.text("20", insert=(131, mm-20*3+2), font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(150, mm-30*3+2), end=(145, mm-30*3+2), stroke_width=1, stroke="black"))
dwg.add(dwg.text("30", insert=(131, mm-30*3+2), font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(150, mm-40*3+2), end=(145, mm-40*3+2), stroke_width=1, stroke="black"))
dwg.add(dwg.text("40", insert=(131, mm-40*3+2), font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(150, mm-50*3+2), end=(145, mm-50*3+2), stroke_width=1, stroke="black"))
dwg.add(dwg.text("50", insert=(131, mm-50*3+2), font_size=12, style="text-anchor:center"))
dwg.add(dwg.text("organism", insert=(tr/2+75, 250), color="red", font_size=12, style="text-anchor:center"))

# right info box
dwg.add(dwg.line(start=(tr+70, 45), end=(tr+70, 105), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(tr+180, 45), end=(tr+180, 105), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(tr+180, 45), end=(tr+70, 45), stroke_width=1, stroke="black"))
dwg.add(dwg.line(start=(tr+180, 105), end=(tr+70, 105), stroke_width=1, stroke="black"))
dwg.add(dwg.text("sses type", insert=(tr+80, 60), color="red", font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(tr+80, 72), end=(tr+95, 72), stroke_width=4, stroke="black"))
dwg.add(dwg.line(start=(tr+81, 72), end=(tr+94, 72), stroke_width=3, stroke=blue))
dwg.add(dwg.text("strands", insert=(tr+100, 75), color="red", font_size=12, style="text-anchor:center"))
dwg.add(dwg.line(start=(tr+80, 87), end=(tr+95, 87), stroke_width=4, stroke="black"))
dwg.add(dwg.line(start=(tr+81, 87), end=(tr+94, 87), stroke_width=3, stroke=orange))
dwg.add(dwg.text("helicess", insert=(tr+100, 90), color="red", font_size=12, style="text-anchor:center"))

# left text
left = 123
top = 190
dwg.add(dwg.text("number of sses in protein", insert=(left, top), transform=f'rotate(270, {left}, {top})'))


print("save")
dwg.save()
