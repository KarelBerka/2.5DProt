#!/usr/bin/env python3
# tool which generates diagram like:
# http://147.251.21.23/pub/2023_12_05_piny/up_and_down/diagram-up2.png
import argparse
import time
import os
import shutil
import svgwrite
import json
import sys
import pandas as pd

HEAD_TEMPLATE = """
<html>
    <head>
        <title>2DProt status page</title>
        <style>

            table {
                border-collapse: collapse;
                width: 90%;
            }

            td {
                padding: .5em;
            }
        </style>
    </head>
    <body>
    <img src="vzor.png">
    """


def print_template_table_to_file(fd, dir_name, pin_set, short):
    if short:
        print(f"<br>Pin directory {name}<br> pin set {pin_set}", file=fd)
        print("<table style='width:220px'><tr>", file=fd)
        print(f"<td rowspan='2'>1<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H6-avg.svg' width='10'> </td>", file=fd)
        print(f"<td rowspan='2'>2<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H7-avg.svg' width='10'> </td>", file=fd)
        print(f"<td rowspan='2'>3<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H9-avg.svg' width='10'> </td>", file=fd)
        print(f"<td >4b<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H13-avg.svg' width='10'> </td>", file=fd)
        print(f"<td rowspan='2'>5<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H15-avg.svg' width='10'> </td>", file=fd)
        print(f"<td rowspan='2'>10<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H59-avg.svg' width='10'> </td>", file=fd)
        print(f"<td >9a<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H53-avg.svg' width='10'> </td>", file=fd)
        print(f"<td rowspan='2'>8<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H51-avg.svg' width='10'> </td>", file=fd)
        print(f"<td rowspan='2'>7<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H48-avg.svg' width='10'> </td>", file=fd)
        print(f"<td>6b<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H45-avg.svg' width='10'> </td>", file=fd)
        print("</tr>", file=fd)
        print(f"<tr><td>4a<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H11-avg.svg' width='10'></dt>", file=fd)
        print(f"<td>9b<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H56-avg.svg' width='10'></dt>", file=fd)
        print(f"<td>6a<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H44-avg.svg' width='10'></dt></tr>", file=fd)
        print("</table>", file=fd)
    else:
        print(f"<br>Pin directory {name}<br> pin set {pin_set}", file=fd)
        print("<table  style='width:420px'><tr>", file=fd)
        print(f"<td rowspan='2'>1<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H6-avg.svg' width='60'> </td>", file=fd)
        print(f"<td rowspan='2'>2<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H7-avg.svg' width='60'> </td>", file=fd)
        print(f"<td rowspan='2'>3<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H9-avg.svg' width='60'> </td>", file=fd)
        print(f"<td>4b<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H13-avg.svg' width='60'> </td>", file=fd)
        print(f"<td rowspan='2'>5<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H15-avg.svg' width='60'> </td>", file=fd)
        print(f"<td rowspan='2'>10<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H59-avg.svg' width='60'> </td>", file=fd)
        print(f"<td>9a<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H53-avg.svg' width='60'> </td>", file=fd)
        print(f"<td  rowspan='2'>8<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H51-avg.svg' width='60'> </td>", file=fd)
        print(f"<td  rowspan='2'>7<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H48-avg.svg' width='60'> </td>", file=fd)
        print(f"<td>6b<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H45-avg.svg' width='60'> </td>", file=fd)
        print("</tr>", file=fd)
        print(f"<tr><td>4a<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H11-avg.svg' width='60'></dt>", file=fd)
        print(f"<td>9b<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H56-avg.svg' width='60'></dt>", file=fd)
        print(f"<td>6a<br><img src='{pin_set}/firm-{dir_name}-{pin_set}-H44-avg.svg' width='60'></dt></tr>", file=fd)
        print("</table>", file=fd)

#    Pin directory {name}<br> pin set {args.pin_set}


def draw_avg_sse_diagram(prob_data, occurence, sse_name, fdir, html_file, up_down, name, pin_set, cs, short):
    # print sse
    sse_max = 0
    while sse_max < len(prob_data[sse_name]):
        if occurence[sse_name][sse_max] == 0:
             break
        if occurence[sse_name][0] / occurence[sse_name][sse_max] > 10:
             break
        sse_max = sse_max+1

    if short:
        width = 6
    else:
        width = 24

    if short:
        h = 3
    else:
        h = 6
    height = h*(sse_max)+2

    dwg = svgwrite.Drawing(f"{fdir}/{pin_set[0]}/firm-{name}-{pin_set[0]}-{sse_name}-avg.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(0, 0), size=(width, height), fill="white"))
    c = 0
    for i in prob_data[sse_name]:
        if i >= sse_max:
            break
        c += 1

        shapes = dwg.add(dwg.g(id='shapes'))
        op = prob_data[sse_name][i][1]

        if up_down == "down":
          y_coord = c*(h-1)+1
          if cs == "rg":
           if op > 50:
            color = "green"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=op/50-1, stroke=color))
           else:
            color = "red"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=1-op/50, stroke=color))
          elif cs == "by":
           if op > 50:
            color = "blue"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=op/50-1, stroke=color))
           else:
            color = "yellow"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=1-op/50, stroke=color))
          else:
           if op > 50:
            color = "blue"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=op/50-1, stroke=color))
           else:
            color = "red"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=1-op/50, stroke=color))

        elif up_down == "up":
          y_coord = sse_max*h - c*h+1
          if cs == "rg":
           if op > 50:
            color = "green"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=op/50-1, stroke=color))
           else:
            color = "red"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=1-op/50, stroke=color))
          elif cs == "by":
           if op > 50:
            color = "blue"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=op/50-1, stroke=color))
           else:
            color = "yellow"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=1-op/50, stroke=color))
          else:
           if op > 50:
            color = "blue"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=op/50-1, stroke=color))
           else:
            color = "red"
            shapes.add(dwg.rect(insert=(0, y_coord), size=(width, h-1), fill=color, opacity=1-op/50, stroke=color))

        if not short:
            shapes.add(dwg.text(f"{prob_data[sse_name][i][1]} ({prob_data[sse_name][i][0]})", insert=(3, y_coord + 3), font_size=4))
    dwg.save()


def enum_list(r_list, r_all):
    cc = {}
    for i in r_list:
        if i not in cc:
            cc[i] = 0
        cc[i] += 1
    dd = {}
    for i in cc:
        cc[i] = int(cc[i]/r_all*100)
        if cc[i] >= 10:
            dd[i] = cc[i]
    mmax = 0
    mrec = ""
    for i in dd:
         if dd[i] > mmax:
             mmax = dd[i]
             mrec = i
    return cc, mrec, mmax, dd


def l_diff(l1, l2):
    ld = {}
    for i in l1:
        if i in l2:
            ld[i] = l1[i]-l2[i]
        else:
            ld[i] = l1[i]
    for i in l2:
        if i not in l1:
            ld[i] = -l2[i]

    return {k: v for k, v in sorted(ld.items(), key=lambda item:item[1])}


def draw_avg_sse_diagram_with_limit(pin_set, new_dir, html_file, sses_list, name, occurence, short):
    aux = {}
    oc = {}
    
    aaa = 0
    acount = 0

    for sse in sses_list:
        aux[sse] = {}
        oc[sse] = {}
        if sse in data:
            len_max = 0
            for pin in pin_set:
                if pin in occurence[sse][0]:
                    len_max += occurence[sse][0][pin]

            for ind in data[sse]:
                aux[sse][ind] = {}
                oc[sse][ind] = 0
                rec = []
                for pin in pin_set:
                    if pin in data[sse][ind]:
                        rec += data[sse][ind][pin]
                        oc[sse][ind] += occurence[sse][ind][pin]
                d, r, m, c = enum_list(rec, len_max)
                if pin_set[0] in occurence[sse][ind]:
                    aaa += occurence[sse][ind][pin_set[0]]*m
                    acount += occurence[sse][ind][pin_set[0]]
                aux[sse][ind] = [r,m]
            if sses_dict[sse] in ["1", "3", "5", "7", "9a", "9b"]:
                up_down = "down"
            elif sses_dict[sse] in ["2", "4a", "4b", "6a", "6b", "8", "10"]:
                up_down = "up"
            draw_avg_sse_diagram(aux, oc, sse, new_dir, html_file, up_down, name, pin_set, "rg", short)
    print("aaaaaaaaaaa",pin_set,  aaa, acount, aaa/acount)

def draw_dif_sse_diagram_with_limit(pin_set1, pin_set2, new_dir, html_file, sses_list, name, occurence, short):
    aux1 = {}
    aux2 = {}
    oc1 = {}
    oc2 = {}
    aux = {}
    diffs = []
    sssum = 0
    sscount = 0

    for sse in sses_list:
        aux1[sse] = {}
        oc1[sse] = {}
        aux2[sse] = {}
        oc2[sse] = {}
        aux[sse] = {}
        if sse in data:

            len_max1 = 0
            for pin in pin_set1:
                if pin in occurence[sse][0]:
                    len_max1 += occurence[sse][0][pin]

            len_max2 = 0
            for pin in pin_set2:
                if pin in occurence[sse][0]:
                    len_max2 += occurence[sse][0][pin]
            len_max = min(len_max1, len_max2)

            for ind in data[sse]:
                aux1[sse][ind] = {}
                oc1[sse][ind] = 0
                rec = []
                for pin in pin_set1:
                    if pin in data[sse][ind]:
                        rec += data[sse][ind][pin]
                        oc1[sse][ind] += occurence[sse][ind][pin]
                d, r, m,c  = enum_list(rec, len_max1)

                aux1[sse][ind] = d

            for ind in data[sse]:
                aux2[sse][ind] = {}
                oc2[sse][ind] = 0
                rec = []
                for pin in pin_set2:
                    if pin in data[sse][ind]:
                        rec += data[sse][ind][pin]
                        oc2[sse][ind] += occurence[sse][ind][pin]
                d, r, m, c = enum_list(rec, len_max2)
                aux2[sse][ind] = d

                ssum = 0
                txt = ""
                for i in aux1[sse][ind]:
                    if i in aux2[sse][ind]:
                        ssum += min(aux1[sse][ind][i], aux2[sse][ind][i])
                        if min(aux1[sse][ind][i], aux2[sse][ind][i]) >= 1:
                            txt += f"{int(min(aux1[sse][ind][i], aux2[sse][ind][i]))},"
            
                #print(f"{sse}, {ind}, {pin_set1[0]} {pin_set2[0]} {data[sse][ind].keys()}")
                if (pin_set1[0] in data[sse][ind]) and (pin_set2[0] in data[sse][ind]):
                    d = l_diff(enum_list(data[sse][ind][pin_set1[0]], len(data[sse][ind][pin_set1[0]]))[3], enum_list(data[sse][ind][pin_set2[0]], len(data[sse][ind][pin_set2[0]]))[3])
#                        print("2 chyba sse", ind, sse, pin_set1[0], pin_set2[0], "ssum", aux2[sse][ind][amin_acid],aux2[sse][ind][amin_acid])
                    print("rozdil", d)
                    diffs.append(d)
                    if (occurence[sse][ind][pin_set1[0]]*100/occurence[sse][0][pin_set1[0]] > 10) and (occurence[sse][ind][pin_set2[0]]*100/occurence[sse][0][pin_set2[0]] > 10):
                        print("ssum", ssum)
                        sssum += ssum
                        sscount += 1
#                        print(f"{enum_list(data[sse][ind][pin_set1[0]], len(data[sse][ind][pin_set1[0]]))[3]} bude {enum_list(data[sse][ind][pin_set2[0]], len(data[sse][ind][pin_set2[0]]))[3]}")

                aux[sse][ind] = [txt, ssum]

            if sses_dict[sse] in ["1", "3", "5", "7", "9a", "9b"]:
                up_down = "down"
            elif sses_dict[sse] in ["2", "4a", "4b", "6a", "6b", "8", "10"]:
                up_down = "up"
            draw_avg_sse_diagram(aux, oc1, sse, new_dir, html_file, up_down, name, [pin_set1[0]+"_"+pin_set2[0]], "by",short)
    asum = {}
    aasum = {}
    print("******************", sssum, sscount, sssum/sscount)
    for i in diffs:
        for j in i:
            if j not in asum:
                asum[j] = 0
                aasum[j] = 0
#            print("ij", i,j)
            asum[j] += i[j]
            aasum[j] += int(abs(i[j]))
    for j in aasum:
        aasum[j] = (aasum[j] - asum[j]) / 2
    
    print("asum -rozdil", {k: v for k, v in sorted(asum.items(), key=lambda item:item[0])})
    print("aasum - rozdil", {k: v for k, v in sorted(aasum.items(), key=lambda item:item[0])})



def draw_dif_sse_diagram_with_limit_with_aa(pin_set1, pin_set2, new_dir, html_file, sses_list, name, occurence, amin_acid, short):
    aux1 = {}
    aux2 = {}
    oc1 = {}
    oc2 = {}
    aux = {}
    diffs = []

    for sse in sses_list:
        aux1[sse] = {}
        oc1[sse] = {}
        aux2[sse] = {}
        oc2[sse] = {}
        aux[sse] = {}
        if sse in data:

            len_max1 = 0
            for pin in pin_set1:
                if pin in occurence[sse][0]:
                    len_max1 += occurence[sse][0][pin]

            len_max2 = 0
            for pin in pin_set2:
                if pin in occurence[sse][0]:
                    len_max2 += occurence[sse][0][pin]
            len_max = min(len_max1, len_max2)

            for ind in data[sse]:
                aux1[sse][ind] = {}
                oc1[sse][ind] = 0
                rec = []
                for pin in pin_set1:
                    if pin in data[sse][ind]:
                        rec += data[sse][ind][pin]
                        oc1[sse][ind] += occurence[sse][ind][pin]
                d, r, m,c  = enum_list(rec, len_max1)
                print("chyba sse", ind, sse, amin_acid, pin_set1[0], "d ---", d, "r", r, "m", m, "--")

                aux1[sse][ind] = d

            for ind in data[sse]:
                aux2[sse][ind] = {}
                oc2[sse][ind] = 0
                rec = []
                for pin in pin_set2:
                    if pin in data[sse][ind]:
                        rec += data[sse][ind][pin]
                        oc2[sse][ind] += occurence[sse][ind][pin]
                d, r, m, c = enum_list(rec, len_max2)
#                print("chyba sse", ind, sse, amin_acid, pin_set2[0], "d ---", d, "r", r, "m", m, "--")

                aux2[sse][ind] = d

                txt = ""
                if amin_acid in aux1[sse][ind]:
                    if amin_acid in aux2[sse][ind]:
#                        print("0 chyba sse", ind, sse, pin_set1[0], pin_set2[0], "ssum", aux1[sse][ind][amin_acid] - aux2[sse][ind][amin_acid], aux1[sse][ind][amin_acid],aux2[sse][ind][amin_acid])
                        
                        ssum = aux1[sse][ind][amin_acid] - aux2[sse][ind][amin_acid]
                    else:
#                        print("1 chyba sse", ind, sse, pin_set1[0], pin_set2[0], "ssum", aux1[sse][ind][amin_acid], aux1[sse][ind][amin_acid])
                        ssum = aux1[sse][ind][amin_acid]
                else:
                    if amin_acid in aux2[sse][ind]:
#
                        ssum = -aux2[sse][ind][amin_acid]
                    else:
#                        print("3 chyba sse", ind, sse, pin_set1[0], pin_set2[0], "ssum", 0)
                        ssum = 0
#                print("hhhh", sse, ind, ssum, end=":")
#                if amin_acid in aux1[sse][ind]:
#                    print("1 - ", aux1[sse][ind][amin_acid], end=", ")
#                if amin_acid in aux2[sse][ind]:
#                    print("2 - ", aux2[sse][ind][amin_acid], end="")
#                print()
                aux[sse][ind] = [amin_acid[0] + str(ssum), (ssum+100)/2]
#                print("x chyba sse", ind, sse, pin_set1[0], pin_set2[0], "res", aux[sse][ind])
#                print("ssum", ssum)
#                draw_avg_sse_diagram(ssum, oc1, sse, new_dir, html_file, up_down, name, [pin_set1[0]+"_"+pin_set2[0]], "by",short)

            if sses_dict[sse] in ["1", "3", "5", "7", "9a", "9b"]:
                up_down = "down"
            elif sses_dict[sse] in ["2", "4a", "4b", "6a", "6b", "8", "10"]:
                up_down = "up"
            draw_avg_sse_diagram(aux, oc1, sse, new_dir, html_file, up_down, name, [pin_set1[0]+"_"+pin_set2[0]+"_"+amin_acid], "br",short)


def stats_with_limit(pin_set, sses_list):
    r_acids = {}
    r_sum = 0

    for sse in sses_list:
        if sse in data:
            for ind in data[sse]:
                for pin in pin_set:
                    if pin in data[sse][ind]:
                        for aa in data[sse][ind][pin]:
                            if aa not in r_acids:
                                r_acids[aa] = 0
                            r_acids[aa] += 1
                            r_sum += 1
    print("rrr", r_sum)
    for aa in r_acids:
        r_acids[aa] = r_acids[aa] * 100 / r_sum
    print("rrr", r_acids)


sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]   # BIG DATASET
sses_dict = {"H6": "1", "H7": "2", "H9": "3", "H11": "4a", "H13": "4b", "H15": "5", "H44": "6a", "H45": "6b", "H48": "7", "H51": "8", "H53": "9a", "H56": "9b", "H59": "10"}  # BIG DATASETsses_dict = {"H6": "1", "H7": "2", "H9": "3", "H11": "4a", "H13": "4b", "H15": "5", "H44": "6a", "H45": "6b", "H48": "7", "H51": "8", "H53": "9a", "H56": "9b", "H59": "10"}  # BIG DATASET

# read arguments
parser = argparse.ArgumentParser(
    description='generated residues per sses diagram')
parser.add_argument('--dir', required=True, action="store",
                    help="name of t")
parser.add_argument('--charges', required=True, action="store",
                    help="name of charges dir")
parser.add_argument('--pin_set', required=True, action="store",
                    help="name of pin subdir")
parser.add_argument('--cmp_pin_set', required=False, action="store",
                    help="name of pin subdir to comarison")
parser.add_argument('--short', required=False, action="store_true",
                    help="short output")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = os.path.basename(os.path.normpath(args.dir))

if args.short:
    short = True
else:
    short = False


if args.cmp_pin_set:
    new_dir = f"res_graph-types-{name}-{args.pin_set}-{args.cmp_pin_set}-{timestamp}"
    cmp_pin_set = args.cmp_pin_set
else:
    new_dir = f"res_graph-types-{name}-{args.pin_set}-{timestamp}"
    cmp_pin_set = ""

os.makedirs(new_dir, exist_ok=True)
print("new_dir", new_dir)

fn_order = []
pin_mapping = {}

data = {}        # data[sse][res_number][pin]
occurence = {}   # occurence[sse][res_numer][pin]

# read input data
for dirname in os.listdir(args.dir):
    if os.path.isdir(dirname) and ((dirname == args.pin_set) or (dirname == args.cmp_in_set)):
        continue
    print("pin directory name:", dirname[0:4], ": ", dirname)
    for filename in os.listdir(f"{args.dir}/{dirname}"):
        sfilename = filename.split("-")
        if dirname[0:4] not in pin_mapping:
            pin_mapping[dirname[0:4]] = [filename]
        else:
            pin_mapping[dirname[0:4]].append(filename)
        # go through all PINs in PIN directory and create record for them

        # read PIN json and read residues charges in SSEs
        with open(f"{args.dir}/{dirname}/{filename}", "r") as fd:
          print("read", f"{args.dir}/{dirname}/{filename}")
          jdata = json.load(fd)
          fn_order.append(filename)
          f_charges = f"{args.charges}/{sfilename[1]}.alphacharges"

          with open(f"{args.dir}/{dirname}/{filename}", "r") as fd:
            for sse in jdata["helices_residues"]:
                if sse in sses_list:
                    if sse not in data:
                        data[sse] = {}
                        occurence[sse] = {}
                else:
                    continue

                with open(f_charges, "rb") as fc:
                    line = fc.readline()
                    counter = 0
                    processed = []
                    while line:
                        sline = str(line)
                        counter += 1
                        if sline.startswith("TER"):
                            break
                        elif (sline[2:6] == "ATOM"):
                            ind = int(sline[24:28])
                            if (ind not in processed) and (ind >= jdata["helices_residues"][sse][1]) and (ind <= jdata["helices_residues"][sse][0]):
                                processed.append(ind)
                                if (ind - jdata["helices_residues"][sse][1]) not in data[sse]:
                                    data[sse][ind - jdata["helices_residues"][sse][1]] = {}
                                    occurence[sse][ind - jdata["helices_residues"][sse][1]] = {}

                                if dirname[0:4] not in data[sse][ind - jdata["helices_residues"][sse][1]]:
                                   data[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]] = []
                                   occurence[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]] = 0

                                if sline[19:22] in ["ALA", "VAL", "LEU", "ILE", "MET", "TYR", "PHE", "TRP" ]:
                                    #nonpolar
                                    txt = "hydro"
                                elif sline[19:22] in ["SER", "THR", "ASN", "GLN"]:
                                    #polar
                                    txt = "polar"
                                elif sline[19:22] in ["SEC", "GLY", "PRO", "CYS"]:
                                    # special cases
                                    txt = "special"
                                elif sline[19:22] in ["ASP", "GLU"]:
                                    #negativly charged
                                    txt = "neg"
                                elif sline[19:22] in ["LYS", "ARG", "HIS"]:
                                    #positively charged
                                    txt = "pos"
                                else:
                                   print("neni uveden ", sline[19:22])#, i , name)
                                   sys.exit()

                                data[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]].append(txt)
                                occurence[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]] += 1
                            if (ind > jdata["helices_residues"][sse][0]):
                                break

                        line = fc.readline()

# do stats graphs
shutil.copyfile("vzor.png",f"{new_dir}/vzor.png")
html_file = f"{new_dir}/index.html"

# red-green stability graph for pin_set
with open(html_file, "w") as fd:
    print(HEAD_TEMPLATE, file=fd)
#    print(f"Pin directory {name}<br> pin set {args.pin_set}", file=fd)
    print("<table border='1' colspan='2' ><tr><td> each residuum is represented by number - how much percent of proteins in PIN group have 'the major' amino-acid type.<br>", file=fd)
    print(" E.g. if for protein set PIN1 the first res in SSE '1' is in 40% LEU, 35% PHE and 25% THR - on the first res of SSE '1' will be '40% (LEU)'<br>", file=fd)
    print("more stable res (with the major amino-acid probability >50% is dark green, unstable <50% is <style='background-color:red;'>red</style>)<br>", file=fd)
    print("0% usually means the res is only in small percentage of proteins </td></tr></table>", file=fd)
    print("<table border='1'><tr><td>", file=fd)
    print_template_table_to_file(fd, name, args.pin_set,short)
os.makedirs(f"{new_dir}/{args.pin_set}", exist_ok=True)
draw_avg_sse_diagram_with_limit([args.pin_set], new_dir, html_file, sses_list, name, occurence, short)


if args.cmp_pin_set:
    # red-green stability graph for cmp_pin_set
    with open(html_file, "a") as fd:
        print("</td><td>", file=fd)
        print_template_table_to_file(fd, name, args.cmp_pin_set, short)
    os.makedirs(f"{new_dir}/{args.cmp_pin_set}", exist_ok=True)
    draw_avg_sse_diagram_with_limit([args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, short)

    # diff graph for pin_Set and cmp_pin_set
    with open(html_file, "a") as fd:
        if short:
            print("</td><td>", file=fd)
        else:
            print("</td></tr>", file=fd)
            print("<table border='1' colspan='2' ><tr><td> the difference of the first pin group and second pin group<br>", file=fd)
            print(f"the number in res means how much are {args.pin_set} and {args.cmp_pin_set} similar - in percentage<br>", file=fd)
            print (" big number (blue)100% means both have the same amino-acids - e.g. 50%ILE and 50%ALA<br>", file=fd)
            print("the low number (yellow) means they are different (e.g. 10%  - for example situation in first pin group there is 10% ALA and 90% VAL, the second PIN group have 100% ALA on the same position<br>", file=fd)
            print("</td></tr>", file=fd)
            print("<tr align='center'><td colspan='2'>", file=fd)
            
        print_template_table_to_file(fd, name, args.pin_set + "_"+args.cmp_pin_set, short)
    os.makedirs(f"{new_dir}/{args.pin_set}_{args.cmp_pin_set}", exist_ok=True)
    draw_dif_sse_diagram_with_limit([args.pin_set], [args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, short)
    with open(html_file, "a") as fd:
        if not short:
            short = True
            print("</td></tr></table>", file=fd)
            print("<table border='1'bgcolor='grey'><tr><td colspan='5' >", file=fd)
            print(f"each graph is if different amino-acid. If the res is blue residuum is in {args.pin_set} and not in {args.cmp_pin_set} ", file=fd)
            print("it is red if it is otherwise", file=fd)
            print(f"it is white if the probability in {args.pin_set} and {args.cmp_pin_set} is the same", file=fd)
            print("</td></tr>", file=fd)
            print("<tr align='center'><td colspan='2'>", file=fd)
            
#            print("<tr align='center'><td colspan='2'>", file=fd)
            amin_acid = "special"
            print_template_table_to_file(fd, name, args.pin_set + "_"+args.cmp_pin_set + "_" + amin_acid, short)
            os.makedirs(f"{new_dir}/{args.pin_set}_{args.cmp_pin_set}_{amin_acid}", exist_ok=True)
            draw_dif_sse_diagram_with_limit_with_aa([args.pin_set], [args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, amin_acid, short)
            print("</td><td>", file=fd)
#            print("<tr align='center'><td colspan='2'>", file=fd)
            amin_acid = "hydro"
            print_template_table_to_file(fd, name, args.pin_set + "_"+args.cmp_pin_set + "_" + amin_acid, short)
            os.makedirs(f"{new_dir}/{args.pin_set}_{args.cmp_pin_set}_{amin_acid}", exist_ok=True)
            draw_dif_sse_diagram_with_limit_with_aa([args.pin_set], [args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, amin_acid, short)
            print("</td><td>", file=fd)
#            print("<tr align='center'><td colspan='2'>", file=fd)
            amin_acid = "polar"
            print_template_table_to_file(fd, name, args.pin_set + "_"+args.cmp_pin_set + "_" + amin_acid, short)
            os.makedirs(f"{new_dir}/{args.pin_set}_{args.cmp_pin_set}_{amin_acid}", exist_ok=True)
            draw_dif_sse_diagram_with_limit_with_aa([args.pin_set], [args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, amin_acid, short)
            print("</td><td>", file=fd)
#            print("<tr align='center'><td colspan='2'>", file=fd)
            amin_acid = "neg"
            print_template_table_to_file(fd, name, args.pin_set + "_"+args.cmp_pin_set + "_" + amin_acid, short)
            os.makedirs(f"{new_dir}/{args.pin_set}_{args.cmp_pin_set}_{amin_acid}", exist_ok=True)
            draw_dif_sse_diagram_with_limit_with_aa([args.pin_set], [args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, amin_acid, short)
            print("</td><td>", file=fd)
#            print("<tr align='center'><td colspan='2'>", file=fd)
            amin_acid = "pos"
            print_template_table_to_file(fd, name, args.pin_set + "_"+args.cmp_pin_set + "_" + amin_acid, short)
            os.makedirs(f"{new_dir}/{args.pin_set}_{args.cmp_pin_set}_{amin_acid}", exist_ok=True)
            draw_dif_sse_diagram_with_limit_with_aa([args.pin_set], [args.cmp_pin_set], new_dir, html_file, sses_list, name, occurence, amin_acid, short)

else:
    with open(html_file, "a") as fd:
        print("</td></tr></table>", file=fd)

print("")

#stats_with_limit([args.pin_set], sses_list)
#
#stats_with_limit_cmp([args.pin_set],[args.cmp_pin_set], sses_list)



