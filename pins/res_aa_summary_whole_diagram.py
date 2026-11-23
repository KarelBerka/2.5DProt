#!/usr/bin/env python3
# tool which generates diagram like:
# http://147.251.21.23/pub/2023_12_05_piny/up_and_down/diagram-up2.png

import argparse
import time
import os
import svgwrite
import json


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


def draw_scale(output_dir):
    print("draw_scale")
    width = 60
    height = 130

    dwg = svgwrite.Drawing(f"{output_dir}/scale_rg.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(0, 0), size=(width, height), fill="white"))

    for i in range(11):
        shapes = dwg.add(dwg.g(id='shapes'))
        print(i*10)
        op = i*10
        if op > 50:
            color = "green"
            shapes.add(dwg.rect(insert=(10, i*10+10), size=(20, 9), fill=color, opacity=op/50-1, stroke=color))
            print(i*10, op/50-1)

        else:
            color = "red"
            shapes.add(dwg.rect(insert=(10, i*10+10), size=(20, 9), fill=color, opacity=1-op/50, stroke=color))
            print(i*10, 1-op/50)
        shapes.add(dwg.text(f"{i*10}%", insert=(33, i*10+18), font_size=8))
    dwg.save()

    dwg = svgwrite.Drawing(f"{output_dir}/scale_by.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(0, 0), size=(width, height), fill="white"))

    for i in range(11):
        shapes = dwg.add(dwg.g(id='shapes'))
        print(i*10)
        op = i*10
        if op > 50:
            color = "blue"
            shapes.add(dwg.rect(insert=(10, i*10+10), size=(20, 9), fill=color, opacity=op/50-1, stroke=color))
            print(i*10, op/50-1)

        else:
            color = "yellow"
            shapes.add(dwg.rect(insert=(10, i*10+10), size=(20, 9), fill=color, opacity=1-op/50, stroke=color))
            print(i*10, 1-op/50)
        shapes.add(dwg.text(f"{i*10}%", insert=(33, i*10+18), font_size=8))
    dwg.save()


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
            draw_avg_sse_diagram(aux, oc1, sse, new_dir, html_file, up_down, name, [pin_set1[0]+"_"+pin_set2[0]], "by", short)
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

    print("asum -rozdil", {k: v for k, v in sorted(asum.items(), key=lambda item: item[0])})
    print("aasum - rozdil", {k: v for k, v in sorted(aasum.items(), key=lambda item: item[0])})


def draw_dif_sse_diagram_with_limit_with_aa(pin_set1, pin_set2, new_dir, html_file, sses_list, name, occurence, amin_acid, short):
    aux1 = {}
    aux2 = {}
    oc1 = {}
    oc2 = {}
    aux = {}

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

            for ind in data[sse]:
                aux1[sse][ind] = {}
                oc1[sse][ind] = 0
                rec = []
                for pin in pin_set1:
                    if pin in data[sse][ind]:
                        rec += data[sse][ind][pin]
                        oc1[sse][ind] += occurence[sse][ind][pin]
                d, r, m, c = enum_list(rec, len_max1)
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
                aux2[sse][ind] = d

                if amin_acid in aux1[sse][ind]:
                    if amin_acid in aux2[sse][ind]:
                        ssum = aux1[sse][ind][amin_acid] - aux2[sse][ind][amin_acid]
                    else:
                        ssum = aux1[sse][ind][amin_acid]
                else:
                    if amin_acid in aux2[sse][ind]:
                        ssum = -aux2[sse][ind][amin_acid]
                    else:
                        ssum = 0
                aux[sse][ind] = [amin_acid[0] + str(ssum), (ssum+100)/2]

            if sses_dict[sse] in ["1", "3", "5", "7", "9a", "9b"]:
                up_down = "down"
            elif sses_dict[sse] in ["2", "4a", "4b", "6a", "6b", "8", "10"]:
                up_down = "up"
            draw_avg_sse_diagram(aux, oc1, sse, new_dir, html_file, up_down, name, [pin_set1[0]+"_"+pin_set2[0]+"_"+amin_acid], "br", short)


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


def compute_maximal_occurence(data_list):
#    print("---1", data_list)
    aa_sum = len(data_list)
#    print("---2", aa_sum)
    occ_list = {}
    for aa in data_list:
        if aa not in occ_list:
            occ_list[aa] = 1
        else:
            occ_list[aa] += 1
    max_aa = 0
    most_aa = ""
#    print("---3", occ_list)
    for aa in occ_list:
       if occ_list[aa] > max_aa:
           max_aa = occ_list[aa]
           most_aa = aa
#    print("---4", most_aa, max_aa)
#    sys.exit()
    return most_aa, int(max_aa/aa_sum*100)
    

sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]   # BIG DATASET
sses_dict = {"H6": "1", "H7": "2", "H9": "3", "H11": "4a", "H13": "4b", "H15": "5", "H44": "6a", "H45": "6b", "H48": "7", "H51": "8", "H53": "9a", "H56": "9b", "H59": "10"}  # BIG DATASETsses_dict = {"H6": "1", "H7": "2", "H9": "3", "H11": "4a", "H13": "4b", "H15": "5", "H44": "6a", "H45": "6b", "H48": "7", "H51": "8", "H53": "9a", "H56": "9b", "H59": "10"}  # BIG DATASET

# read arguments
parser = argparse.ArgumentParser(
    description='generated residues per sses diagram')
parser.add_argument('--dir', required=True, action="store",
                    help="name of subdir with PINx direftories from --(cmp_)pin_set e.g. layouts/pins_diminik")
parser.add_argument('--charges', required=True, action="store",
                    help="name of charges dir e.g. charges")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = os.path.basename(os.path.normpath(args.dir))
new_dir = f"amino_acid_summary_graph-{name}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)
print("new_dir", new_dir)

fn_order = []
pin_mapping = {}

data = {}        # data[sse][res_number][pin]
occurence = {}   # occurence[sse][res_numer][pin]

# read input data
for dirname in os.listdir(args.dir):
    print("dirname", dirname)
    if os.path.isdir(dirname) and ((dirname == args.pin_set) or (dirname == args.cmp_in_set)):
        continue
    print("pin directory name:", dirname[0:4], ": ", dirname)
    cc = 0
    for filename in os.listdir(f"{args.dir}/{dirname}"):
        print("counter", dirname, cc, end = ": ")
        cc += 1
#        if cc > 5:
#            break
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

                                if sline[19:22] == "ARG":
                                     txt="R"
                                elif sline[19:22] == "HIS":
                                     txt="H"
                                elif sline[19:22] == "LYS":
                                     txt="K"
                                elif sline[19:22] == "ASP":
                                     txt="D"
                                elif sline[19:22] == "GLU":
                                     txt="E"
                                elif sline[19:22] == "SER":
                                     txt="S"
                                elif sline[19:22] == "THR":
                                     txt="T"
                                elif sline[19:22] == "ASN":
                                     txt="N"
                                elif sline[19:22] == "GLN":
                                     txt="Q"
                                elif sline[19:22] == "CYS":
                                     txt="C"
                                elif sline[19:22] == "GLY":
                                     txt="G"
                                elif sline[19:22] == "PRO":
                                     txt="P"
                                elif sline[19:22] == "ALA":
                                     txt="A"
                                elif sline[19:22] == "VAL":
                                     txt="V"
                                elif sline[19:22] == "ILE":
                                     txt="I"
                                elif sline[19:22] == "LEU":
                                     txt="L"
                                elif sline[19:22] == "MET":
                                     txt="M"
                                elif sline[19:22] == "PHE":
                                     txt="F"
                                elif sline[19:22] == "TYR":
                                     txt="Y"
                                elif sline[19:22] == "TRP":
                                     txt="W"
                                
#                                if sline[19:22] in ["ALA", "VAL", "LEU", "ILE", "MET", "TYR", "PHE", "TRP"]:
#                                    # nonpolar
#                                    txt = "hydro"
#                                elif sline[19:22] in ["SER", "THR", "ASN", "GLN"]:
#                                    # polar
#                                    txt = "polar"
#                                elif sline[19:22] in ["SEC", "GLY", "PRO", "CYS"]:
#                                    # special cases
#                                    txt = "special"
#                                elif sline[19:22] in ["ASP", "GLU"]:
#                                    # negativly charged
##                                    txt = "neg"
#                                elif sline[19:22] in ["LYS", "ARG", "HIS"]:
#                                    # positively charged
#                                    txt = "pos"
#                                else:
#                                    print("neni uveden ", sline[19:22])
#                                    sys.exit()

#                                data[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]].append(sline[19:22])
                                data[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]].append(txt)
                                occurence[sse][ind - jdata["helices_residues"][sse][1]][dirname[0:4]] += 1
                            if (ind > jdata["helices_residues"][sse][0]):
                                break

                        line = fc.readline()
ratio={}
for sse  in  data:
    ratio[sse] = {}
    for number in data[sse]:
        ratio[sse][number] = {}
        for pin in data[sse][number]:
            ratio[sse][number][pin] = len(data[sse][number][pin])

#print(ratio['H6'])
#print(data['H6'])
#print(data.keys())

#sys.exit()

print("data", data.keys())
#for i in data.keys():
i="H6"

#print("i", i, data[i])
#print("occurence", i, occurence[i])

width = 800
height = 500
dwg = svgwrite.Drawing(f"{new_dir}/graph.svg", size=(width, height))
dwg.add(dwg.rect(insert=(0,0), size=(width, height), fill="white"))

y = 0
limit = 90
limit1 = 89
limit2 = 173
limit3 = 260
sameness = 75
for sse in data.keys():
    y+=2
    for res in data[sse]:
         print(f"res: {sse} {res:3} ", end=" ")
         print("---", data[sse][res].keys())
         # display only sses residues, in which at least one pin have >10% domains
         should_be_shown = False
         for pin in ratio[sse][res]:
             if ratio[sse][res][pin]*10 > ratio[sse][0][pin]:
                 should_be_shown = True
         if not should_be_shown:
             continue

         if "PIN1" in data[sse][res]:
             max_aa, aa = compute_maximal_occurence(data[sse][res]['PIN1'])
             if aa > sameness:
                 val = max_aa
                 the_same = True
             else:
                 val = "."
                 the_same = False
         else:
             val = "."
             the_same = False
#         print("the same", the_same)
         aa_list = []
         for pin_type in data[sse][res]:
             max_aa, aa = compute_maximal_occurence(data[sse][res][pin_type])
             if (aa > sameness):
                 aa_list.append(max_aa)
             if (max_aa != val) or (aa < sameness):
#                  print("zmena", val, max_aa, max_aa ==  data[sse][res], sse, res)
                 the_same = False
         if len(aa_list) > 0:
             cc_max, cc = compute_maximal_occurence(aa_list)
         else:
             cc = 0
             cc_max = []
         print("the same", the_same, cc_max, cc, aa_list, len(aa_list))
         if ((cc > 85) and len(aa_list)==8) or ((cc >90 ) and (len(aa_list)>=7)):
             print("ruzovy", end=" ")
             if (y > limit1) and (y <= limit2):
                 print("2")
                 dwg.add(dwg.rect(insert=(10+y*7+10-limit1*7-1,30+11*8), size=(7,8*8), fill="pink"))
             elif (y > limit2) and (y<= limit3):
                 print("3")
                 dwg.add(dwg.rect(insert=(10+y*7+10-limit2*7-1,30+11*8*2), size=(7,8*8), fill="pink"))
             elif (y > limit3):
                 print("4")
                 dwg.add(dwg.rect(insert=(10+y*7+10-limit3*7-1,30+11*8*3), size=(7,8*8), fill="pink"))
             else:
                 print("1")
                 dwg.add(dwg.rect(insert=(10+y*7+10-1,30), size=(7,8*8), fill="pink"))

         if the_same:
             print("cerveny", end=" ")
             if (y > limit1) and (y <= limit2):
                 print("2")
                 dwg.add(dwg.rect(insert=(10+y*7+10-limit1*7-1,30+11*8), size=(7,8*8), fill="red"))
             elif (y > limit2) and (y<= limit3):
                 print("3")
                 dwg.add(dwg.rect(insert=(10+y*7+10-limit2*7-1,30+11*8*2), size=(7,8*8), fill="red"))
             elif (y > limit3):
                 print("4")
                 dwg.add(dwg.rect(insert=(10+y*7+10-limit3*7-1,30+11*8*3), size=(7,8*8), fill="red"))
             else:
                 print("1")
                 dwg.add(dwg.rect(insert=(10+y*7+10-1,30), size=(7,8*8), fill="red"))
         y += 1
#sys.exit()
x = 0
y = 0

for sse in data.keys():
    y+=2
    for res in data[sse]:
#         print(f"res: {res:3}", end=" ")
         # display only sses residues, in which at least one pin have >10% domains
         should_be_shown = False
         for pin in ratio[sse][res]:
             if ratio[sse][res][pin]*10 > ratio[sse][0][pin]:
                 should_be_shown = True
         if not should_be_shown:
             continue
         for pin_type in data[sse][res]:
             x += 1
             max_aa, aa = compute_maximal_occurence(data[sse][res][pin_type])
             print(sse, res, pin_type,":", max_aa, aa)
             if (y > limit1) and (y <= limit2):
                 if aa > sameness:
#                     print(f"{max_aa}", end="" )
                     dwg.add(dwg.text(max_aa, insert=(10+y*7+10-limit1*7,x*8+30+11*8), font_size=8))
                 else:
#                     print(f".", end="" )
                     dwg.add(dwg.text(".", insert=(10+y*7+10-limit1*7,x*8+30+11*8), font_size=8))
             elif (y > limit2) and (y<= limit3):
                 if aa > sameness:
#                     print(f"{max_aa}", end="" )
                     dwg.add(dwg.text(max_aa, insert=(10+y*7+10-limit2*7,x*8+30+11*8*2), font_size=8))
                 else:
#                     print(f".", end="" )
                     dwg.add(dwg.text(".", insert=(10+y*7+10-limit2*7,x*8+30+11*8*2), font_size=8))
             elif (y > limit3):
                 if aa > sameness:
#                     print(f"{max_aa}", end="" )
                     dwg.add(dwg.text(max_aa, insert=(10+y*7+10-limit3*7,x*8+30+11*8*3), font_size=8))
                 else:
#                     print(f".", end="" )
                     dwg.add(dwg.text(".", insert=(10+y*7+10-limit3*7,x*8+30+11*8*3), font_size=8))
             else:
                 if aa > sameness:
#                     print(f"{max_aa}", end="" )
                     dwg.add(dwg.text(max_aa, insert=(10+y*7+10,x*8+30), font_size=8))
                 else:
#                     print(f".", end="" )
                     dwg.add(dwg.text(".", insert=(10+y*7+10,x*8+30), font_size=8))
#         print("")
         x = 0
         y += 1

dwg.add(dwg.rect(insert=(10+30,25), size=(182, 2), fill="orange"))
dwg.add(dwg.text("TM 1", insert=(10+120, 24), font_size=8))

dwg.add(dwg.rect(insert=(10+240,25), size=(165, 2), fill="orange"))
dwg.add(dwg.text("TM 2", insert=(10+320, 24), font_size=8))

dwg.add(dwg.rect(insert=(10+435,25), size=(190, 2), fill="pink"))
dwg.add(dwg.text("TM 3", insert=(10+530, 24), font_size=8))

dwg.add(dwg.rect(insert=(10+29,25+11*8), size=(98, 2), fill="pink"))
dwg.add(dwg.text("TM 4a", insert=(10+75, 24+11*8), font_size=8))

dwg.add(dwg.rect(insert=(10+157,25+11*8), size=(78, 2), fill="pink"))
dwg.add(dwg.text("TM 4b", insert=(10+190, 24+11*8), font_size=8))

dwg.add(dwg.rect(insert=(10+260,25+11*8), size=(333, 2), fill="pink"))
dwg.add(dwg.text("TM 5", insert=(10+420, 24+11*8), font_size=8))

dwg.add(dwg.rect(insert=(10+29,25+11*8*2), size=(220, 2), fill="orange"))
dwg.add(dwg.text("TM 6a", insert=(10+150, 24+11*8*2), font_size=8))

dwg.add(dwg.rect(insert=(10+275,25+11*8*2), size=(120, 2), fill="orange"))
dwg.add(dwg.text("TM 6b", insert=(10+320, 24+11*8*2), font_size=8))

dwg.add(dwg.rect(insert=(10+425,25+11*8*2), size=(180, 2), fill="orange"))
dwg.add(dwg.text("TM 7", insert=(10+510, 24+11*8*2), font_size=8))

dwg.add(dwg.rect(insert=(10+29,25+11*8*3), size=(168, 2), fill="pink"))
dwg.add(dwg.text("TM 8", insert=(10+110, 24+11*8*3), font_size=8))

dwg.add(dwg.rect(insert=(10+227,25+11*8*3), size=(70, 2), fill="pink"))
dwg.add(dwg.text("TM 9a", insert=(10+260, 24+11*8*3), font_size=8))

dwg.add(dwg.rect(insert=(10+322,25+11*8*3), size=(58, 2), fill="pink"))
dwg.add(dwg.text("TM 9b", insert=(10+345, 24+11*8*3), font_size=8))

dwg.add(dwg.rect(insert=(10+407,25+11*8*3), size=(190, 2), fill="pink"))
dwg.add(dwg.text("TM 10", insert=(10+500, 24+11*8*3), font_size=8))

dwg.add(dwg.text("PIN1", insert=(2, 30+8), font_size=8))
dwg.add(dwg.text("PIN2", insert=(2, 30+8*2), font_size=8))
dwg.add(dwg.text("PIN3", insert=(2, 30+8*3), font_size=8))
dwg.add(dwg.text("PIN4", insert=(2, 30+8*4), font_size=8))
dwg.add(dwg.text("PIN5", insert=(2, 30+8*5), font_size=8))
dwg.add(dwg.text("PIN6", insert=(2, 30+8*6), font_size=8))
dwg.add(dwg.text("PIN7", insert=(2, 30+8*7), font_size=8))
dwg.add(dwg.text("PIN8", insert=(2, 30+8*8), font_size=8))

dwg.add(dwg.text("PIN1", insert=(2, 30+8+11*8*1), font_size=8))
dwg.add(dwg.text("PIN2", insert=(2, 30+8*2+11*8*1), font_size=8))
dwg.add(dwg.text("PIN3", insert=(2, 30+8*3+11*8*1), font_size=8))
dwg.add(dwg.text("PIN4", insert=(2, 30+8*4+11*8*1), font_size=8))
dwg.add(dwg.text("PIN5", insert=(2, 30+8*5+11*8*1), font_size=8))
dwg.add(dwg.text("PIN6", insert=(2, 30+8*6+11*8*1), font_size=8))
dwg.add(dwg.text("PIN7", insert=(2, 30+8*7+11*8*1), font_size=8))
dwg.add(dwg.text("PIN8", insert=(2, 30+8*8+11*8*1), font_size=8))

dwg.add(dwg.text("PIN1", insert=(2, 30+8+11*8*2), font_size=8))
dwg.add(dwg.text("PIN2", insert=(2, 30+8*2+11*8*2), font_size=8))
dwg.add(dwg.text("PIN3", insert=(2, 30+8*3+11*8*2), font_size=8))
dwg.add(dwg.text("PIN4", insert=(2, 30+8*4+11*8*2), font_size=8))
dwg.add(dwg.text("PIN5", insert=(2, 30+8*5+11*8*2), font_size=8))
dwg.add(dwg.text("PIN6", insert=(2, 30+8*6+11*8*2), font_size=8))
dwg.add(dwg.text("PIN7", insert=(2, 30+8*7+11*8*2), font_size=8))
dwg.add(dwg.text("PIN8", insert=(2, 30+8*8+11*8*2), font_size=8))

dwg.add(dwg.text("PIN1", insert=(2, 30+8+11*8*3), font_size=8))
dwg.add(dwg.text("PIN2", insert=(2, 30+8*2+11*8*3), font_size=8))
dwg.add(dwg.text("PIN3", insert=(2, 30+8*3+11*8*3), font_size=8))
dwg.add(dwg.text("PIN4", insert=(2, 30+8*4+11*8*3), font_size=8))
dwg.add(dwg.text("PIN5", insert=(2, 30+8*5+11*8*3), font_size=8))
dwg.add(dwg.text("PIN6", insert=(2, 30+8*6+11*8*3), font_size=8))
dwg.add(dwg.text("PIN7", insert=(2, 30+8*7+11*8*3), font_size=8))
dwg.add(dwg.text("PIN8", insert=(2, 30+8*8+11*8*3), font_size=8))
#dwg.add(dwg.circle(center=(29+5*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+9*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+12*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+13*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+16*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+17*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+20*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+37*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+40*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+43*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
#dwg.add(dwg.circle(center=(29+44*7, 30+8*8+4), r=2, fill="#00AA00", stroke="#00AA00"))
xx = 46; yy = 0; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
xx = 49; yy = 0; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
xx = 19+1; yy = 1; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
xx = 43+1; yy = 1; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
xx = 46+1; yy = 1; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
xx = 51+1; yy = 1; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
#xx = 46; yy = 0; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="FFAA1D", stroke="#FFAA1D"))
#xx = 46; yy = 0; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 56; yy = 1; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#FFAA1D", stroke="#FFAA1D"))
#xx = 59; yy = 1; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#FFAA1D", stroke="#FFAA1D"))
#xx = 66; yy = 1; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#FFAA1D", stroke="#FFAA1D"))
#xx = 67; yy = 1; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#FFAA1D", stroke="#FFAA1D"))
#dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 56; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 57; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 61; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 64; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 65; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 68; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 69; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
xx = 70+1; yy = 2; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
#xx = 72; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
#xx = 75; yy = 2; dwg.add(dwg.circle(center=(29+xx*7, 30+8*8+11*8*yy+4), r=2, fill="#00AA00", stroke="#00AA00"))
xx = 42+1; yy = 3; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))
xx = 43+1; yy = 3; dwg.add(dwg.polygon(points=[(29+xx*7+2, 30+8*8+11*8*yy+4+2), (29+xx*7-2, 30+8*8+11*8*yy+4+2), (29+xx*7, 30+8*8+11*8*yy+4-1)], fill="#0000AA", stroke="#0000AA"))



dwg.save()
