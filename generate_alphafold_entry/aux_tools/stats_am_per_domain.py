#!/usr/bin/env python3
# tool which generate for given alphafold family
# its diagrams


import sys
import argparse
import time
import os
import subprocess
import shutil
import json

import utils_common
import svgwrite
sys.path.append("..")
import utils_generate_svg


def sse_avg_am(am_fn, interval):
    res_sum = 0
    am_avg = {}
    am_sum = 0
    am_count = 0
    res = 0
    with open(f"{am_fn}", "r") as fd:
         line = fd.readline()
         while line:
            sline = line.split("\t")
            n = int(sline[1][1:-1])
#            print(n, interval[0], interval[1])
            #if  n < interval[0] and n >= interval[1] and not n in res_list:
#            print("ccc", n, res)
            if n == res:
                am_count += 1
                am_sum += float(sline[2])
#                print("bbbb ", sline[2], am_fn, n)
            elif am_count > 0:
                am_avg[res] = am_sum / am_count
                res = n
                am_sum = float(sline[2])
                am_count = 1
#                print("nove", n, res)
            else:
                res = n
                am_sum = float(sline[2])
                am_count = 1
#                print("start")

            line = fd.readline()
    print("interval", interval, end=": ")
    for n in am_avg.keys():
        if  n <= interval[0] and n >= interval[1]:
            res_sum += am_avg[n]
#            print(f"{n} - {am_avg[n]}", end=", ")

    print(f"sum - {int(res_sum*100/(interval[0]-interval[1]+1))}") 
    return (int(res_sum*100/(interval[0]-interval[1]+1)))

def cl1(start, sset, tab, max_dest):
    sses_set = set([start])
    sses_all = sset
    while True:
        sses_not_set = sses_all - sses_set
        is_chosen = False
        mmin = max_dest + 1
        ch = None
        for sse in sses_not_set:
            dist = 0
#            for h in sses_set:
#                dist += tab[sse][h]
            for h in sses_set:
                if dist < tab[sse][h]:
                    dist = tab[sse][h]
#           print("dist je", dist, len(sses_set), h, sse,  tab[sse])
            if dist < max_dest and mmin > dist:
#            if dist/len(sses_set) < max_dest and mmin > dist/len(sses_set): #+max_dest*0.05*len(sses_set):
                is_chosen = True
                mmin = dist/len(sses_set)
                ch = sse
                #sses_set.add(sse)
#                break

        if not is_chosen:
            break
        else:
            sses_set.add(ch)
    return sses_set


def draw_am_diagrams(lf, prot, chain, domain, family, directory):
    svg_args = argparse.Namespace(
        add_descriptions="chain",
        layouts=[lf],
        ligands=False,
        opac='1',
        debug=False,
        irotation=None,
        itemplate=None,
        iligands=None,
        orotation=f"{directory}/rotation.json",
        otemplate=f"{directory}/{family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{new_dir}/image-{prot}_{chain}:{domain}_am.svg",
        channels = False,
        am = True,
        data=None,
    )
    r = utils_generate_svg.compute(svg_args)

def draw_multi_am_diagrams(ld, family, directory):
    svg_args = argparse.Namespace(
        add_descriptions="multi, template",
        layouts=ld,
        ligands=False,
        opac='0',
        debug=False,
        irotation=f"{new_dir}/rotation.json", #None,
        itemplate=None,
        iligands=None,
        only="sses",
        orotation=None, #
        otemplate=f"{directory}/{family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{directory}/{family}_am.svg",
        channels=False,
        am=True,
        data=None,
        fill=False,
    )
    r = utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{directory}")
    svg_args = argparse.Namespace(
        add_descriptions="multi, template",
        layouts=ld,
        ligands=False,
        opac='0',
        only="sses",
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        itemplate=f"{directory}/{family}.json",
        iligands=None,
        orotation=None, #f"{directory}/rotation.json",
        otemplate=None,
        oligands=f"{directory}/ligands.json",
        output=f"{directory}/{family}_am2.svg",
        channels = False,
        am = True,
        data=None,
        fill=False,
    )
    r = utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{directory}")


# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--layout_dir', required=True, action="store",
                    help="name of the directory, which contains domain layout")
parser.add_argument('--name', required=True, action="store",
                    help="name of the dataset")
parser.add_argument('--am_dir', required=True, action="store",
                    help="name of the dairectory which contains am data")

args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
if args.layout_dir[-1] =="/":
    name = args.layout_dir[:-1]
else:
    name = args.layout_dir

dom_data = {}
fn_order = []
sses_order = []
for filename in os.listdir(args.layout_dir):
     if filename.startswith("layout-"):
         fnsplit = filename.split("-")
         protein = fnsplit[1]
#         print("file", filename)
         if not os.path.isfile(f"{args.layout_dir}/new_layout/{filename}"):
             continue
         with open(f"{args.layout_dir}/{filename}", "r") as fd:
             data = json.load(fd)
#         print("filename", filename, data["helices_residues"])
         fn_order.append(filename)
         dom_data[filename] = {}
         for sse in data["helices_residues"]:
             if int(data["h_type"][sse]) == 1:
                 if sse not in sses_order:
                     sses_order.append(sse)
                 am_fn = f"{args.am_dir}/{protein}.am"
                 if os.path.isfile(am_fn):
                     dom_data[filename][sse] = sse_avg_am(am_fn, data['helices_residues'][sse])
#                     print("aaaa", filename, sse,  dom_data[filename][sse], data['helices_residues'][sse])
#         print(filename, dom_data[filename])
#         sys.exit()

val = 15*len(sses_order)
new_dir = f"am-{args.name}-{val}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)
print("ssss new_dir", new_dir)

diff_matrix = {}
ssum = {}
for f in  fn_order:
    fnsplit = f.split("-")
    fp = fnsplit[1]
#    print(fp, end=" ")
#print()

for f in  fn_order:
    fnsplit = f.split("-")
    fp = fnsplit[1]
    diff_matrix[f] = {}
#   print("bla ", f, end=" ")
    
    for g in fn_order:
        gnsplit = g.split("-")
        gp = gnsplit[1]
#       print("dsum je to 0")
        dsum = 0
        s = []
        for sse in (dom_data[f] | dom_data[g]): 
            if sse in dom_data[f] and sse in dom_data[g]:
                s.append(dom_data[f][sse])
                dsum += abs(dom_data[f][sse] - dom_data[g][sse])
#               print("dsum po", dsum, dom_data[f][sse], dom_data[g][sse], f, g
            elif sse not in dom_data[f] and sse not in dom_data[g]:
                dsum += 0
            else:
                dsum += 50
#               print("dsum po", 50)
        ssum[f] = s
        diff_matrix[f][g] = dsum 
#       print(dsum, end=" ")
#   print()

rest = set(fn_order)
s = []
print("--------------------------")

c = 0
while len(rest)>0:
    for e in rest:
        break
    p = cl1(e, rest, diff_matrix, val)
    s.append(p)
    rest = rest - p 
    sub_dir = f"{new_dir}/{c}"
    os.makedirs(sub_dir, exist_ok=True)
    cou = 0
    shutil.copyfile(f"{args.layout_dir}/rotation.json", f"{new_dir}/rotation.json")
    for i in p:
#        print("ssss ", i)
        if os.path.isfile(f"{args.layout_dir}/new_layout/{i}"):
            shutil.copyfile(f"{args.layout_dir}/new_layout/{i}", f"{sub_dir}/{i}")
#            print("i je", i, "---")
            si = i.split("-")
#            print("co", f"https://alphafold.ebi.ac.uk/files/AF-{si[1]}-F1-model_v4cif")
#            print("kam", sub_dir, f"{si[1]}.cif")
            utils_common.wget_file(f"https://alphafold.ebi.ac.uk/files/AF-{si[1]}-F1-model_v4.cif", sub_dir, f"{sub_dir}/{si[1]}.cif")
            cou += 1
#    print("counter", cou)
    if cou > 0:
        draw_multi_am_diagrams([sub_dir], args.name, sub_dir)

    for f in p:
        fnsplit = f.split("-")
        fp = fnsplit[1]
#        print(fp, end=" ")
        for g in p:
            gnsplit = g.split("-")
            gp = gnsplit[1]
#        print("|", ssum[f])
#    print(len(p))
    c += 1

print("ssss", len(s), end=": ")
for p in s:
    print(len(p), end=", ")

#    shutil.copyfile(f"{args.layout_dir}/{p}", f"{}")
print()

print("tabulka")
c = 0
for p in s:
    c += 1
    for f in p:
        print(f"{c} - {f} :", end=" ")
        for sse in sses_order:
            if sse in dom_data[f]:
                print(f"{dom_data[f][sse]}", end=" ")
            else:
                print("-", end=" ")
        print()


