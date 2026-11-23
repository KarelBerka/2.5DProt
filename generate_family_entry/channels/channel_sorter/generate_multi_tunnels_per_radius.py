#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import argparse
import json
import os
import sys
import utils_tunnel_pdb
import utils_read_cif
import time
import svgwrite
import math


def has_chain_borders(protein_name, chain_id, main_dir, tunnel_type_name):
    # find whether output file does not contain log about not completed borders
    # - like: TODO - exapmle
    output_filename = f"{main_dir}/output"
    tunnel_type_str = f" {tunnel_type_name}:"
    protein_chain_str = f"{protein_name}{chain_id}"

    # we go through output file and try to find record about problems with uncomplete borders
    with open(output_filename, "r") as fd:
        line = fd.readline()
        while line:
            if (tunnel_type_str in line) and (protein_chain_str in line) and ("missing borders" in line):
                # record is there -> borders are uncomplete
#                print("line", line)
                return False
            line = fd.readline()

    # record was not there -> no problem
    return True


def find_max_tunnel(tunnel_list):
    width = 0
    file_name = ""
    for fn in tunnel_list:
        if tunnel_list[fn] > width:
            width = tunnel_list[fn]
            file_name = fn
    return file_name, width


def on_list_or_higher(number, n_list):
    counter = number
    while counter not in n_list:
        counter += 1
    return counter


def on_list_or_lower(number, n_list):
    counter = number
    while counter not in n_list:
        counter -= 1
    return counter


def print_tunnel_radiuses(file_name):
    # create tunnel radius graph
    points = utils_tunnel_pdb.get_all_tunnel_points(file_name)
    i = 0
    end = points[-1]["dist"]
    last = 0
    while i < end:
        # for each interval [0,1], [1,2], .. we have to find radius
        found = False
        for r in points:
            if (r["dist"] > i) and (r["dist"] < i+1):
                last = r["radius"]
                found = True
                break
        print(last, end=" ")
        i += 1
    print()


def dist2d(a, b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


def print_tunnel_radiuses_to_list(file_name, cif_dir, pc_name):
    # create tunnel radius graph
    if cif_dir == None:
        cif_file = None
    else:
        cif_file = f"{cif_dir}/{pc_name[:-2]}_rot.cif"
    points = utils_tunnel_pdb.get_all_tunnel_points(file_name)
    ret_list = []
    i = 0
    end = points[-1]["dist"]

    start = [points[0]['x'], points[0]['y'], points[0]['z']]
    last = 0
    if cif_file:
        # we have the ligands positions so test them
        hetatm_dict = utils_read_cif.find_all_ligands_coords(cif_file, "ALL", ["HOH"])
    else:
        hetatm_dict = {}
    while i < end:
        # for each interval [0,1], [1,2], .. we have to find radius
        found = False
        for r in points:
            r_dist = (dist2d(start,[r['x'], r['y'], r['z']]))
            if (r_dist > i) and (r_dist < i+1):
                last = r["radius"]
                found = True
                break
        if found == True:
            # we have to found out, weather is the ligand in tunnel and set it to the list
            if hetatm_dict == {}:
                # no info about ligands -> there is no one
                ret_list.append((last, False))
            else:
                is_near = False
                for ha in hetatm_dict:
                    for coord in hetatm_dict[ha]["coords"]:
                        if dist2d([r['x'], r['y'], r['z']], coord) < last:
                            is_near = True
                            break
                if is_near == True:
                    # ligand is in the tunnel
                    ret_list.append((last, True))
                else:
                    # ligand is not in the tunnel
                    ret_list.append((last, False))
        i += 1
    return ret_list


def compute_coords(i_p, x, y):
   dwg_x = i_p['borders']+(i_p['width']-2*i_p['borders'])/i_p['xmax']*x
   dwg_y = i_p['height']-i_p['borders']-(i_p['height']-2*i_p['borders'])/i_p['ymax']*y
   return [dwg_x, dwg_y]


def draw_axes(dwg, i_p):
    # draw coordinates system
    dwg.add(dwg.line(start=compute_coords(i_p, 0, i_p['ymax']), end=compute_coords(i_p, 0, 0), stroke_width=3, stroke="black"))
    dwg.add(dwg.line(start=compute_coords(i_p, 0, 0), end=compute_coords(i_p, i_p['xmax'], 0), stroke_width=3, stroke="black"))
    # draw line for radius 1,2,3
    dwg.add(dwg.line(start=compute_coords(i_p, 0, 1), end=compute_coords(i_p, i_p['xmax'], 1), stroke_width=3, stroke="#BBBBBB"))
    dwg.add(dwg.line(start=compute_coords(i_p, 0, 2), end=compute_coords(i_p, i_p['xmax'], 2), stroke_width=3, stroke="#BBBBBB"))
    dwg.add(dwg.line(start=compute_coords(i_p, 0, 3), end=compute_coords(i_p, i_p['xmax'], 3), stroke_width=3, stroke="#BBBBBB"))

    xx = 0
    while xx < i_p['xmax']:
        dwg.add(dwg.line(start=compute_coords(i_p, xx, 0.1), end=compute_coords(i_p, xx, -0.1), stroke_width=3, stroke="black"))
        if xx < 10:
            dwg.add(dwg.text(xx, insert=compute_coords(i_p, xx-0.2, -0.8), font_size=20))
        else:
            dwg.add(dwg.text(xx, insert=compute_coords(i_p, xx-0.6, -0.8), font_size=20))
        xx += 2
    yy = 0
    dwg.add(dwg.text("disance from tunnel start (in A)", insert=compute_coords(i_p, xx-9, 0.3), font_size=15))
    while yy < i_p['ymax']:
        dwg.add(dwg.line(start=compute_coords(i_p, 0.3, yy), end=compute_coords(i_p, -0.3, yy), stroke_width=3, stroke="black"))
        dwg.add(dwg.text(yy, insert=compute_coords(i_p, -1, yy-0.3), font_size=20))
        yy += 1
    dwg.add(dwg.text("radius (in A)", insert=compute_coords(i_p, 0.5, yy-0.3), font_size=15))


def set_color_for_width_1_20_1070_10(width):
    if (width == 0):
        color = "#FF0000"
    elif (width <= 0.8):
        color = "#880088"
    elif (width <= 1.0):
        color = "#0000FF"
    elif (width <= 1.2):
        color = "#008888"
    elif (width <= 1.4):
        color = "#00FF00"
    elif (width <= 1.6):
        color = "#888800"
    else:
        color = "#FF0000"
    return color


def set_color_for_width(width):
    if (width == 0):
        color = "#FF0000"
    elif (width <= 1):
        color = "#FF0088"
    elif (width <= 1.5):
        color = "#FF00FF"
    elif (width <= 2):
        color = "#8800FF"
    elif (width <= 2.5):
        color = "#0000FF"
    else:
        color = "#000000"
    return color


def draw_tunnel_to_graph(dwg, i_p, t_r_list, tunnel_point, color, color_spectrum="1.10.630.10", name="."):
    # for radius list draw a line to result graph using colour based on bottleneck width
    width = tunnel_point[0]
    ligand_presence = tunnel_point[1]
    if color_spectrum == "1.10.630.10":
        sw = 0.5
    else:
        sw = 3

    if color == None:
        if color_spectrum == "1.20.1070.10":
            color = set_color_for_width_1_20_1070_10(width)
            sw = 1
        elif color_spectrum == "1.10.630.10":
            color = set_color_for_width(width)
            sw = 0.5
        elif color_spectrum == "PIN":
            color = set_color_for_width_1_20_1070_10(width)
            sw = 1
        elif args.color_spectrum == "P11169" or args.color_spectrum == "SLC":
            color = set_color_for_width_1_20_1070_10(width)
            print("i_p", t_r_list)
            if t_r_list[10][0]>1.8:
                color="#FF00FF"
            elif t_r_list[-3][0]>2:
                color = "#00FFFF"
            else:
                color = "#0000FF"
            sw = 1
        else:
            color = set_color_for_width(width)
            sw = 1

    else:
        sw = 3
    old_radius = None
    counter = 0
    for radius in t_r_list:
        if old_radius != None:
            ligand_presence = radius[1]
            if ligand_presence:
                color_act = "black"
                sw_act = sw*4
            else:
                color_act = color
                sw_act = sw
            print(f"{radius[0]}", end=" ")
            dwg.add(dwg.line(start=compute_coords(i_p, counter-1, old_radius[0]), end=compute_coords(i_p, counter, radius[0]), stroke_width=sw_act, stroke=color_act))
        old_radius = radius
        counter += 1
    dwg.add(dwg.text(name, insert=compute_coords(i_p, counter-1, old_radius[0]), font_size=10))
    print()


# read arguments
parser = argparse.ArgumentParser(
    description='''print pymol command to show tunnelborders \n
                 eg.:\n ./generate_multi_tunnels_per_radius.py --tunnel_name sss --dir 1.20.1070.10/generated-tunnels-list_1.20.1070.10.small-2026-01-06T11_39_58/ --tunnel s --ommit_start 5 --ommit_end 5 --color_spectrum 1.20.1070.10 \n
                 ''')
parser.add_argument('--tunnel', required=True, action="store",
                    help="id of tunnel (2a,3,...")
parser.add_argument('--dir', required=True, action="store",
                    help="directory containing 'annotation' dir with annotation files, '2e' dir with tunnel files, 'cifs' dir with cif files")
parser.add_argument('--tunnel_uniq', required=False, action="store_true",
                    help="show only one tunnel of given tunnel type for each domain")
parser.add_argument('--tunnel_name', required=False, action="store",
                    help="")
parser.add_argument('--add_domains_name', required=False, action="store_true", default=False, 
                    help="")
parser.add_argument('--ommit_start', required=False, default=0, action="store",
                    help="if you went to ommit a start part of tunnellt write here how many A")
parser.add_argument('--ommit_end', required=False, default=1, action="store",
                    help="if you went to ommit a start part of tunnellt write here how many A")
parser.add_argument('--color_spectrum', required=False, default="1.10.630.10", action="store",
                    help="if you went to ommit a start part of tunnellt write here how many A")
parser.add_argument('--ligands', required=False, default=False, action="store_true",
                    help="parts of tunnels wgich contains ligand will be black - needs cifs subdirectory of dir with _rot.cofs files")

args = parser.parse_args()

channels_clusters = {}
main_dir = args.dir
tunnels_dir = f"{main_dir}/{args.tunnel}/"

if args.ligands:
    cif_dir = f"{main_dir}/cifs"
else:
    cif_dir = None

tunnel_list = {}
for file_name in os.listdir(tunnels_dir):
    # sort tunnels based on domain to which they belongs
    if file_name.startswith(f"tunnel_") or file_name.startswith(f"path_"):
        tunnel_file = f"{tunnels_dir}/{file_name}"
        sfile_name = file_name.split("_")
        if len(sfile_name[1])<7:
             continue
        tunnel_pc_name = sfile_name[1]
        width = utils_tunnel_pdb.get_tunnel_min_radius(tunnel_file, 5)
        if tunnel_pc_name not in tunnel_list:
            tunnel_list[tunnel_pc_name] = {}
        tunnel_list[tunnel_pc_name][file_name] = width
print("tunnel_list", tunnel_list)

# create the file for output graphics
if args.color_spectrum == "1.20.1070.10":
    width = 1000
    xmax = 25
elif args.color_spectrum == "PIN":
    width = 1000
    xmax = 50
elif args.color_spectrum == "1.10.630.10":
    if args.tunnel == "2b":
        xmax = 26
    else:
        xmax = 28
    width = 1000
elif args.color_spectrum == "P11169" or args.color_spectrum == "SLC":
    xmax = 50
    width = 1000
else:
    width = 1000
    xmax = 45

image_param = {'width':width,
              'height':200,
              'borders':30,
              'xmax':xmax,     # tunnel can have start <-> end dist 30 A
              'ymax':5,      # tunnel can be upto 5A width
             }

dirname = os.path.basename(os.path.dirname(main_dir))
if dirname == "":
    dirname = main_dir
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
output_file_name = f"radius-{dirname}-{args.tunnel}-{timestamp}.svg"
dwg = svgwrite.Drawing(output_file_name, size=(image_param['width'], image_param['height']))
dwg.add(dwg.rect(insert=(1, 1), size=(image_param['width'], image_param['height']), fill="white"))
draw_axes(dwg, image_param)

# draw a line for each tunnel
if args.tunnel_uniq:
    # show only the maximal tunnel for each domain
    rad = {}
    # at first find only maximal tunnels - ad put it to rad dict
    for pc_name in tunnel_list:
        file_name, width = find_max_tunnel(tunnel_list[pc_name])
        if file_name == "":
            continue
        tunnel_radiuses_list = print_tunnel_radiuses_to_list(f"{tunnels_dir}/{file_name}", cif_dir, pc_name)
        rad[file_name] = min(tunnel_radiuses_list[int(args.ommit_start):-int(args.ommit_end)])
    # do rad2 dict - tunnels are sorted by radiuses (the widest have to be draw at the end (not to be overwritten by the narrow)
    rad2=dict(sorted(rad.items(), key=lambda item: item[1]))
    # draw all tunnels from rad2
    for pc_name in rad2:
        if args.tunnel_name and args.tunnel_name in pc_name:
            color = "black"
        else:
            color = None
        spc_name = pc_name.split("_")
        tunnel_radiuses_list = print_tunnel_radiuses_to_list(f"{tunnels_dir}/{pc_name}", cif_dir, spc_name[1])
        print("pc_name", pc_name, end = ": ")
        if args.add_domains_name:
            draw_tunnel_to_graph(dwg, image_param, tunnel_radiuses_list, rad2[pc_name], color, args.color_spectrum, pc_name)
        else:
            draw_tunnel_to_graph(dwg, image_param, tunnel_radiuses_list, rad2[pc_name], color, args.color_spectrum, "")
else:
    rad = {}
    # compute all tunnels and put it to rad dict
    for pc_name in tunnel_list:
        for file_name in tunnel_list[pc_name]:
            tunnel_radiuses_list = print_tunnel_radiuses_to_list(f"{tunnels_dir}/{file_name}", cif_dir, pc_name)
            if len(tunnel_radiuses_list)> int(args.ommit_start) + int(args.ommit_end):
                rad[file_name] = min(tunnel_radiuses_list[int(args.ommit_start):-int(args.ommit_end)])
    # do rad2 dict - tunnels are sorted by radiuses (the widest have to be draw at the end (not to be overwritten by the narrow)
    rad2=dict(sorted(rad.items(), key=lambda item: item[1]))
    # draw all tunnels from rad2
    for file_name in rad2:
        if args.tunnel_name and args.tunnel_name in file_name:
            color = "black"
        else:
            color = None
        sfile_name = file_name.split("_")
        tunnel_radiuses_list = print_tunnel_radiuses_to_list(f"{tunnels_dir}/{file_name}", cif_dir, sfile_name[1])
        name=file_name.split("_")
        if args.add_domains_name:
            draw_tunnel_to_graph(dwg, image_param, tunnel_radiuses_list, rad2[file_name], color, args.color_spectrum, name[1][:-3])
        else:
            draw_tunnel_to_graph(dwg, image_param, tunnel_radiuses_list, rad2[file_name], color, args.color_spectrum, "")

if args.tunnel_name:
  # draw chosen proteins once omore to have their graph on the top
  if args.tunnel_uniq:
    for pc_name in rad2:
        # show only the maximal tunnel 
        if args.tunnel_name and args.tunnel_name in pc_name:
            # only if it is set in tunnel_name
            tunnel_radiuses_list = print_tunnel_radiuses_to_list(f"{tunnels_dir}/{file_name}", cif_dir, pc_name)
            color = "black"
            draw_tunnel_to_graph(dwg, image_param, tunnel_radiuses_list, width, color)
  else:
    for file_name in rad2:
        if args.tunnel_name and args.tunnel_name in file_name:
            # only if it is set in tunnel_name
             color = "black"
             tunnel_radiuses_list = print_tunnel_radiuses_to_list(f"{tunnels_dir}/{file_name}", cif_dir, pc_name)
             draw_tunnel_to_graph(dwg, image_param, tunnel_radiuses_list, rad2[file_name], color)

dwg.save()
print("output_file_name", output_file_name)
