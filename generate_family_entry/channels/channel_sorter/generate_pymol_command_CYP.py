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

MODEL = ["1akd", "A", "00"]


def sses_model_mapping_1akd(model_annotation_file):
    # secondary structures mapping to names based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 34 and sse["end"] == 36:
                mapping["A'"] = sse["label"]
            if sse["start"] == 38 and sse["end"] == 46:
                mapping["A"] = sse["label"]
            if sse["start"] == 53 and sse["end"] == 56:
                mapping["1_1"] = sse["label"]
            if sse["start"] == 68 and sse["end"] == 77:
                mapping["B"] = sse["label"]
            if sse["start"] == 90 and sse["end"] == 95:
                mapping["B'"] = sse["label"]
            if sse["start"] == 105 and sse["end"] == 119:
                mapping["C"] = sse["label"]
            if sse["start"] == 121 and sse["end"] == 145:
                mapping["D"] = sse["label"]
            if sse["start"] == 150 and sse["end"] == 167:
                mapping["E"] = sse["label"]
            if sse["start"] == 171 and sse["end"] == 186:
                mapping["F"] = sse["label"]
            if sse["start"] == 193 and sse["end"] == 214:
                mapping["G"] = sse["label"]
            if sse["start"] == 219 and sse["end"] == 225:
                mapping["H"] = sse["label"]
            if sse["start"] == 235 and sse["end"] == 267:
                mapping["I"] = sse["label"]
            if sse["start"] == 268 and sse["end"] == 277:
                mapping["J"] = sse["label"]
            if sse["start"] == 278 and sse["end"] == 292:
                mapping["K"] = sse["label"]
            if sse["start"] == 298 and sse["end"] == 301:
                mapping["1_4"] = sse["label"]
            if sse["start"] == 322 and sse["end"] == 327:
                mapping["K'"] = sse["label"]
            if sse["start"] == 360 and sse["end"] == 378:
                mapping["L"] = sse["label"]
            if sse["start"] == 391 and sse["end"] == 392:
                mapping["4_1"] = sse["label"]
    return mapping


def set_color_for_width3(width):
    if (width == 0):
        color = '"grey"'
    elif (width <= 1):
        color = '"tv_green"'
    elif (width <= 1.5):
        color = '"tv_yellow"'
    elif (width <= 2):
        color = '"tv_orange"'
    elif (width <= 2.5):
        color = '"tv_red"'
    else:
        color = '"tv_blue"'
    return color


def define_colors_in_python2():
    print("set_color green1, [1,1,0]")
    print("set_color green2, [1,0.75,0]")
    print("set_color green3, [1,0.5,0]")
    print("set_color green4, [1,0.25,0]")
    print("set_color green5, [1,0,0]")
    print("set_color green6, [1,0,0.25]")
    
#    print("set_color green1, [0.16,0.2,0.16]")
#    print("set_color green2, [0.32,0.4,0.32]")
#    print("set_color green3, [0.48,0.6,0.48]")
#    print("set_color green4, [0.64,0.8,0.64]")
#    print("set_color green5, [0.8,1,0.8]")


def define_colors_in_python():
    print("set_color red1, [0.5,0,0]")
    print("set_color red2, [0.5,0,0.25]")
    print("set_color red3, [0.5,0,0.5]")
    print("set_color red4, [0.25,0,0.5]")
    print("set_color red5, [0.0,0,0.5]")
    print("set_color red6, [1,1,1]")


def set_tunnel_color_for_width(width):
    if (width == 0):
        color = '"red1"'
    elif (width <= 1):
        color = '"red2"'
    elif (width <= 1.5):
        color = '"red3"'
    elif (width <= 2):
        color = '"red4"'
    elif (width <= 2.5):
        color = '"red5"'
    else:
        color = '"red6"'
    return color


def set_border_color_for_width(width):
    if (width == 0):
        color = '"green1"'
    elif (width <= 1):
        color = '"green2"'
    elif (width <= 1.5):
        color = '"green3"'
    elif (width <= 2):
        color = '"green4"'
    elif (width <= 2.5):
        color = '"green5"'
    else:
        color = '"green6"'
    return color


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
                return False
            line = fd.readline()

    # record was not there -> no problem
    return True


def set_border_color(tunnels_dir, protein_name, chain_id, main_dir):
    # if given protein do not have all borders - set dark gray - there can be tunnel
    # but can't be detected
    tunnel_type_name = os.path.basename(os.path.dirname(tunnels_dir))
    has_all_borders = has_chain_borders(protein_name, chain_id, main_dir, tunnel_type_name)
    if not has_all_borders:
        return '"grey30"', None

    # protein has all bordest - thus we can look for the width
    # of the maximal bottleneck of of its tunnels and set the color based on it
    width = 0
    for file_name in os.listdir(tunnels_dir):
        if file_name.startswith(f"tunnel_{protein_name}{chain_id}"):
            tunnel_file = f"{tunnels_dir}/{file_name}"
            width_loc = utils_tunnel_pdb.get_tunnel_min_radius(tunnel_file, 5)
            if width_loc > width:
                width = width_loc
    return set_border_color_for_width(width), width


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
#        if counter > 5000:
#            return None
    return counter


def on_list_or_lower(number, n_list):
    counter = number
    while counter not in n_list:
        counter -= 1
#        if counter < 0
#            return None
    return counter


###########################################################################
## extract tunnels borders
###########################################################################

def extract_1_borders(mapping, sses_list, la_diff):

    if mapping["C"] in sses_list:
        start = sses_list[mapping["C"]]["start"]
    else:
        start = on_list_or_lower(123, la_diff)
    if mapping["E"] in sses_list:
        end = sses_list[mapping["E"]]["end"]
    else:
        end = on_list_or_higher(189, la_diff)
    print(f"extract {protein_name}{chain_id}_1_1, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")

    if mapping["L"] in sses_list:
        start = sses_list[mapping["L"]]["start"]
    else:
        start = on_list_or_lower(445, la_diff)
    if mapping["L"] in sses_list:
        end = sses_list[mapping["L"]]["end"]
    else:
        end = on_list_or_higher(462, la_diff)
    print(f"extract {protein_name}{chain_id}_1_2, {protein_name}{chain_id}_rot and resi {la_diff[start]-10}-{la_diff[end]}")


def extract_2a_borders(mapping, sses_list, la_diff):
    # 2a tunnel has 2 parts
    # upper part (from 1_1 to the middle of B-C loop)
    if mapping["1_1"] in sses_list:
        start = sses_list[mapping["1_1"]]["end"]
    else:
        start = on_list_or_lower(56, la_diff)
    if mapping["B"] in sses_list:
        end1 = sses_list[mapping["B"]]["end"]
    else:
        end1 = on_list_or_higher(85, la_diff)
    if mapping["C"] in sses_list:
        end2 = sses_list[mapping["C"]]["start"]
    else:
        end2 = on_list_or_higher(101, la_diff)
    print(f"extract {protein_name}{chain_id}_2a_1, {protein_name}{chain_id}_rot and resi {la_diff[start]-10}-{int((la_diff[end1]+la_diff[end2]+10)/2)}")

    # F- G loop (more the G part)
    if mapping["F"] in sses_list:
        start = sses_list[mapping["F"]]["end"]
    else:
        start = on_list_or_lower(188, la_diff)
    if mapping["G"] in sses_list:
        end = sses_list[mapping["G"]]["start"]
    else:
        end = on_list_or_higher(229, la_diff)
    print(f"extract {protein_name}{chain_id}_2a_2, {protein_name}{chain_id}_rot and resi {la_diff[start]-4}-{la_diff[end]+4}")


def extract_2b_borders(mapping, sses_list, la_diff):
    # 2b tunnel has 2 parts
    if mapping["1_1"] in sses_list:
        start = sses_list[mapping["1_1"]]["end"]
    else:
        start = on_list_or_lower(54, la_diff)
    if mapping["C"] in sses_list:
        end = sses_list[mapping["C"]]["start"]
    else:
        end = on_list_or_higher(101, la_diff)
    print(f"extract {protein_name}{chain_id}_2b_1, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")

    if mapping["1_4"] in sses_list:
        start = sses_list[mapping["1_4"]]["start"]
    else:
        start = on_list_or_lower(351, la_diff)
    if mapping["K'"] in sses_list:
        end = sses_list[mapping["K'"]]["start"]
    else:
        end = on_list_or_higher(378, la_diff)
    print(f"extract {protein_name}{chain_id}_2b_2, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")


def extract_2c_borders(mapping, sses_list, la_diff):
    # 2c tunnel has 3 parts
    if mapping["G"] in sses_list:
        start = sses_list[mapping["G"]]["start"]
    else:
        start = on_list_or_lower(222, la_diff)

    if mapping["G"] in sses_list:
        end = sses_list[mapping["G"]]["end"]
    else:
        end = on_list_or_higher(237, la_diff)
    print(f"extract {protein_name}{chain_id}_2c_1, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")

    if mapping["B"] in sses_list:
        start1 = sses_list[mapping["B"]]["end"]
    else:
        start1 = on_list_or_lower(76, la_diff)
    if mapping["C"] in sses_list:
        start2 = sses_list[mapping["C"]]["start"]
    else:
        start2 = on_list_or_lower(102, la_diff)
    start = int((start1+start2)/2)
    if mapping["C"] in sses_list:
        end = sses_list[mapping["C"]]["end"]
    else:
        end = on_list_or_higher(117, la_diff)
    print(f"extract {protein_name}{chain_id}_2c_2, {protein_name}{chain_id}_rot and resi {int((la_diff[start1]+la_diff[start2])/2)}-{la_diff[end]}")

    if mapping["I"] in sses_list:
        start = sses_list[mapping["I"]]["start"]
    else:
        start = on_list_or_lower(264, la_diff)
    if mapping["I"] in sses_list:
        end2 = sses_list[mapping["I"]]["end"]
    else:
        end2 = on_list_or_higher(302, la_diff)
    end = int((start+end2)/2)+6
    print(f"extract {protein_name}{chain_id}_2c_3, {protein_name}{chain_id}_rot and resi {la_diff[start]-6}-{int((la_diff[start]+la_diff[end2])/2)+6}")


def extract_2e_borders(mapping, sses_list, la_diff):
    if mapping["B"] in sses_list:
        start = sses_list[mapping["B"]]["start"]
    else:
        start = on_list_or_lower(67, la_diff)

    if mapping["C"] in sses_list:
        end = sses_list[mapping["C"]]["end"]
    else:
        end = on_list_or_higher(117, la_diff)
    print(f"extract {protein_name}{chain_id}_2e_1, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")


def extract_2f_borders(mapping, sses_list, la_diff):
    # 2f tunnel has 3 parts:
    # F-G loop
    if mapping["F"] in sses_list:
        start = sses_list[mapping["F"]]["start"]
    else:
        start = on_list_or_lower(178, la_diff)
    if mapping["G"] in sses_list:
        end = sses_list[mapping["G"]]["end"]
    else:
        end = on_list_or_higher(250, la_diff)
    print(f"extract {protein_name}{chain_id}_2f_1, {protein_name}{chain_id}_rot and resi {la_diff[start]+6}-{la_diff[end]-9}")

    # 1_1 surrounding
    if mapping["A"] in sses_list:
        start = sses_list[mapping["A"]]["start"]
    else:
        start = on_list_or_lower(33, la_diff)
    if mapping["1_1"] in sses_list:
        end = sses_list[mapping["1_1"]]["end"]
    else:
        end = on_list_or_higher(68, la_diff)
    print(f"extract {protein_name}{chain_id}_2f_2, {protein_name}{chain_id}_rot and resi {la_diff[start]-6}-{la_diff[end]+11}")

    # last loop from L to the end
    if mapping["L"] in sses_list:
        start = sses_list[mapping["L"]]["end"]
    else:
        start = on_list_or_lower(449, la_diff)
    end = pdb_end-6
    print(f"extract {protein_name}{chain_id}_2f_3, {protein_name}{chain_id}_rot and resi {la_diff[start]+11}-{end-3}")


def extract_3_borders(mapping, sses_list, la_diff):
    if mapping["E"] in sses_list:
        start = sses_list[mapping["E"]]["start"]
    else:
        start = on_list_or_lower(152, la_diff)
    if mapping["G"] in sses_list:
        end = sses_list[mapping["G"]]["end"]
    else:
        end = on_list_or_higher(237, la_diff)

    start = sses_list[mapping["E"]]["start"]
    print(f"extract {protein_name}{chain_id}_3_1, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")

    if mapping["I"] in sses_list:
        start = sses_list[mapping["I"]]["start"]
    else:
        start = on_list_or_lower(264, la_diff)
    if mapping["I"] in sses_list:
        end = sses_list[mapping["I"]]["end"]
    else:
        end = on_list_or_higher(302, la_diff)
    print(f"extract {protein_name}{chain_id}_3_2, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")


def extract_4_borders(mapping, sses_list, la_diff):
    if mapping["F"] in sses_list:
        start = sses_list[mapping["F"]]["end"]
    else:
        start = on_list_or_lower(206, la_diff)
    if mapping["G"] in sses_list:
        end = sses_list[mapping["G"]]["start"]
    else:
        end = on_list_or_higher(223, la_diff)
    print(f"extract {protein_name}{chain_id}_4_1, {protein_name}{chain_id}_rot and resi {la_diff[start]-6}-{la_diff[end]+6}")


def extract_5_borders(mapping, sses_list, la_diff):
    if mapping["I"] in sses_list:
        start = sses_list[mapping["I"]]["start"]
    else:
        start = on_list_or_lower(264, la_diff)
    if mapping["K"] in sses_list:
        end = sses_list[mapping["K"]]["end"]
    else:
        end = on_list_or_higher(100, la_diff)
    print(f"extract {protein_name}{chain_id}_5_1, {protein_name}{chain_id}_rot and resi {la_diff[start]+10}-{la_diff[end]+5}")
    if mapping["K'"] in sses_list:
        start = sses_list[mapping["K'"]]["start"]
    else:
        start = on_list_or_lower(100, la_diff)
    if mapping["L"] in sses_list:
        end = sses_list[mapping["L"]]["start"]
    else:
        end = on_list_or_higher(100, la_diff)
    print(f"extract {protein_name}{chain_id}_5_2, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")
    if mapping["A"] in sses_list:
        end = sses_list[mapping["A"]]["start"]
    else:
        end = on_list_or_higher(100, la_diff)
    print(f"extract {protein_name}{chain_id}_5_3, {protein_name}{chain_id}_rot and resi {la_diff[end]-5}-{la_diff[end]+5}")


def extract_s_borders(mapping, sses_list, la_diff):
    if mapping["L"] in sses_list:
        start = sses_list[mapping["L"]]["end"]
    else:
        start = on_list_or_lower(449, la_diff)
    end = pdb_end
    print(f"extract {protein_name}{chain_id}_s_1, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{end}")

    if mapping["E"] in sses_list:
        start = sses_list[mapping["E"]]["end"]
    else:
        start = on_list_or_lower(167, la_diff)
    if mapping["F"] in sses_list:
        end1 = sses_list[mapping["F"]]["end"]
    else:
        end1 = on_list_or_lower(190, la_diff)
    if mapping["G"] in sses_list:
        end2 = sses_list[mapping["G"]]["end"]
    else:
        end2 = on_list_or_higher(242, la_diff)
    print(f"extract {protein_name}{chain_id}_s_2, {protein_name}{chain_id}_rot and resi {la_diff[start]-2}-{(la_diff[end1]-4+la_diff[end2]-4)/2}")

    if mapping["I"] in sses_list:
        start = sses_list[mapping["I"]]["start"]
    else:
        start = on_list_or_lower(264, la_diff)
    if mapping["I"] in sses_list:
        end = sses_list[mapping["I"]]["end"]
    else:
        end = on_list_or_higher(302, la_diff)
    print(f"extract {protein_name}{chain_id}_s_3, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]-10}")

def extract_s2f_borders(mapping, sses_list, la_diff):
    # s2f tunnel has 3 parts:

    # F-G loop
    if mapping["E"] in sses_list:
        start = sses_list[mapping["E"]]["end"]
    else:
        start = on_list_or_lower(178, la_diff)
    if mapping["G"] in sses_list:
        end = sses_list[mapping["G"]]["end"]
    else:
        end = on_list_or_higher(250, la_diff)
    print(f"extract {protein_name}{chain_id}_s2f_1, {protein_name}{chain_id}_rot and resi {la_diff[start]+6}-{la_diff[end]-9}")

    # 1_1 surrounding
    if mapping["A"] in sses_list:
        start = sses_list[mapping["A"]]["start"]
    else:
        start = on_list_or_lower(33, la_diff)
    if mapping["1_1"] in sses_list:
        end = sses_list[mapping["1_1"]]["end"]
    else:
        end = on_list_or_higher(68, la_diff)
    print(f"extract {protein_name}{chain_id}_s2f_2, {protein_name}{chain_id}_rot and resi {la_diff[start]-6}-{la_diff[end]+11}")

    # last loop from L to the end
    if mapping["L"] in sses_list:
        start = sses_list[mapping["L"]]["end"]
    else:
        start = on_list_or_lower(449, la_diff)
    end = pdb_end-6
    print(f"extract {protein_name}{chain_id}_s2f_3, {protein_name}{chain_id}_rot and resi {la_diff[start]+11}-{end-3}")

    if mapping["I"] in sses_list:
        start = sses_list[mapping["I"]]["start"]
    else:
        start = on_list_or_lower(264, la_diff)
    if mapping["I"] in sses_list:
        end = sses_list[mapping["I"]]["end"]
    else:
        end = on_list_or_higher(302, la_diff)
    print(f"extract {protein_name}{chain_id}_s2f_4, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]-10}")



def extract_w_borders(mapping, sses_list, la_diff):

    if mapping["L"] in sses_list:
        end = sses_list[mapping["L"]]["start"]
    else:
        end = on_list_or_lower(450, la_diff)
    print(f"extract {protein_name}{chain_id}_w_1, {protein_name}{chain_id}_rot and resi {la_diff[end]-13}-{la_diff[end]}")
    if mapping["B"] in sses_list:
        start = sses_list[mapping["B"]]["start"]
    else:
        start = on_list_or_lower(67, la_diff)

    if mapping["C"] in sses_list:
        end = sses_list[mapping["C"]]["end"]
    else:
        end = on_list_or_higher(117, la_diff)
    print(f"extract {protein_name}{chain_id}_w_2, {protein_name}{chain_id}_rot and resi {la_diff[start]}-{la_diff[end]}")


# read arguments
parser = argparse.ArgumentParser(
    description='print pymol command to show tunnelborders')
parser.add_argument('--tunnel', required=True, action="store",
                    help="id of tunnel (2a,3,...")
parser.add_argument('--dir', required=True, action="store",
                    help="directory containing 'annotation' dir with annotation files, '2e' dir with tunnel files, 'cifs' dir with cif files")
parser.add_argument('--protein', required=False, action="store",
                    help="e.g. 1tqnA00 protein which output we want - otherwise we will see all of them")
parser.add_argument('--tunnel_uniq', required=False, action="store_true",
                    help="show only one tunnel of given tunnel type for each domain")
parser.add_argument('--tunnel_type', required=False, action="store",
                    help="type of tunnel output licorice|wire|spheres (default spheres)")
parser.add_argument('--chain_a', required=False, action="store_true",
                    help="show only chain A domains")
parser.add_argument('--chain_type', required=False, action="store",
                    help="type of chain output licorice|wire|cartoon (default cartoon)")
parser.add_argument('--show', required=True, action="store",
                    help="could be set to chains|tunnels|all")

args = parser.parse_args()


channels_clusters = {}
main_dir = args.dir
model_annotation_file = f"{main_dir}/annotation/1akdA00-annotated.sses.json"
tunnels_dir = f"{main_dir}/{args.tunnel}/"

mapping = sses_model_mapping_1akd(model_annotation_file)

define_colors_in_python()
define_colors_in_python2()


if (args.show != "tunnels") and (args.show != "chains") and (args.show != "all"):
    print("show argument must be set to chains|tunnels|all")

# tunnels part
if (args.show == "tunnels") or (args.show == "all"):
    tunnel_list = {}
    for file_name in os.listdir(tunnels_dir):
        # sort tunnels based on domain to which they belongs
        if file_name.startswith("tunnel_"):
            tunnel_file = f"{tunnels_dir}/{file_name}"
            sfile_name = file_name.split("_")
            if len(sfile_name[1]) < 7:
                # not the right format probably bogus file
                continue
            tunnel_pc_name = sfile_name[1]
            width = utils_tunnel_pdb.get_tunnel_min_radius(tunnel_file, 5)
            if tunnel_pc_name not in tunnel_list:
                tunnel_list[tunnel_pc_name] = {}
            tunnel_list[tunnel_pc_name][file_name] = width

    for pc_name in tunnel_list:
        # show only the maximal tunnel for each domain
        if args.tunnel_uniq:
            file_name, width = find_max_tunnel(tunnel_list[pc_name])
            color = set_tunnel_color_for_width(width)
#            print(int(width*100), "width")
            print(f"load {tunnels_dir}/{file_name}")
            print(f"cmd.color_deep({color}, '{file_name[:-4]}', 0)")
            if args.tunnel_type:
                print(f"cmd.show_as('{args.tunnel_type}'  ,'{file_name[:-4]}')")
        else:
            # show all tunnels
            for file_name in tunnel_list[pc_name]:
                color = set_tunnel_color_for_width(tunnel_list[pc_name][file_name])
                print(f"load {tunnels_dir}/{file_name}")
                print(f"cmd.color_deep({color}, '{file_name[:-4]}', 0)")
                if args.tunnel_type:
                    print(f"cmd.show_as('{args.tunnel_type}'  ,'{file_name[:-4]}')")
if args.show == "tunnels":
    # all was done tunnels are shown, no need to shou cains as well
    sys.exit()

# chains part
if args.protein:
    protein = args.protein
    pcds = [[protein[0:4], protein[4:5], protein[5:]]]
else:
    pcds = []
    for file_name in os.listdir(f"{main_dir}/annotation/"):
        if file_name.endswith("-annotated.sses.json") and file_name != "1akdA00-annotated.sses.json":
            pcd = [file_name[0:4], file_name[4:5], file_name[5:7]]
            pcds.append(pcd)

first_pcd = True
for pcd in pcds:
    protein_name = pcd[0]
    chain_id = pcd[1]
    domain = pcd[2]
    if (chain_id != "A") and (args.chain_a):
        continue

    pcd_rot_file = f"{main_dir}/cifs/{protein_name}{chain_id}_rot.cif"
    la_diff = utils_read_cif.transform_label_to_auth_asym_id(pcd_rot_file, chain_id)
    pdb_end = utils_read_cif.find_auth_asym_id_end(pcd_rot_file, chain_id)

    color, width = set_border_color(tunnels_dir, protein_name, chain_id, main_dir)
#    if  width <2:
#        continue
    print(f"load {main_dir}/cifs/{protein_name}{chain_id}_rot.cif")

    #  from tunnel pdb file in {args.tunnel} dir find the width of tunnel and set relevant color
    print(f"cmd.color_deep({color}, '{protein_name}{chain_id}_rot', 0)")
    if args.chain_type:
        print(f"cmd.show_as('{args.chain_type}'  ,'{protein_name}{chain_id}_rot')")
    sses_list = {}
    # draw parts of chains which defines {args.tunnel} tunnel
    # TODO - filtr only main sses
    with open(f"{main_dir}/annotation/{protein_name}{chain_id}{domain}-annotated.sses.json", "r") as fa:
        data = json.load(fa)
        pcd = list(data.keys())[0]
        for sse in data[pcd]["secondary_structure_elements"]:
            sses_list[sse["label"]] = {"start": sse["start"], "end": sse["end"]}

    if args.tunnel == "1":
        # 1 tunnel
        extract_1_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "2a":
        # 2a tunnel
        extract_2a_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "2b":
        # 2b tunnel
        extract_2b_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "2c":
        # 2c tunnel
        extract_2c_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "2e":
        # 2e tunnel
        extract_2e_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "2f":
        # 2f tunnel
        extract_2f_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "3":
        # 3 tunnel
        extract_3_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "4":
        # 4 tunnel
        extract_4_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "5":
        # 5 tunnel
        extract_5_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "s":
        # solvent tunnel
        extract_s_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "w":
        # water tunnel
        extract_w_borders(mapping, sses_list, la_diff)
    elif args.tunnel == "s2f":
        # water tunnel
        extract_s2f_borders(mapping, sses_list, la_diff)

    print(f"delete {protein_name}{chain_id}_rot")
    if first_pcd:
        # we want to show only one HEM (not to overload image with too much data) - so we show the HEM for the first pcd
        print(f"load {pcd_rot_file}")
        color = '"white"'
        print(f"cmd.color_deep({color}, '{protein_name}{chain_id}_rot', 0)")
        print(f"extract {protein_name}{chain_id}_all, {protein_name}{chain_id}_rot and resi {0}-{pdb_end}")
        print(f"delete {protein_name}{chain_id}_all")
        first_pcd = False

