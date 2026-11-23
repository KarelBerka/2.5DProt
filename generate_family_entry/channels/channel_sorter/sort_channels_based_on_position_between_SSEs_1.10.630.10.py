#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import argparse
import json
import os
import numpy as np
import math
import sys
import utils_tunnel_pdb
import utils_read_cif
import utils_1_10_630_10
from gemmi import cif

#MODEL = ["1akd", "A", "00"]


def line_cross_triangle(triangle, line):
    # returns True if line segment started in l crosses triangle t
    p1 = np.array([triangle[0][0], triangle[0][1], triangle[0][2]])  # point 1
    p2 = np.array([triangle[1][0], triangle[1][1], triangle[1][2]])  # point 2
    p3 = np.array([triangle[2][0], triangle[2][1], triangle[2][2]])  # point 3

    s1 = np.array([line[0][0], line[0][1], line[0][2]])
    s2 = np.array([line[1][0], line[1][1], line[1][2]])

    p21 = np.subtract(p2, p1)              # vector from p1 to p2
    p31 = np.subtract(p3, p1)              # vector from p1 to p3

    triangle_norm = np.cross(p21, p31)    # triangle norm

    N12 = np.cross((np.subtract(p2, p1)), triangle_norm)  # normal vector to side p21
    N23 = np.cross((np.subtract(p3, p2)), triangle_norm)  # normal vector to side p32
    N31 = np.cross((np.subtract(p1, p3)), triangle_norm)  # normal vector to side p13

    D = -np.dot(triangle_norm, p1)  # distance of plane N . x to point p1 (to define plane of triangle p1,p2,p3)

    p = (np.dot(triangle_norm, (s2-s1)))
    if p != 0:
        t = -(D + np.dot(triangle_norm, s1)) / (np.dot(triangle_norm, (s2-s1)))   # intersection of triangle and line
    else:
        return False

    if t < 0 or t > 1:             # it is segment not line
        return False
    else:
        R = s1 + t * (s2 - s1)     # intersection of tringle and segmet
        dist_r_n12 = np.dot((R-p1), N12)/np.linalg.norm(N12)
        dist_r_n23 = np.dot((R-p2), N23)/np.linalg.norm(N23)
        dist_r_n31 = np.dot((R-p3), N31)/np.linalg.norm(N31)
        if dist_r_n12 < 0 and dist_r_n23 < 0 and dist_r_n31 < 0:
            return True
        else:
            return False


def line_cross_polygon(sses, line):
    # test weather line segment "line" cross polygon defined by points - sses
    if len(sses) < 3:
        print(f"polygon {sses} has <3 points")
        sys.exit()
    else:
        t0 = sses[0]
        t2 = sses[1]
        counter = 2
        cross = False

        while (counter < len(sses)):
            t1 = t2
            t2 = sses[counter]
            if line_cross_triangle([t0, t1, t2], line):
                cross = True
                counter = len(sses)
            else:
                counter += 1
        return cross


def channel_cross_polygon(sses, channel):
    # test weather set of line segments cross polygon defined by points - sses
    start = None
    end = None
    fin_ret = 0
    counter = 0
    for l in channel:
        if start is None:
            start = l
            end = l
        else:
            start = end
            end = l
            ret = line_cross_polygon(sses, [start, end])
            if ret:
                fin_ret += 1
        counter += 1
    # there were no crossing
    return fin_ret


def test_line_cross_triangle():
    # line cross triangle test function
    line = [[0, 0, 0], [2, 0, 0]]
    triangle = [[1, 1, 1], [1, 1, -1], [1, -1, 0]]
    ret = line_cross_triangle(triangle, line)

    line = [[-2, 0, 0], [2, 0, 0]]
    triangle = [[1, 1, 1], [1, 1, -1], [-1, -1, 0]]
    ret = line_cross_triangle(triangle, line)

    line = [[0, 0, 0], [-2, 0, 0]]
    triangle = [[1, 1, 1], [1, 1, -1], [1, -1, 0]]
    ret = line_cross_triangle(triangle, line)


def short_ok_long_by_length(tunnel_dir, min_len, max_len, types="tunnels"):
    # filter tunnels which are too short < min_len
    # and long > max_len
    short_tunnels = []
    long_tunnels = []
    ok_tunnels = []
    for file_name in os.listdir(tunnel_dir):
        if (types == "tunnels" and file_name.startswith("tunnel_")) or \
           (types == "pores" and file_name.startswith("pore_")) or \
           (types == "all" and (file_name.startswith("tunnel_") or (file_name.startswith("pore_")))):
            length = utils_tunnel_pdb.read_tunnel_length(f"{tunnel_dir}/{file_name}")
            if length < min_len:
                short_tunnels.append(file_name)
            elif length < max_len:
                ok_tunnels.append(file_name)
            else:
                long_tunnels.append(file_name)
    return short_tunnels, ok_tunnels, long_tunnels


def dist2d(a, b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


def near_hem(tunnel_dir, tunnel_list, cif_file, chain):
    # filter channels which starts far from hem
    ok = []

    # find FE position in HEM in protein in cif file:
    auth_chain = utils_read_cif.transform_label_asym_id_to_auth_asym_id(cif_file, chain)
    hem_coords = utils_read_cif.find_fe_hem_coords_auth(cif_file, chain_id=auth_chain)
#    print("hledame hem na ", chain, auth_chain, hem_coords, cif_file)
#    if hem_coords == 0:
#        print(f"no HEM (HEC, HEB) found in {cif_file} {chain} - label")
#        return []
    # go through tunnels
    for tunnel in tunnel_list:  # ch is e.g. 'gnuplot-4wnv00B.txt_7'
        ok.append(tunnel)
    return ok

#
#def sses_model_mapping_1akd(model_annotation_file):
#    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
#    mapping = {}
#    with open(model_annotation_file, "r") as fr:
#        data = json.load(fr)
#        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
#        for sse in data[pcd]["secondary_structure_elements"]:
#            if sse["start"] == 34 and sse["end"] == 36:
#                mapping["A'"] = sse["label"]
#            if sse["start"] == 38 and sse["end"] == 46:
#                mapping["A"] = sse["label"]
#            if sse["start"] == 53 and sse["end"] == 56:
#                mapping["1_1"] = sse["label"]
#            if sse["start"] == 68 and sse["end"] == 77:
#                mapping["B"] = sse["label"]
#            if sse["start"] == 90 and sse["end"] == 95:
#                mapping["B'"] = sse["label"]
#            if sse["start"] == 105 and sse["end"] == 119:
#                mapping["C"] = sse["label"]
#            if sse["start"] == 121 and sse["end"] == 145:
#                mapping["D"] = sse["label"]
#            if sse["start"] == 150 and sse["end"] == 167:
#                mapping["E"] = sse["label"]
#            if sse["start"] == 171 and sse["end"] == 186:
#                mapping["F"] = sse["label"]
#            if sse["start"] == 193 and sse["end"] == 214:
#                mapping["G"] = sse["label"]
#            if sse["start"] == 219 and sse["end"] == 225:
#                mapping["H"] = sse["label"]
#            if sse["start"] == 235 and sse["end"] == 267:
#                mapping["I"] = sse["label"]
#            if sse["start"] == 268 and sse["end"] == 277:
#                mapping["J"] = sse["label"]
#            if sse["start"] == 278 and sse["end"] == 292:
#                mapping["K"] = sse["label"]
#            if sse["start"] == 298 and sse["end"] == 301:
#                mapping["1_4"] = sse["label"]
#            if sse["start"] == 322 and sse["end"] == 327:
#                mapping["K'"] = sse["label"]
#            if sse["start"] == 360 and sse["end"] == 378:
#                mapping["L"] = sse["label"]
#            if sse["start"] == 382 and sse["end"] == 383:
#                mapping["3_3"] = sse["label"]
#            if sse["start"] == 403 and sse["end"] == 405:
#                mapping["3_2"] = sse["label"]
##            sses_list[sse["label"]] = {"start": sse["start_vector"], "end": sse["end_vector"]}
#    return mapping


def get_sses_set(annot_fn):
    sses_list = {}
    # TODO - filtr only main sses
    with open(f"{annot_fn}", "r") as fa:
        data = json.load(fa)
        pcd = list(data.keys())[0]
        for sse in data[pcd]["secondary_structure_elements"]:
            sses_list[sse["label"]] = {"start_vector": sse["start_vector"], "end_vector": sse["end_vector"], "start": sse["start"], "end": sse["end"]}
    return sses_list


def test_is_sse_in_protein( protein_annotation_file, mapping, name):
    # read protein sses coordinates
    sses_set = get_sses_set(protein_annotation_file)
    if mapping[name] in sses_set:
        return True
    else:
        return False


def is_channels_in_cluster(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id, border_points):
    # find wheather given channel goes through set polygon

    # read rotation files - to rotate tunnels
    with open(f"{rotation_file}", "r") as fr:
        data_rot = json.load(fr)
        r = np.linalg.inv(data_rot["rotation"])
        t = data_rot["translation"]

    # read protein sses coordinates
    sses_set = get_sses_set(protein_annotation_file)
    polygon_points = []

    # choose the polygon borders from sses
    for res in border_points:
        if res[0] not in sses_set:
            return res
            break
        if res[1] == 0:
            polygon_points.append(sses_set[res[0]]["start_vector"])
        elif res[1] == 1:
            polygon_points.append(sses_set[res[0]]["end_vector"])

    #  read and transform tunnel points
    # SSEs in annotation are rotated (ubertemplate)
    # we have to rotate tunnels as well - using ubertemplate -rotation file
    tunnel_points_orig = utils_tunnel_pdb.get_tunnel_points(tunnel_file)
    tunnel_points = []
    for sv in tunnel_points_orig:
        transf_sv = [r[0][0]*sv[0]+r[1][0]*sv[1]+r[2][0]*sv[2]+t[0][0], \
                     r[0][1]*sv[0]+r[1][1]*sv[1]+r[2][1]*sv[2]+t[1][0], \
                     r[0][2]*sv[0]+r[1][2]*sv[1]+r[2][2]*sv[2]+t[2][0]]
        tunnel_points.append(transf_sv)

    # test weather tunnel cross sses polygon
#    print("tunnel_points", tunnel_points)
    res = channel_cross_polygon(polygon_points, tunnel_points)
    if res == 1:
        return 0

    if len(polygon_points) > 3:
        # try to the second possible rotation
        polygon_points2 = polygon_points[1:] + [polygon_points[0]]
        res = channel_cross_polygon(polygon_points2, tunnel_points)
        if res == 1:
            return 0

    if len(polygon_points) > 4:
        # try to the third possible rotation
        polygon_points3 = polygon_points[1:] + [polygon_points[0]]
        res = channel_cross_polygon(polygon_points3, tunnel_points)
        if res == 1:
            return 0

    return 1


def sses_res(rec, sses_set):
    # return res number of the given point ["H22", 1] or ["H36", 0]]
    if rec[0] not in sses_set:
        return -1
    else:
        if rec[1] == 0:
            return sses_set[rec[0]]["start"]
        elif rec[1] == 1:
            return sses_set[rec[0]]["end"]


def sses_coords(rec, sses_set):
    if rec[0] not in sses_set:
        return -1
    else:
        if rec[1] == 0:
            return sses_set[rec[0]]["start_vector"]
        elif rec[1] == 1:
            return sses_set[rec[0]]["end_vector"]

def avg_coords(coord_set):
    x_sum = 0
    y_sum = 0
    z_sum = 0
    for i in coord_set:
#        print("i[0]", i[0])
        x_sum += float(i[0])
        y_sum += float(i[1])
        z_sum += float(i[2])
    return [x_sum/len(coord_set), y_sum/len(coord_set), z_sum/len(coord_set)]




def find_res_coords(res, protein_cif_file, chain_id,  r_data):
    # find coords of residue res and from the given cif file
    # r_data are rotation data

    res_values = []    # avg values of res
    doc = cif.read_file(protein_cif_file)  # copy all the data from mmCIF file
    block = doc.sole_block()  # mmCIF has exactly one block
    a = dict(block.get_mmcif_category("_atom_site"))
    i = 0
    while i < len(a['id']):
        if a['group_PDB'][i] =='ATOM' and a['label_asym_id'][i] == chain_id and int(a['label_seq_id'][i]) == res:
            orig_coords = [float(a['Cartn_x'][i]), float(a['Cartn_y'][i]), float(a['Cartn_z'][i])]
            new_coords = list(np.array(np.dot(r_data["rotation"], np.array(orig_coords)) + np.array([r_data["translation"][0][0], r_data["translation"][1][0], r_data["translation"][2][0]])))
            res_values.append(new_coords)
        i += 1
    if len(res_values)> 0:
        return avg_coords(res_values)
    else:
        return -1


def find_sses_point_res(rec, sses_set, protein_cif_file, chain_id,  r_data):
    # return coords in 3D to compute line polygon
    # [["H22", 1], 0.5, ["H36", 0]] - tested
    if (len(rec) == 3) and (type(rec[2]) != int):
        ret1 = sses_res(rec[0], sses_set)
        ret2 = sses_res(rec[2], sses_set)
        if (ret1 == -1) or (ret2 == -1):
            return -1
        diff = ret2-ret1
        return  ret1 + int(diff*rec[1])
    # ["H94", 0, -13]
    elif (len(rec) == 3) and (type(rec[2]) == int):
        ret = sses_res(rec, sses_set)
        if (ret == -1):
            return -1
        else:
            return ret + rec[2]
    # or record like ["H59", 0] - which means start of SSE "H59"
    # or record like ["H59", 1] - which means end of SSE "H59" - tested
    elif len(rec) == 2:
        ret =sses_res(rec, sses_set)
        return ret


def find_sses_point_coords(rec, sses_set, protein_cif_file, chain_id, r_data):
    # return coords in 3D to compute line polygon
    # there can be record  [["H22", 1], 0.5, ["H36", 0]]
    #             - which mean point which is in the middle
    #               of chain from H22 end to H36 start

    # [["H22", 1], 0.5, ["H36", 0]] - tested
    if (len(rec) == 3)  and (type(rec[2]) != int):
        ret1 = sses_res(rec[0], sses_set)
        ret2 = sses_res(rec[2], sses_set)
        if (ret1 == -1) or (ret2 == -1):
            return -1

        diff = ret2 - ret1
        ret = ret1 + int(diff*rec[1])
        x = find_res_coords(ret, protein_cif_file, chain_id, r_data)
        # if the wanted res is not defned try to test the residues around it as well
        if x == -1:
            x = find_res_coords(ret-1, protein_cif_file, chain_id, r_data)
        if x == -1:
            x = find_res_coords(ret+1, protein_cif_file, chain_id, r_data)
        return x
    # ["H94", 0, -13]
    elif (len(rec) == 3) and (type(rec[2]) == int):
        ret = sses_res(rec, sses_set)
        return find_res_coords(ret + rec[2], protein_cif_file, chain_id, r_data)
    # [["H22", 1]] - tested
    elif len(rec) == 2:
        ret = sses_coords(rec, sses_set)
        return ret


def find_res_list(start_res, end_res, protein_cif_file, chain_id, r_data):
    # find coords of residues between start_res and end_res from the given cif file
    if end_res < start_res:
        down = True
        aux = start_res
        start_res = end_res
        end_res = aux
    else:
        down = False

    residues_set = []       # all residues points
    residues_set_half = []  # only half of residues
    res_values = []    # avg values of res
    act_res = start_res
    r = -1
    doc = cif.read_file(protein_cif_file)  # copy all the data from mmCIF file
    block = doc.sole_block()  # mmCIF has exactly one block
    a = dict(block.get_mmcif_category("_atom_site"))
    while act_res <= end_res:
        i = 0
        while i < len(a['id']):
            if a['group_PDB'][i] == 'ATOM' and a['label_asym_id'][i] == chain_id and int(a['label_seq_id'][i]) == act_res:
                orig_coords = [float(a['Cartn_x'][i]), float(a['Cartn_y'][i]), float(a['Cartn_z'][i])]
                new_coords = list(np.array(np.dot(r_data["rotation"], np.array(orig_coords)) + np.array([r_data["translation"][0][0], r_data["translation"][1][0], r_data["translation"][2][0]])))
                res_values.append(new_coords)
            i += 1
        if len(res_values) > 0:
            residues_set.append(avg_coords(res_values))

        res_values = []
        act_res += 1
    rcounter = 0
    for r in residues_set:
        rcounter += 1
        if (r == residues_set[0]) or (r == residues_set[-1]) or (rcounter % 2 == 1):
            residues_set_half.append(r)
    if down:
        residues_set.reverse()
        residues_set_half.reverse()
    return residues_set_half


def is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id, border_points, rotation=0):
    # new method
    # sses_points   - polygon limits (sses points, loop points) e.g. [['LOOP', ['H22', 1], ['H36', 0]], ['LINE', ['H36', 0], ['H22', 1]]]
    # find wheather given channel goes through set polygon

    # read rotation files - to rotate tunnels
    with open(f"{rotation_file}", "r") as fr:
        r_data = json.load(fr)
        r = np.linalg.inv(r_data["rotation"])
        t = r_data["translation"]

    # read protein sses coordinates
    sses_set = get_sses_set(protein_annotation_file)
    polygon_points = []

    for record in border_points:
        if record[0] == "LINE":
            # ["LINE", [sses_mapping["I"], 0], [sses_mapping["G"], 1]]
            # the border is line between two SSEs ends - add two points
            # first point
            start_res = find_sses_point_coords(record[1], sses_set, protein_cif_file, chain_id, r_data)
            end_res = find_sses_point_coords(record[2], sses_set, protein_cif_file, chain_id, r_data)
            if (start_res == -1):
                return record[1]
            elif (end_res == -1):
                return record[2]

            else:
                p = [start_res, end_res]
                polygon_points += p

        if record[0] == "LOOP":
            # ["LOOP", [sses_mapping["I"], 0], [sses_mapping["G"], 1]]
            # the border is curve between two SSEs ends - add two points
            # first point - find the res number
            start_res = find_sses_point_res(record[1], sses_set, protein_cif_file, chain_id, r_data)
            # second point - fin the res nimber
            end_res = find_sses_point_res(record[2], sses_set, protein_cif_file, chain_id, r_data)
            if (start_res == -1):
                return record[1]
            elif (end_res == -1):
                return record[2]
            try:
                p = find_res_list(start_res, end_res, protein_cif_file, chain_id, r_data)
            except:
                return 5
            if len(p) == 0:
                return 6

            polygon_points += p

#        sys.exit()
    if len(polygon_points) > 1:
            polygon_points.append(polygon_points[0])

    if rotation != 0:
        max_len = len(polygon_points)
        rotation_len = int(max_len * rotation)
        aux_polygon_points = polygon_points[rotation_len:]+polygon_points[:rotation_len]
        polygon_points = aux_polygon_points

    #  read and transform tunnel points
    # SSEs in annotation are rotated (ubertemplate)
    # we have to rotate tunnels as well - using ubertemplate -rotation file
    tunnel_points_orig = utils_tunnel_pdb.get_tunnel_points(tunnel_file)

    # we add last point just to get tunnel end to the surface a bit
    # TODO - tune the length of the rest
    rest = [tunnel_points_orig[-1][0]+(tunnel_points_orig[-1][0]-tunnel_points_orig[0][0])*0.2,
            tunnel_points_orig[-1][1]+(tunnel_points_orig[-1][1]-tunnel_points_orig[0][1])*0.2,
            tunnel_points_orig[-1][2]+(tunnel_points_orig[-1][2]-tunnel_points_orig[0][2])*0.2]
    tunnel_points_orig.append(rest)

    tunnel_points = []
    for sv in tunnel_points_orig:
        transf_sv = [r[0][0]*sv[0]+r[1][0]*sv[1]+r[2][0]*sv[2]+t[0][0],
                     r[0][1]*sv[0]+r[1][1]*sv[1]+r[2][1]*sv[2]+t[1][0],
                     r[0][2]*sv[0]+r[1][2]*sv[1]+r[2][2]*sv[2]+t[2][0]]
        tunnel_points.append(transf_sv)

#    print("tunnel_points", tunnel_points)
    # test weather tunnel cross sses polygon
    res = channel_cross_polygon(polygon_points, tunnel_points)
    if res == 1:
        return 0
    else:
        return 1


def create_reverse_mapping(mapping):
    rev = {}
    for i in mapping:
        rev[mapping[i]] = i
    return rev


def is_tunnel_of_known_type(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id, sses_mapping, bugs_list, part_only):

    min_radius = utils_tunnel_pdb.get_tunnel_min_radius(tunnel_file)
    tunnel_len = utils_tunnel_pdb.read_tunnel_length(tunnel_file)
    sten = utils_tunnel_pdb.get_sten_distance(tunnel_file)
    rev_m = create_reverse_mapping(sses_mapping)
    bugs = set()

    # main -------------------------------------------------------------------------------------------
    # s
    if  test_is_sse_in_protein( protein_annotation_file, sses_mapping, "F"):
        # if SSE F is present - define s - channel polygon based on F
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["E"], 0], [sses_mapping["E"], 1]],
                                 ["LINE", [sses_mapping["E"], 1], [sses_mapping["F"], 0]],
                                 ["LINE", [sses_mapping["F"], 0], [sses_mapping["F"], 1]],
                                 ["LINE", [sses_mapping["F"], 1], [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]]],
                                 ["LOOP", [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]], [sses_mapping["3_3"], 0]],
                                 ["LINE", [sses_mapping["3_3"], 0], [sses_mapping["E"], 0]],])
    else:
        # if F is not present try to substitute  it by F'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["E"], 0], [sses_mapping["E"], 1]],
                                 ["LINE", [sses_mapping["E"], 1], [sses_mapping["F'"], 0]],
                                 ["LINE", [sses_mapping["F'"], 0], [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]]],
                                 ["LOOP", [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]], [sses_mapping["3_3"], 0]],
                                 ["LINE", [sses_mapping["3_3"], 0], [sses_mapping["E"], 0]],])
    if ret == 0:
        return "s"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders s: {ret} {rev_m}")

    # 2f ./
    if  test_is_sse_in_protein( protein_annotation_file, sses_mapping, "F"):
        # if SSE F is present - define s - channel polygon based on F
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["F"], 1], [sses_mapping["1_1"], 1]],
                                 ["LINE", [sses_mapping["1_1"], 1], [sses_mapping["A"], 0]],
                                 ["LINE", [sses_mapping["A"], 0], [sses_mapping["3_2"], 1]],
                                 ["LOOP", [sses_mapping["3_2"], 1], [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]]],
                                 ["LINE", [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]], [sses_mapping["F"], 1]]])
    else:
        # if F is not present try to substitute  it by F'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["F'"], 0], [sses_mapping["1_1"], 1]],
                                 ["LINE", [sses_mapping["1_1"], 1], [sses_mapping["A"], 0]],
                                 ["LINE", [sses_mapping["A"], 0], [sses_mapping["3_2"], 1]],
                                 ["LOOP", [sses_mapping["3_2"], 1], [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]]],
                                 ["LINE", [[sses_mapping["3_2"], 1], 0.5, [sses_mapping["3_3"], 0]], [sses_mapping["F'"], 0]]])
    if ret == 0:
        return "2f"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 2f: {ret} {rev_m}")

    # 2e ./
    # surrounding of tunnel 2e has a specific shape - we have to start counting in 25% of the loop
    # otherwise we will count short tunnels which ends below 2e
    ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LOOP", [sses_mapping["B"], 1], [sses_mapping["C"], 0]],
                                 ["LINE", [sses_mapping["C"], 0], [sses_mapping["B"], 1]]], rotation=0.25)
    if ret == 0:
        return "2e"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 2e: {ret} {rev_m}")

    # 4 ./ should be before 2a - protein P08684 region 2a and 4 overlaps and result tunnels should belong to 4
    if  test_is_sse_in_protein( protein_annotation_file, sses_mapping, "F"):
        # if SSE F is present - define s - channel polygon based on F
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LOOP", [sses_mapping["F"], 1], [sses_mapping["G"], 0]],
                                 ["LINE", [sses_mapping["G"], 0], [sses_mapping["F"], 1]]])
    else:
        # if F is not present try to substitute  it by F'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LOOP", [sses_mapping["F'"], 0], [sses_mapping["G"], 0]],
                                 ["LINE", [sses_mapping["G"], 0], [sses_mapping["F'"], 0]]])
    if ret == 0:
        return "4"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 4: {ret} {rev_m}")

    # 2a ./
    if (test_is_sse_in_protein( protein_annotation_file, sses_mapping, "B") and test_is_sse_in_protein( protein_annotation_file, sses_mapping, "C")):
        # we have SSE B and C - just use it
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["1_1"], 1], [[sses_mapping["B"], 1], 0.5, [sses_mapping["C"], 0]]],
                                 ["LINE", [[sses_mapping["B"], 1], 0.5, [sses_mapping["C"], 0]], [sses_mapping["G"], 0]],
                                 ["LOOP", [sses_mapping["G"], 0], [sses_mapping["F"], 1]],
                                 ["LINE", [sses_mapping["F"], 1], [sses_mapping["1_1"], 1]]])
    else:
        # B or C or both are not present try to use B'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["1_1"], 1], [[sses_mapping["B'"], 0], 0.5, [sses_mapping["B'"], 1]]],
                                 ["LINE", [[sses_mapping["B'"], 0], 0.5, [sses_mapping["B'"], 1]], [sses_mapping["G"], 0]],
                                 ["LOOP", [sses_mapping["G"], 0], [sses_mapping["F"], 1]],
                                 ["LINE", [sses_mapping["F"], 1], [sses_mapping["1_1"], 1]]])
    if ret == 0:
        return "2a"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 2a: {ret} {rev_m}")

    # 2b ./
    if (test_is_sse_in_protein( protein_annotation_file, sses_mapping, "B") and test_is_sse_in_protein( protein_annotation_file, sses_mapping, "C")):
        # we have SSE B and C - just use it
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["1_1"], 1], [sses_mapping["1_4"], 1]],
                                 ["LINE", [sses_mapping["1_4"], 1], [sses_mapping["B"], 1]],
                                 ["LOOP", [sses_mapping["B"], 1], [[sses_mapping["B"], 1], 0.5, [sses_mapping["C"], 0]]],
                                 ["LINE", [[sses_mapping["B"], 1], 0.5, [sses_mapping["C"], 0]], [sses_mapping["1_1"], 1]]])

    else:
        # B or C or both are not present try to use B'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["1_1"], 1], [sses_mapping["1_4"], 1]],
                                 ["LINE", [sses_mapping["1_4"], 1], [sses_mapping["B'"], 1]],
                                 ["LOOP", [sses_mapping["B"], 1], [[sses_mapping["B'"], 0], 0.5, [sses_mapping["B'"], 1]]],
                                 ["LINE", [[sses_mapping["B'"], 0], 0.5, [sses_mapping["B'"], 1]], [sses_mapping["1_1"], 1]]])

    if ret == 0:
        return "2b"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 2b: {ret} {rev_m}")

    # 2c ./
    if  test_is_sse_in_protein( protein_annotation_file, sses_mapping, "B"):
        # we have SSE B and C - just use it
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["I"], 0], [sses_mapping["G"], 1]],
                                 ["LINE", [sses_mapping["G"], 1], [sses_mapping["G"], 0]],
                                 ["LINE", [sses_mapping["G"], 0], [[sses_mapping["B"], 1], 0.5, [sses_mapping["C"], 0]]],
                                 ["LOOP", [[sses_mapping["B"], 1], 0.5, [sses_mapping["C"], 0]], [sses_mapping["C"], 0]],
                                 ["LINE", [sses_mapping["C"], 0], [sses_mapping["I"], 0]]])
    else:
        # B or C or both are not present try to use B'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["I"], 0], [sses_mapping["G"], 1]],
                                 ["LINE", [sses_mapping["G"], 1], [sses_mapping["G"], 0]],
                                 ["LINE", [sses_mapping["G"], 0], [[sses_mapping["B'"], 0], 0.5, [sses_mapping["B'"], 1]]],
                                 ["LOOP", [[sses_mapping["B'"], 0], 0.5, [sses_mapping["B'"], 1]], [sses_mapping["C"], 0]],
                                 ["LINE", [sses_mapping["C"], 0], [sses_mapping["I"], 0]]])
    if ret == 0:
        return "2c"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 2c: {ret} {rev_m}")

    # small ------------------------------------------------------------------------------------------
    # 5 ./
    ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["K'"], 1], [sses_mapping["K"], 0]],
                                 ["LOOP", [sses_mapping["K"], 0], [sses_mapping["K'"], 1]]])
    if ret != 0:
        ret = is_channels_in_cluster(tunnel_file, protein_annotation_file, rotation_file,
                                 protein_cif_file, chain_id,
                                 [[sses_mapping["K"], 0], [sses_mapping["K"], 1], [sses_mapping["A"], 0],
                                  [sses_mapping["A"], 1], [sses_mapping["K'"], 0], [sses_mapping["K'"], 1]])
    if ret == 0:
        return "5"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 5: {ret} {rev_m}")

    # 3
    if test_is_sse_in_protein( protein_annotation_file, sses_mapping, "F"):
        # if SSE F is present - define s - channel polygon based on F
        ret = is_channels_in_cluster(tunnel_file, protein_annotation_file, rotation_file,
                                 protein_cif_file, chain_id,
                                 [[sses_mapping["F"], 1], [sses_mapping["F"], 0], [sses_mapping["G"], 1],
                                  [sses_mapping["G"], 0]])
    else:
        # if F is not present try to substitute  it by F'
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file,
                                 protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["E"], 1], [sses_mapping["G"], 1]],
                                 ["LINE", [sses_mapping["G"], 1], [sses_mapping["G"], 0]],
                                 ["LINE", [sses_mapping["G"], 0], [sses_mapping["F'"], 0]],
                                 ["LINE", [sses_mapping["F'"], 0], [sses_mapping["E"], 1]]])
    if ret == 0:
        return "3"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 3: {ret} {rev_m}")

    # 2d
    if test_is_sse_in_protein( protein_annotation_file, sses_mapping, "A") :
        # if A is present use it (def: egress between the N-terminus and helices a/A' and A)
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file,
                                protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["1_1"], 0], [sses_mapping["1_1"], 1]],
                                 ["LINE", [sses_mapping["1_1"], 1], [sses_mapping["A"], 0]],
                                 ["LOOP", [sses_mapping["A"], 0], [sses_mapping["1_1"], 0]]])
    elif test_is_sse_in_protein( protein_annotation_file, sses_mapping, "A'") :
        # if A is not present try to use A' instead
        ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file,
                                protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["1_1"], 0], [sses_mapping["1_1"], 1]],
                                 ["LINE", [sses_mapping["1_1"], 1], [sses_mapping["A'"], 0]],
                                 ["LINE", [sses_mapping["A'"], 0], [sses_mapping["A'"], 1]],
                                 ["LINE", [sses_mapping["A'"], 1], [sses_mapping["1_1"], 0]]])
    if ret == 0:
        return "2d"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 2d: {ret} {rev_m}")

    # opposite -------------------------------------------------------------------------------------
    # Water ./
    ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["B"], 0], [sses_mapping["B"], 1]],
                                 ["LINE", [sses_mapping["B"], 1], [sses_mapping["C"], 0]],
                                 ["LINE", [sses_mapping["C"], 0], [sses_mapping["C"], 1]],
                                 ["LINE", [sses_mapping["C"], 1], [sses_mapping["L"], 0]],
                                 ["LOOP", [sses_mapping["L"], 0], [sses_mapping["L"], 0, -13]],
                                 ["LINE", [sses_mapping["L"], 0, -13], [sses_mapping["B"], 0]]])
    if ret == 0:
        return "w"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders w: {ret} {rev_m}")

    # 1 ./
    ret = is_channels_in_cluster(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                 [[sses_mapping["C"], 0], [sses_mapping["L"], 0], [sses_mapping["L"], 1],
                                  [sses_mapping["H"], 0], [sses_mapping["H"], 1]])
    if ret == 0:
        return "1"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders 1: {ret} {rev_m}")

    # Water extended ./
    ret = is_tunnel_in_cluster2(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, chain_id,
                                [["LINE", [sses_mapping["B"], 0], [sses_mapping["B"], 1]],
                                 ["LINE", [sses_mapping["B"], 1], [sses_mapping["C"], 0]],
                                 ["LINE", [sses_mapping["C"], 0], [sses_mapping["C"], 1]],
                                 ["LINE", [sses_mapping["C"], 1], [sses_mapping["L"], 0]],
                                 ["LINE", [sses_mapping["L"], 0], [sses_mapping["L"], 1]],
                                 ["LINE", [sses_mapping["L"], 1], [sses_mapping["B"], 0]]])
    if ret == 0:
        return "w2"
    elif ret != 1:
        bugs.add(f"{tunnel_file} missing borders w: {ret} {rev_m}")

    if len(bugs) > 0:
        for bug in bugs:
            bugs_list.add(bug)
    return None


def add(ch1, ch2):
    r = []
    for c in ch1:
        r.append(c)
    for c in ch2:
        if c not in r:
            r.append(c)
    return r


# read arguments
parser = argparse.ArgumentParser(
    description='model is 1akdA00')
parser.add_argument('--tunnels_dir', required=True, action="store",
                    help="directory containing tunnel_*.pdb file (pdb/profile)")
parser.add_argument('--protein', required=True, action="store",
                    help="which channels should be detected")
parser.add_argument('--domain', required=True, action="store",
                    help="which channels should be detected")
parser.add_argument('--chain', required=True, action="store",
                    help="which channels should be detected")
parser.add_argument('--annotation_dir', required=True, action="store",
                    help="annotation file for wanted protein annotated with the same label as annotated model")
parser.add_argument('--cif_dir', required=True, action="store",
                    help="contains cif file for model and wanted protein all for wanted protein annotated with the same label as annotated model")
parser.add_argument('--output_fn', required=True, action="store",
                    help="output filename store the channel list")
parser.add_argument('--type', required=True, action="store",
                    help="database type af/pdb")
parser.add_argument('--part_only', action="store", default="all",
                    help="we want to compute tunnels only on 'solvent'/'water' tunnel part of HEM") 
#parser.add_argument('--pt', required=False, action="store",
#                    help="pores/tunnels/both") 
                        # used because ligands on the left part have to go directly to HEM center
                        # on water side can go to arbitrary part of HEM

args = parser.parse_args()

part_only = args.part_only
protein_cif_file = f"{args.cif_dir}/{args.protein}.cif"
cd_list = [[args.chain, args.domain]]

# filter too short and too long channels
if args.part_only != "other":
    short_tunnels, in_interval_tunnels, long_tunnels = short_ok_long_by_length(args.tunnels_dir, 10, 40)
else:
    short_tunnels, in_interval_tunnels, long_tunnels = short_ok_long_by_length(args.tunnels_dir, 10, 40, "all")

##############################################################")
# filter tunnels which do not start near <10 HEM
good_tunnels = {}
if args.part_only != "other":
  for cd in cd_list:
    if args.type == "af":
        good_tunnels[cd[0]] = in_interval_tunnels
    else:
        res = near_hem(args.tunnels_dir, in_interval_tunnels, protein_cif_file, cd[0])
        if res != []:
            good_tunnels[cd[0]] = res
else:
  for cd in cd_list:
    good_tunnels[cd[0]] = in_interval_tunnels

if len(good_tunnels) == 0:
    print(f"BUG = {protein_cif_file} {cd[0]} has 0 good tunnels")

# now we have to read domain data
# creating SSEs mapping - with model {MODEL}"
if args.part_only != "other":
    model_annotation_file = f"{args.annotation_dir}/{utils_1_10_630_10.MODEL[0]}{utils_1_10_630_10.MODEL[1]}{utils_1_10_630_10.MODEL[2]}-annotated.sses.json"
    sses_mapping = utils_1_10_630_10.sses_model_mapping_1og2(model_annotation_file)

protein_annotation_file = f"{args.annotation_dir}/{args.protein}{args.chain}{args.domain}-annotated.sses.json"
rotation_file = f"{args.annotation_dir}/{args.protein}{args.chain}{args.domain}-rotation.json"
bugs = set()

# rotate protein cif to be in tha same rotation as tunnels
out_dir = os.path.dirname(os.path.dirname(os.path.abspath(protein_annotation_file)))
os.makedirs(f"{out_dir}/cifs", exist_ok=True)
filename = f"{out_dir}/cifs/{args.protein}{args.chain}_rot.cif"
utils_read_cif.change_atoms_coords_with_rotation_file_label(protein_cif_file, rotation_file, filename, args.chain)

# for all tunnels find out whether it fits to any known type or not
for chain in cd_list:
    if chain[0] in good_tunnels:
        for tunnel in good_tunnels[chain[0]]:
            tunnel_file = f"{args.tunnels_dir}/{tunnel}"
            if args.part_only != "other":
                res = is_tunnel_of_known_type(tunnel_file, protein_annotation_file, rotation_file, protein_cif_file, cd[0], sses_mapping, bugs, part_only)
                if res == None:
                   # change the name of not sorted tunnels
                   res = "unclassified"
                if (part_only=="solvent") and (res in ["1","w"]):
                    # we are interested only in ["2a","2b", "2c", "2d", "2e", "2f", "3", "4", "5", "s"] tunnels 1 and water ones are sortes but unused for us
                    res="garbage"
                if (part_only=="water") and (res in ["2a","2b", "2c", "2d", "2e", "2f", "3", "4", "5", "s"]):
                    # we are interested only in ["1", "w"] other tunnels are sorted but unused for us
                    res="garbage"
                os.makedirs(f"{out_dir}/{res}", exist_ok=True)
                filename = f"{out_dir}/{res}/tunnel_{args.protein}{args.chain}{args.domain}_{tunnel[7:]}"
                utils_tunnel_pdb.change_tunnel_coords(tunnel_file, rotation_file, filename)
            else:
                res = "NO"
                print(tunnel, tunnel.split("_"))
                stunnel = tunnel.split("_")
                if len(stunnel) > 1:
                   os.makedirs(f"{out_dir}/{res}", exist_ok=True)
                   filename = f"{out_dir}/{res}/{stunnel[0]}_{args.protein}{args.chain}{args.domain}_{stunnel[1]}"
                   utils_tunnel_pdb.change_tunnel_coords(tunnel_file, rotation_file, filename)

# if all necessary protein SSEs (area borders) was not found
# and it can affects the counting (all tunnels was not detected), write a bug report about it
if len(bugs) > 0:
    for i in bugs:
        print(f"{args.protein} {i}")
