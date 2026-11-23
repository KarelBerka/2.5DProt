#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import json
import utils_read_cif
import numpy as np
import math

#MODEL = ["1ap9", "A", "00"]
MODEL = ["2vt4", "A", "00"]

def sses_model_mapping_1at9(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1ap9
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 8 and sse["end"] == 30:
                mapping["A"] = sse["label"]
            if sse["start"] == 37 and sse["end"] == 62:
                mapping["B"] = sse["label"]
            if sse["start"] == 80 and sse["end"] == 101:
                mapping["C"] = sse["label"]
            if sse["start"] == 105 and sse["end"] == 127:
                mapping["D"] = sse["label"]
            if sse["start"] == 134 and sse["end"] == 159:
                mapping["E"] = sse["label"]
            if sse["start"] == 165 and sse["end"] == 191:
                mapping["F"] = sse["label"]
            if sse["start"] == 201 and sse["end"] == 224:
                mapping["G"] = sse["label"]
    return mapping


def sses_model_mapping_2vt4(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-2vt4
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 21 and sse["end"] == 39:
                mapping["A"] = sse["label"]
            if sse["start"] == 45 and sse["end"] == 74:
                mapping["B"] = sse["label"]
            if sse["start"] == 81 and sse["end"] == 115:
                mapping["C"] = sse["label"]
            if sse["start"] == 125 and sse["end"] == 149:
                mapping["D"] = sse["label"]
            if sse["start"] == 175 and sse["end"] == 207:
                mapping["E"] = sse["label"]
            if sse["start"] == 225 and sse["end"] == 256:
                mapping["F"] = sse["label"]
            if sse["start"] == 262 and sse["end"] == 286:
                mapping["G"] = sse["label"]
    return mapping


def predict_1_20_1070_10_tunnel_start_end(annotation_dir, pcd):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_2vt4(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    s1, s2 = None, None
    e1, e2 = None, None
    l = list(data.keys())[0]
    for sse in data[l]["secondary_structure_elements"]:
#        if sse['label'] == mapping["C"]:
#            s1 = sse["start"]
#            e1 = sse["end"]
#        if sse['label'] == mapping["G"]:
#            s2 = sse["start"]
#            e2 = sse["end"]
        if sse['label'] == mapping["C"]:
            s1 = sse["start"]
        if sse['label'] == mapping["B"]:
            e1 = sse["start"]
        if sse['label'] == mapping["G"]:
            s2 = sse["start"]
        if sse['label'] == mapping["F"]:
            e2 = sse["start"]

    if s1 != None and s2 != None and e1 != None and e2 != None:
        # transform start and end residue to coordinate in the original (not rotated) cif file
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
        start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
        start = [(start1[0]+start2[0])/2, (start1[1]+start2[1])/2, (start1[2]+start2[2])/2]

        end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
        end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
        end = [(end1[0]+end2[0])/2, (end1[1]+end2[1])/2, (end1[2]+end2[2])/2]
        return start, end
    else:
        return None, None


def two_d_dist(a,b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


def center_coord(s):
    counter = 0
    coord = [0,0,0]
    for p in s:
        coord[0] += p[0]
        coord[1] += p[1]
        coord[2] += p[2]
        counter += 1
    return [coord[0]/counter, coord[1]/counter, coord[2]/counter ]


def predict_1_20_1070_10_ec_binding_point_from_ligand(annotation_dir, pcd):
    # ligand C90 is in binding place - use it
    cif_file = f"{annotation_dir}/{pcd[0]}.cif"
    model_ligands = utils_read_cif.find_all_ligands_coords(cif_file, "ALL", "")

    for i in model_ligands:
        # find coordination of ligand "P32" or "C90", ... all of them denotes binding places
        if model_ligands[i]["name"] in ["P32", "G90", "CVD", "Y00", "WHJ", "5FW", "68H", "XF5", "H98", "3WC", "RET", "ERM", "CAU", "7LD", "H8G", "H8D", "H8J", "H8M", "A1EB7", "A8T", "SG8", "VIS", "V0S", "2GM", "89F", "EP5", "836", "C01",  "0HK", "3C0", "9EC", "1KS", "SNT", "CY8", "I32", "XTK", "FVK", "ANR"]:
            # find all res from protein pcd[0], chain[1]
            pcd_res = utils_read_cif.read_all_atoms_coords(cif_file, None, pcd[1])
            # test if the ligand model_ligand[i] dist from res is <5 -- the ligand is in this domain
            mmin  = 200
            for atom in pcd_res:
              for point in model_ligands[i]["coords"]:
                dist = two_d_dist(pcd_res[atom], point)
                if dist < mmin:
                    mmin = dist
                if dist < 6:
                     return center_coord(model_ligands[i]["coords"])
    return None


def predict_1_20_1070_10_ec_binding_point_from_ligand_p07700(annotation_dir, pcd):
    # find ligand P32 from P07700 protein pvt4
    # this protein is in binding place
    # use its center as a start of tunnel 

    model_cif = f"{annotation_dir}/{MODEL[0]}.cif"
    model_ligands = utils_read_cif.find_all_ligands_coords(model_cif, MODEL[1], "")
    lig_coords = [0, 0, 0]
    counter = 0
    for i in model_ligands:
        # find coordination of ligand P32 in P07700
        if model_ligands[i]["name"] == "P32":
            for j in model_ligands[i]["coords"]:
                 lig_coords = [lig_coords[0]+j[0], lig_coords[1]+j[1], lig_coords[2]+j[2]]
                 counter += 1
    lig_coords = [lig_coords[0]/counter, lig_coords[1]/counter, lig_coords[2]/counter]

    # rotate it to unified domain rotation from overprot
    p07700_rotation=f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-rotation.json"
    with open(p07700_rotation, "r") as fr:
        r_data = json.load(fr)
        m_r = np.array(r_data["rotation"])
        m_t = np.array(r_data["translation"])
    transf_p07700_ligand = np.dot(m_r, np.array(lig_coords)) + m_t.T

    # rotate it to pcd protein rotation
    rotation_json=f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-rotation.json"
    with open(rotation_json, "r") as fr:
        r_data = json.load(fr)
        r = np.array(r_data["rotation"])
        t = np.array(r_data["translation"])

    ret = np.dot(np.linalg.inv(r), transf_p07700_ligand.T-t)
    rfl = list(ret.flatten())

    # return this value
    return [float(rfl[0]), float(rfl[1]), float(rfl[2])]


def compute3d_point_from_two_points(p1, p2, ratio):
    x = p1[0]+(p2[0]-p1[0])*ratio
    y = p1[1]+(p2[1]-p1[1])*ratio
    z = p1[2]+(p2[2]-p1[2])*ratio
    return [x, y, z]


def is_near_atom(p, cif_file, limit = 1.2):
    lst = utils_read_cif.read_all_atoms_coords(cif_file)
#    print(lst)
    for i in lst:
#       print("l", l)
       l = lst[i]
       dif = (float(l[0])-p[0])*(float(l[0])-p[0])+(float(l[1])-p[1])*(float(l[1])-p[1])+(float(l[2])-p[2])*(float(l[2])-p[2])
#       print(int(dif*10), p, l)
       if dif < limit*limit:
           return True
    return False

def predict_1_20_1070_10_ligand_position(annotation_dir, pcd):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_2vt4(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    point1, point2, point3 = None, None, None
    cif_file = f"{annotation_dir}/{pcd[0]}.cif"

    l = list(data.keys())[0]
    for sse in data[l]["secondary_structure_elements"]:
#        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        if sse['label'] == mapping["F"]:
            p1_start = sse["end"]
            point1_start = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p1_start)
            p1_end = sse["start"]
            point1_end = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p1_end)
            ratio = 10/27
            point1 = compute3d_point_from_two_points(point1_start, point1_end, ratio)

        if sse['label'] == mapping["G"]:
            p2_start = sse["start"]
            point2_start = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p2_start)
            p2_end = sse["end"]
            point2_end = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p2_end)
            ratio = 14/33
            point2 = compute3d_point_from_two_points(point2_start, point2_end, ratio)

        if sse['label'] == mapping["C"]:
            p3_start = sse["start"]
            point3_start = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p3_start)
            p3_end = sse["end"]
            point3_end = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p3_end)
            ratio = 7/22
            point3 = compute3d_point_from_two_points(point3_start, point3_end, ratio)
    if point1==None or point2==None or (point3 == None):
        return None
    else:
        # transform start and end residue to coordinate in the original (not rotated) cif file
#        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
#        point1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p1)
#        point2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p2)
#        point3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], p3)
        point = [(point1[0]+point2[0]+point3[0])/3, (point1[1]+point2[1]+point3[1])/3, (point1[2]+point2[2]+point3[2])/3]
        while is_near_atom(point, f"{annotation_dir}/{pcd[0]}.cif"):
            pointfg = [(point1[0]+point2[0])/2, (point1[1]+point2[1])/2, (point1[2]+point2[2])/2]
            vect = [point3[0]-pointfg[0], point3[1]-pointfg[1], point3[2]-pointfg[2]]
            point = [point[0]+vect[0]*0.02, point[1]+vect[1]*0.02, point[2]+vect[2]*0.02]
        return point


def test_domain_completeness(domain_annotation_file, mapping, percentage=80):
    # test whether input annotation file contains at least "percentage" percent of 1.20.1070.10 main SSEs
    count_found = 0
    with open(domain_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = list(data.keys())[0]
        wanted_list = ["A", "B", "C", "D", "E", "F", "G"]
        for wanted in wanted_list:
            for sse in data[pcd]["secondary_structure_elements"]:
                if sse["label"] == mapping[wanted]:
                    count_found += 1
    if count_found > len(wanted_list)* percentage/100:
        return True
    else:
        return False
