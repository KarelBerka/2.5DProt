#!/usr/bin/env python3
# create imput for molespecypic tool

import json
import os
import utils_tunnel_pdb
import utils_read_cif
import math
import numpy as np
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')


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
        for r in points:
            if (r["dist"] > i) and (r["dist"] < i+1):
                last = r["radius"]
                break
        print(last, end=" ")
        i += 1
    print()


def rotate_and_translate_position(r, t, p):
    return [float(r[0][0]*p[0]+r[1][0]*p[1]+r[2][0]*p[2]+t[0][0]),
            float(r[0][1]*p[0]+r[1][1]*p[1]+r[2][1]*p[2]+t[1][0]),
            float(r[0][2]*p[0]+r[1][2]*p[1]+r[2][2]*p[2]+t[2][0])]


def dist2d(a, b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


def resolve_cif_file(main_dir, pc_name):
    # resolve cif file for domain encoded as protein + chain + domain
    # PDB: 4-char protein id + 1-char chain + 2-char domain tag
    # AF: 36-char uuid protein + 1-char chain, empty domain tag
    candidates = [pc_name[:4], pc_name[:-3], pc_name[:-1], pc_name]
    for base in candidates:
        cif_file = f"{main_dir}/annotation/{base}.cif"
        if os.path.isfile(cif_file):
            return cif_file
    return None


def add_1d_data_to_tunnel_list(file_name, points, main_dir, pc_name):
    # create tunnel radius 1d graph data

    cif_file = resolve_cif_file(main_dir, pc_name)
    if cif_file is None:
        # cif not found - no HEM/ligand data to add
        return

    start = points[0]["coords"]
    chain = pc_name[4:-2]

    # read rotation and translation matrix
    rotation = f"{main_dir}/annotation/{pc_name}-rotation.json"
    with open(rotation, "r") as fr:
        r_data = json.load(fr)
        r = np.linalg.inv(r_data["rotation"])
        t = r_data["translation"]

    # reaf fe and transform it
    # read all HEMs from cif file
    fe_list = utils_read_cif.find_all_fe_hem_coords(cif_file, chain)
    if len(fe_list) == 0:
        # there is no HEM in this domain, just skip it
        return

    # transform all HEMs from cif file to the right coordinations
    transfer_fe_list = []
    for fe in fe_list:
        transf_fe = rotate_and_translate_position(r, t, fe)
        transfer_fe_list.append(transf_fe)

    # choose the nearest HEM to the point[0]
    transf_fe = None
    min_dist = 100
    for fe in transfer_fe_list:
        fe_point_dist = dist2d(points[0]['coords'], fe)
        if fe_point_dist < min_dist:
            min_dist = fe_point_dist
            transf_fe = fe

    if cif_file is not None:
        # read hetatom dict
        hetatm_dict = utils_read_cif.find_all_ligands_coords(cif_file, "ALL", ["HOH", "HEM"])
        transf_hetatm_dict = {}
        # rotate and translate it
        for h_id in hetatm_dict:
            transf_hetatm_dict[h_id] = {'name': hetatm_dict[h_id]['name'], 'coords': []}
            for c in hetatm_dict[h_id]['coords']:
                transf_c = rotate_and_translate_position(r, t, c)
                transf_hetatm_dict[h_id]['coords'].append(transf_c)

    for p in points:
        # add dist from HEM and ligand presence to tunnel point record
        p["dist_from_start"] = dist2d(p['coords'], start)
        try:
            p["dist_from_fe"] = dist2d(p['coords'], transf_fe)
        except:
            # cif file do not have HEM nearby to point[0] - transf_fe was empyt
            p["dist_from_fe"] = "-"

        p["ligand"] = []
        if cif_file is not None:
            # we have the ligands positions so we test weather they are in tunnel
            for ha in hetatm_dict:
                for coord in transf_hetatm_dict[ha]["coords"]:
                    if dist2d([p['coords'][0], p['coords'][1], p['coords'][2]], coord) < p['radius']:
                        # ligand is in tunnel, add its name to ligand list
                        p["ligand"].append(transf_hetatm_dict[ha]["name"])
                        break
    return


def count_min_radius(tunnel):
    t_min = 100
    for point in tunnel:
        if point['radius'] < t_min:
            t_min = point['radius']
    return t_min


def round_floats(o):
    if isinstance(o, float):
        return round(o, 2)
    if isinstance(o, dict):
        return {k: round_floats(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [round_floats(x) for x in o]
    return o


def found_is_similar(t1, t2):
    # test what is the maximal difference of tunnel t1 from t2
    tunnel_max_diff = 0
    for p1 in t1:
        c1 = p1["coords"]
        point_min_diff = 100
        # at first we have to find the nearest part of t2 from t1
        for p2 in t2:
            c2 = p2["coords"]
            d = math.sqrt((c1[0]-c2[0])*(c1[0]-c2[0])+(c1[1]-c2[1])*(c1[1]-c2[1])+(c1[2]-c2[2])*(c1[2]-c2[2]))
            if d < point_min_diff:
                point_min_diff = d
        # if is the nearest part further then max_dif, then we have to set bigger max_dif
        if tunnel_max_diff < point_min_diff:
            tunnel_max_diff = point_min_diff
    return tunnel_max_diff


def do_molespecypic_output(input_dir, name, uniprot_name):
    all_tunnels = {}               # auxiliary structure - used to remove duplicate None tunnels
                                   # suplicate to arbitrary already added one tunnel or other tunnel
                                   # even non added to data
    domains = []                   # list of all domains for each domain write here its rotation TODO
    data = {}
    data["uniprot_name"] = uniprot_name
    data["domain_rotation"] = {}
    annotation_dir = f"{input_dir}/annotation"
    for file_name in os.listdir(annotation_dir):
        if file_name.endswith("-rotation.json"):
            with open(f"{annotation_dir}/{file_name}") as fd:
                rot_data = json.load(fd)
            pcd = file_name[:-len("-rotation.json")]
            domains.append(pcd)
            data["domain_rotation"][pcd] = rot_data
    data["domains"] = domains

    # for each domain write here all uniq tunnels
    data["tunnels"] = {}
    all_tunnels = []
    for tunnel_type in ["2a", "2b", "2c", "2e", "2f", "s", "1", "w", "2d", "3", "4", "5", "unclassified", "w2"]:
        # we have to  finde tunnels for all tunnel sets 
        data["tunnels"][tunnel_type] = []
        tunnel_type_dir = f"{input_dir}/{tunnel_type}"

        if os.path.isdir(tunnel_type_dir):
            for domain in domains:
                # for all domains - find only the best tunnel
                # go through all and search the best one
                best_tunnel_record = None
                best_min_radius = 0
                for file_name in os.listdir(tunnel_type_dir):
                  if file_name.startswith(f"tunnel_{domain}"):
                    split_file_name = file_name.split("_")
                    tunnel_info = utils_tunnel_pdb.get_tunnel_points_with_radius(f"{tunnel_type_dir}/{file_name}")
                    min_radius = count_min_radius(tunnel_info)
                    add_1d_data_to_tunnel_list(f"{tunnel_type_dir}/{file_name}", tunnel_info, f"{input_dir}", split_file_name[1])
                    tunnel_record = {'domain_name': split_file_name[1],
                                     'min_radius': min_radius,
                                     'points': tunnel_info,
                                     'tt': tunnel_type}
                    if tunnel_type == "unclassified":
                        # all None tunnels have to be added
                        # we have to sort duplicate tunnels (thus which are <2A from another one)
                        is_similar = False
                        if split_file_name[2][0] == "6":
                            # tunnels from round 6 are probably broken (are not in w and w2 and 1, thus should be on wrong membran side
                            # just do not add them to unclassified
                            continue
                        for tu in all_tunnels:
                            is_similar = False
                            if tu["domain_name"] == split_file_name[1]:
                                # already added
                                # check whether they are not similar as the current one:
                                    max_diff = found_is_similar(tunnel_info, tu['points'])
                                    if max_diff < 2.5:
                                        is_similar = True
                                        break
                        if is_similar == False and tunnel_record["min_radius"] >0.95:
                            data["tunnels"][tunnel_type].append(tunnel_record)
                    else:
                        # only the widest tunnel from 1, 2a, 2b, ...  is added
                        if min_radius > best_min_radius:
                            best_tunnel_record = tunnel_record
                            best_min_radius = min_radius
                    # but for comparison we have to save all tunnels
                    all_tunnels.append(tunnel_record)
                if best_tunnel_record != None:
                    data["tunnels"][tunnel_type].append(best_tunnel_record)

    print(f" created: {name}")
    with open(f"{name}", "w") as fd:
         print(json.dump(round_floats(data), fd))


