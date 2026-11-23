#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import math
import json
import numpy as np

def dist3d(a, b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


def read_tunnel_length(pdb_file):
    # return tunnel length
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                aux = line
            line = fpdb.readline()
    return float(aux[54:60])


def get_sten_distance(pdb_file):
    # return tunnel length
    start = None
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                if start == None:
                    start = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                aux = line
            line = fpdb.readline()
    end = [float(aux[30:38]), float(aux[38:46]), float(aux[46:54])]
    return dist3d(start, end)


def get_tunnel_min_radius(pdb_file, ommit_start_part=0, ommit_end_part=0):
    # return tunnel length
    min_r = 100
    tunnel_length = read_tunnel_length(pdb_file)
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                # X,Y,Z coords can be too short
                try:
                    distance = float(line[54:60])
                    radius = float(line[60:66])
                except:
                    distance = float(line[54:58])
                    radius = float(line[58:63])
                if (distance > ommit_start_part) and (distance < tunnel_length - ommit_end_part):
                    # we ommit the frst {ommit_start_part} angstrems of tunnels
                    if radius < min_r:
                        min_r = radius
            line = fpdb.readline()
    return min_r


def get_tunnel_points(pdb_file):
    # return tunnel length
    points = []
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                new_point = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                if (points == []) or (dist3d(points[-1],new_point) > 1):
                    points.append(new_point)
            line = fpdb.readline()
    return points


def get_tunnel_points_with_radius(pdb_file):
    # return tunnel length
    points = []
    difference = 0
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                new_point = {}
                new_point['coords'] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                new_point['radius'] = float(line[62:66])

                # get only one point per in each 1A dist
                if (points == []) or (dist3d(points[0]['coords'], new_point['coords']) > difference):
                    # if there is no point, put here the start one
                    points.append(new_point)
                    difference += 1

                    # if the distance is bigger then 1A put here the point several times
                    while (dist3d(points[0]['coords'], new_point['coords']) > difference):
                        points.append(new_point)
                        difference += 1
            line = fpdb.readline()
    return points


def get_all_tunnel_points(pdb_file):
    # return tunnel length
    points = []
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                new_point = {'x':float(line[30:38]),
                             'y':float(line[38:46]),
                             'z':float(line[46:54]),
                             'dist':float(line[54:60]),
                             'radius':float(line[60:])}
                points.append(new_point)
            line = fpdb.readline()
    return points


def read_tunnel_start(pdb_file):
    # return coordinates of tunnel start
    with open(pdb_file, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                return [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            line = fpdb.readline()
    return None


def change_tunnel_coords(pdb_file, rotation_file, output_file):
    # read rotation files - to rotate tunnels
    with open(f"{rotation_file}", "r") as fr:
        r_data = json.load(fr)
        r = np.linalg.inv(r_data["rotation"])
        t = r_data["translation"]

    points = []
    with open(pdb_file, "r") as fpdb, open(output_file, "w") as fout:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                sv = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                transf_sv = [r[0][0]*sv[0]+r[1][0]*sv[1]+r[2][0]*sv[2]+t[0][0], \
                     r[0][1]*sv[0]+r[1][1]*sv[1]+r[2][1]*sv[2]+t[1][0], \
                     r[0][2]*sv[0]+r[1][2]*sv[1]+r[2][2]*sv[2]+t[2][0]]

                cxs = line.find(line[30:38])
                cxl = len(line[30:38])
                cys = line[cxs+cxl:].find(line[38:46])+cxs+cxl
                cyl = len(line[38:46])
                czs = line[cys+cyl:].find(line[46:54])+cys+cyl
                czl = len(line[46:54])
                new_line = f"{line[:cxs]}{str(transf_sv[0])[:cxl]}{line[cxs+cxl:cys]}{str(transf_sv[1])[:cyl]}{line[cys+cyl:czs]}{str(transf_sv[2])[:czl]}{line[czs+czl:]}"
                fout.write(new_line)

            else:
                fout.write(line)
            line = fpdb.readline()


def test_are_tunnels_almost_same(pdb_file1, pdb_file2, toleration=2):
    # returns true if tunnel in pdb_file1 and pdb_file2 are the same =
    # one can have points at least "toleratoin" distant from the other one
    # read the points for the first tunnel
    points1 = []
    with open(pdb_file1, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                new_point = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                points1.append(new_point)
            line = fpdb.readline()

    # read the points for the second tunnel
    points2 = []
    with open(pdb_file2, "r") as fpdb:
        line = fpdb.readline()
        while line:
            # structure:
            # ATOM    NAM  RES TUNID          X       Y       Z      Distnm RadiusA 
            # HETATM  105  X   TUN H   8      47.487  65.494  20.357 12.87  2.22              
            # HETATM  106  X   TUN H   8      47.479  65.749  20.483 13.16  2.34              
            if line.startswith("HETATM"):
                new_point = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                points2.append(new_point)
            line = fpdb.readline()

    # find wheather the first tunnel pass the condition - if no just return False
    for i in points1:
        is_near = False
        for j in points2:
            if dist3d(i,j) < toleration:
                is_near = True
        if is_near == False:
            return False

    # find wheather the second tunnel pass the condition - if no just return False

    for i in points2:
        is_near = False
        for j in points1:
            if dist3d(i,j) < toleration:
                is_near = True
        if is_near == False:
            return False

    # both passes the condition
    return True


