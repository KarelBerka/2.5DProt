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

    t = -(D + np.dot(triangle_norm, s1)) / (np.dot(triangle_norm, (s2-s1)))   # intersection of triangle and line

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
#                print(f"pass - return {ret}: sses: {sses} channel part: {counter} {l}  {[start, end]}")
        counter += 1
    # there were no crossing
    return fin_ret


def test_line_cross_triangle():
    # line cross triangle test function
    line = [[0, 0, 0], [2, 0, 0]]
    triangle = [[1, 1, 1], [1, 1, -1], [1, -1, 0]]
    ret = line_cross_triangle(triangle, line)
    print(f"line cross triangle line:{line} triangle:{triangle}, ret:{ret} - should be True")

    line = [[-2, 0, 0], [2, 0, 0]]
    triangle = [[1, 1, 1], [1, 1, -1], [-1, -1, 0]]
    ret = line_cross_triangle(triangle, line)
    print(f"line cross triangle line:{line} triangle:{triangle}, ret:{ret} - should be True")

    line = [[0, 0, 0], [-2, 0, 0]]
    triangle = [[1, 1, 1], [1, 1, -1], [1, -1, 0]]
    ret = line_cross_triangle(triangle, line)
    print(f"line cross triangle line:{line} triangle:{triangle}, ret:{ret} - should be False")


def count_file_lines(fn):
    # count file lines
    l_count = 0
    with open(fn, "r") as fr:
        line = fr.readline()
        while line:
            l_count += 1
            line = fr.readline()
    return l_count


def short_ok_long_by_lines(gdir, min_ratio, max_ratio):
    # filter too short < MIN_RATIO * avg
    # and long > MAX_RATIO * avg channels
    # ength defined by number of lines in channel files
    lines_sum = 0
    files_count = 0
    # at first set the average number of lines of channel files
    for file_name in os.listdir(gdir):
        if file_name.startswith("gnuplot-"):
            lines_sum += count_file_lines(f"{gdir}/{file_name}")
            files_count += 1
    # avg channel length
    avg = lines_sum/files_count
    short_channels = []
    long_channels = []
    ok_channels = []
    for file_name in os.listdir(gdir):
        if file_name.startswith("gnuplot-"):
            counter = count_file_lines(f"{gdir}/{file_name}")
            if counter < avg * min_ratio:
                short_channels.append(file_name)
            elif counter < avg * max_ratio:
                ok_channels.append(file_name)
            else:
                long_channels.append(file_name)

    return short_channels, ok_channels, long_channels


def dist2d(a, b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


def near_hem(gdir, adir, channel_list, dist=10):
    # filter channels which starts far from hem
    ok = []
    bad = []
    for ch in channel_list:  # ch is e.g. 'gnuplot-4wnv00B.txt_7'
        chsplit = ch.split(".")   # remove suffix
        chss = chsplit[0].split("-")
        ligand_list = f"{adir}/{chss[1]}-ligands.json"
        hem_pos = None
        with open(ligand_list, "r") as fl:
            chain = chss[1][6:]
            data = json.load(fl)
            for lig in data:
                # go through protein ligands and find the hem
                if (data[lig]["auth_chain"] == chain) and (data[lig]['entity_comp'] == "HEM"):
                    hem_pos = data[lig]['center']

        if hem_pos is None:
            # hem not found
            bad.append(ch)
            continue
        with open(f"{gdir}/{ch}", "r") as fg:
            line = fg.readline()
            sline = line.split(" ")
            start = [float(sline[0]), float(sline[1]), float(sline[2])]

            d = dist2d(start, hem_pos)
            print("dist from hem", d, ch)
            if d < 10:
                ok.append(ch)
            else:
                bad.append(ch)
    return ok, bad


def get_sses_set(annot_fn, pcd):
    sses_list = {}
    # TODO - filtr only main sses
    with open(f"{annot_fn}", "r") as fr:
        data = json.load(fr)
        for sse in data[pcd]["secondary_structure_elements"]:
            sses_list[sse["label"]] = {"start": sse["start_vector"], "end": sse["end_vector"]}
    return sses_list


def find_channels_from_cluster(output_fn, channels_set, channel_dir, annot_dir, ommit_list, limit_max, sses_points, name, output=True):
    # near and ommitd by the vprevious version
    cluster = []
    complement = []
    missing_SSEs = []

    print("find_channels_from_the_same_cluster: start cluster generating with limitations", sses_points, name)
    counter = 0

    for c in channels_set:
        counter += 1
        aux = c.split(".")[0][8:]
        pcd = f"{aux[0:4]}{aux[6:]}{aux[4:6]}"
        annot_f = f"{annot_dir}/{pcd}-annotated.sses.json"
        channel_f = f"{channel_dir}/{c}"
        sses_set = get_sses_set(annot_f, pcd)
        sses_missing = False
        s = []

        for res in sses_points:
            if res[0] not in sses_set:
                sses_missing = True
                break
            if res[1] == 0:
                r = sses_set[res[0]]["start"]
            elif res[1] == 1:
                r = sses_set[res[0]]["end"]
            s.append(r)

        if sses_missing:
            missing_SSEs.append(c)
            if output:
                print(f"{counter} - {c} {sses_points} {name} - missing")
                with open(f"{output_fn}_missing", "a") as f:
                    f.write(f"'{channel_f}' lt rgb '{name}_sses_missing'\n")
            continue

        c_points = []
        with open(channel_f, "r") as fr:
            line = fr.readline()
            while line:
                sline = line.split(" ")
                c_points.append([float(sline[0]), float(sline[1]), float(sline[2])])
                line = fr.readline()

        res = channel_cross_polygon(s, c_points)
        if res == 1:
            if output:
                print(f"{counter} - {c} {sses_points} {name} - proslo")
                with open(f"{output_fn}", "a") as f:
                    f.write(f"'{channel_f}' lt rgb '{name}_ok{res}'\n")
            cluster.append(c)
        else:
            print(f"{counter} - {c} {sses_points} {name} - neproslo")
            complement.append(c)
    print(f"cluster {len(cluster)}, complement {len(complement)}, sses_miss {len(missing_SSEs)}")
    return cluster, complement, missing_SSEs


def write_all(channels_set):
    # near and ommitd by the vprevious version
    cluster = []
    complement = []
    missing_SSEs = []

    for c in channels_set:
        cluster.append(c)
    print(f"cluster {len(cluster)}, complement {len(complement)}, sses_miss {len(missing_SSEs)}")
    return cluster, complement, missing_SSEs





def output(ch, channel_dir, name, suffix):
    for c in ch:
        channel_f = f"{channel_dir}/{c}"
        print(f"'{channel_f}' lt rgb '{name}_{suffix}'")


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
    description='TODO')
parser.add_argument('--dir', required=True, action="store",
                    help="channel files directory (gnuplot-{family_nme}-{timestamp})")
parser.add_argument('--annotation_dir', required=True, action="store",
                    help="annotation files directory generated-{family_nme}-{timestamp}/annotation")
parser.add_argument('--ommit', action="store",
                    help="list of ommited domains format: TODO something like 'bla/bla/name' ommit the name ")
parser.add_argument('--limit_max', action="store",
                    help="maximal channel length")
parser.add_argument('--output_fn', required=True, action="store",
                    help="output filename store the channel list")
parser.add_argument('--chain_name', required=True, action="store",
                    help="like 1,2a,water")
args = parser.parse_args()

ch_list = {}
channels_clusters = {}

# filter too short and too long channels
short_ch, ok_ch, long_ch = short_ok_long_by_lines(args.dir, 0.2, 3.5)

# filter chanels which do not start near <10 HEM
ok2_ch, far2_ch = near_hem(args.dir, args.annotation_dir, ok_ch, dist=10)

print("################################################################")
print("all channels", len(short_ch) + len(ok_ch) + len(long_ch))
print(f"good length {len(ok_ch)} (short - {len(short_ch)}, long - {len(long_ch)}) ")
print(f"starts near hem {len(ok2_ch)} (far - {len(far2_ch)})")
print("################################################################")
ok_ch = ok2_ch


# omit the other needed ones
ommit_list = []
if args.ommit:
    with open(args.ommit, "r") as fr:
        line_counter = 1
        point_counter = 0
        line = fr.readline()
        while line:
            sline = line.split("'")
            ssline = sline[1].split("/")
            ommit_list.append(ssline[2])
            line = fr.readline()

if (args.limit_max):
    limit_max = args.limit_max
else:
    limit_max = 100000

# SSES list in overprot:
# H13 - ?? A
# H22 - B
# H28 - B'
# H36 - C
# H41 - D
# H46 - E
# H49 - F
# H59 - G
# H63 - H
# H68 - I
# H70 - J
# H73 - K
# H94 - L ??

gdir = args.dir
annotation_dir = args.annotation_dir

cluster, complement, missing = write_all(ok_ch)
output(cluster, args.dir, "channel2c", "-")

# Channel 2c

# H28(1) - B', H68(0) - I, (H59 (1)) -G, H59(0) - G,
if args.channel_name == "2c":
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H28", 1], ["H68", 0], ["H59", 1], ["H59", 0]], "channel2c", output=True)
    output(cluster, args.dir, "channel2c", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H68", 0], ["H59", 1], ["H59", 0], ["H28", 1]], "channel2c", output=True)
    output(cluster2, args.dir, "channel2c", "-")

    print("channel2c")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - all, missing", len(missing))

#sys.exit()
# Channel 2e

# H28(0) - B', H28(1) - B', (H22 (1)) - B, H36(0) - C,
if args.channel_name == "2e":
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H22", 1], ["H28", 0], ["H28", 1], ["H36", 0]], "channel2e", output=True)
    output(cluster, args.dir, "channel2e", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H28", 0], ["H28", 1], ["H36", 0], ["H22", 1]], "channel2e", output=True)
    output(cluster2, args.dir, "channel2e", "-")

    print("channel2e")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - all, missing", len(missing))

# H30(0) - B', H30(1) - B', (H22 (1)) - B, H36(0) - C,
if False:
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H22", 1], ["H30", 0], ["H30", 1], ["H36", 0]], "channel2e", output=True)
    output(cluster, args.dir, "channel2e", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H30", 0], ["H30", 1], ["H36", 0], ["H22", 1]], "channel2e", output=True)
    output(cluster2, args.dir, "channel2e", "-")

    print("channel2e")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - all, missing", len(missing))

# channel 3 - mezi F (H49) a G(H59) -> G0, G1, H0, H1
if args.channel_name == "3":
    feloop = True     # add the space between E->F loop as well
    if feloop:
        cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H59", 0], ["H59", 1], ["H46", 1], ["H49", 0], ["H49", 1]], "channel3", output=True)
        output(cluster, args.dir, "channel3", "H59,0-H59,1-H46,1-H49,0-H49,1")

        cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H59", 1], ["H46", 1], ["H49", 0], ["H49", 1], ["H59", 0]], "channel3", output=True)
        output(cluster2, args.dir, "channel3", "H59,1-H46,1-H49,0-H49,1-H59,0")

        cluster3, complement3, missing3 = find_channels_from_cluster(args.output_fn, complement2, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H46", 1], ["H49", 0], ["H49", 1], ["H59", 0], ["H59", 1]], "channel3", output=True)
        output(cluster3, args.dir, "channel3", "H59,1-H46,1-H49,0-H49,1-H59,0")
        print("channel3")
        print("lens - ok ", len(cluster), len(complement), len(missing))
        print("lens - ok ", len(cluster2), len(complement2), len(missing2))
        print("lens - ok ", len(cluster3), len(complement3), len(missing3))
        print("lens - all, missing", len(missing))

    else:
        cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H59", 0], ["H59", 1], ["H49", 0], ["H49", 1]], "channel3", output=True)
        output(cluster, args.dir, "channel3", "1")

        cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H59", 1], ["H49", 0], ["H49", 1], ["H59", 0]], "channel3", output=True)
        output(cluster2, args.dir, "channel3", "2")

        print("channel3")
        print("lens - ok ", len(cluster), len(complement), len(missing))
        print("lens - ok ", len(cluster2), len(complement2), len(missing2))
        print("lens - all, missing", len(missing))

# channel 1 - mezi H (H63) a C(H36) a L(H94) -> H0, H1, C0, L0, L1
if args.channel_name == "1":
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H36", 0], ["H94", 0], ["H94", 1], ["H63", 0], ["H63", 1]], "channel1", output=True)
    output(cluster, args.dir, "channel1", "H36,0-H94,0-H94,1-H63,0-H63,1")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H63", 1], ["H36", 0], ["H94", 0], ["H94", 1], ["H63", 0]], "channel1", output=True)
    output(cluster2, args.dir, "channel1", "H63,1-H36,0-H94,0-H94,1-H63,0")

    cluster3, complement3, missing3 = find_channels_from_cluster(args.output_fn, complement2, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H63", 0], ["H63", 1], ["H36", 0], ["H94", 0], ["H94", 1]], "channel1", output=True)
    output(cluster3, args.dir, "channel1", "H63,0-H63,1-H36,0-H94,0-H94,1")

    print("channel1")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - ok ", len(cluster3), len(complement3), len(missing3))
    print("lens - all, missing", len(missing))
#sys.exit()

# channel 4 - mezi F zadek (H49) a G predek(H59) a H54(??)  -> H49,1; H54,0; H54,1; H59,0
if args.channel_name == "4":
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H49", 1], ["H54", 0], ["H54", 1], ["H59", 0]], "channel4", output=True)
    output(cluster, args.dir, "channel4", "H49,1-H54,0-H54,1-H59,0")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H54", 0], ["H54", 1], ["H59", 0], ["H49", 1]], "channel4", output=True)
    output(cluster2, args.dir, "channel4", "H54,0-H54,1-H59,0-H49,1")

    print("channel4")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - all, missing", len(missing))

    ok_ch = add(complement2, missing)

# channel 4 - mezi F zadek (H49) a G predek(H59) a H28(??) - druha verze  -> H49,1; H28,1; H28,0; H59,0
#if False:
#    # todo upravit
#    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
#                                                      [["H49", 1], ["H28", 1], ["H28", 0], ["H59", 0]], "channel4", output=True)
#    output(cluster, args.dir, "channel4", "H49,1-H28,1-H28,0-H59,0")
#
#    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
#                                                          [["H59", 0], ["H49", 1], ["H28", 1], ["H28", 0]], "channel4", output=True)
#    output(cluster2, args.dir, "channel4", "H59,0-H49,1-H28,1-H28,0")
#
#    print("channel4b")
#    print("lens - ok ", len(cluster), len(complement), len(missing))
#    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
#    print("lens - all, missing", len(missing))
#

# Water  channel  - zacatek 1 (H36) konec H36, konec 73  konec 114 -> konec 355 -  H36,0; H36,1; E78,1; H22,1
if False:
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H36", 1], ["H36", 0], ["H22", 1], ["H22", 0], ["H73", 1]], "channelw", output=True)
    output(cluster, args.dir, "channelw", "H36,1-H36,0-H22,1-H22,0-H73,1")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H73", 1], ["H36", 1], ["H36", 0], ["H22", 1], ["H22", 0]], "channelw", output=True)
    output(cluster2, args.dir, "channelw", "H73,1-H36,1-H36,0-H22,1-H22,0")

    cluster3, complement3, missing3 = find_channels_from_cluster(args.output_fn, complement2, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H22", 0], ["H73", 1], ["H36", 1], ["H36", 0], ["H22", 1]], "channelw", output=True)
    output(cluster3, args.dir, "channelw", "H22,0-H73,1-H36,1-H36,0-H22,1")
    print("channelw")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - ok ", len(cluster3), len(complement3), len(missing3))
    print("lens - all, missing", len(missing))


# Solvent channel - nutno udelat intersection from both
# Ek, Ez, FZ, Fk E105(cokoli)
if False:
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H46", 1], ["H46", 0], ["H49", 0], ["H49", 1], ["E105", 0]], "channels", output=True)
    output(cluster, args.dir, "channels", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H46", 0], ["H49", 0], ["H49", 1], ["E105", 0], ["H46", 1]], "channels", output=True)
    output(cluster2, args.dir, "channels", "-")

    cluster3, complement3, missing3 = find_channels_from_cluster(args.output_fn, complement2, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H49", 0], ["H49", 1], ["E105", 0], ["H46", 1], ["H46", 0]], "channels", output=True)
    output(cluster3, args.dir, "channels", "-")
    print("channels")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - ok ", len(cluster3), len(complement3), len(missing3))
    print("lens - all, missing", len(missing))

# Ik, Iz, FZ, Fk E105(cokoli)
if False:
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H68", 1], ["H68", 0], ["H49", 0], ["H49", 1], ["E105", 0]], "channels", output=True)
    output(cluster, args.dir, "channels", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H68", 0], ["H49", 0], ["H49", 1], ["E105", 0], ["H68", 1]], "channelw", output=True)
    output(cluster2, args.dir, "channels", "-")

    cluster3, complement3, missing3 = find_channels_from_cluster(args.output_fn, complement2, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H49", 0], ["H49", 1], ["E105", 0], ["H68", 1], ["H68", 0]], "channelw", output=True)
    output(cluster3, args.dir, "channels", "-")
    print("channels")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - ok ", len(cluster3), len(complement3), len(missing3))
    print("lens - all, missing", len(missing))

# channel 5
# K(H73) a K'(H83) - it is too small add A(H13),0
if True:
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H73", 0], ["H73", 1], ["H13", 0], ["H83", 0], ["H83", 1]], "channel5", output=True)
    output(cluster, args.dir, "channel4", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H83", 1], ["H73", 0], ["H73", 1], ["H13", 0], ["H83", 0]], "channel5", output=True)
    output(cluster2, args.dir, "channel4", "-")

    cluster3, complement3, missing3 = find_channels_from_cluster(args.output_fn, complement2, gdir, annotation_dir, ommit_list, limit_max,
                                                          [["H83", 0], ["H83", 1], ["H73", 0], ["H73", 1], ["H13", 0]], "channel5", output=True)
    output(cluster3, args.dir, "channel4", "-")

    print("channel5")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - ok ", len(cluster3), len(complement3), len(missing3))
    print("lens - all, missing", len(missing))


# Channel 2b
# H28(0) - B' (H22 - B - podle me blbe), E20 (0), E78(1)
# H56(0) (nebo H54,1)
if False:
    # todo upravit
    cluster, complement, missing = find_channels_from_cluster(args.output_fn, ok_ch, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["H28", 0], ["H59", 0], ["E20", 0], ["E78", 1]], "channel2b", output=True)
    output(cluster, args.dir, "channel2b", "-")

    cluster2, complement2, missing2 = find_channels_from_cluster(args.output_fn, complement, gdir, annotation_dir, ommit_list, limit_max,
                                                      [["E78", 1], ["H28", 0], ["H59", 0], ["E20", 0]], "channel2b", output=True)
    output(cluster2, args.dir, "channel2b", "-")

    print("channel2b")
    print("lens - ok ", len(cluster), len(complement), len(missing))
    print("lens - ok ", len(cluster2), len(complement2), len(missing2))
    print("lens - all, missing", len(missing))


