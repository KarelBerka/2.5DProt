#!/usr/bin/env python3
# go through unannotated channels and write which channel is the nearest one

import argparse
import json
import os
import numpy
import time
import math


def pnt2line(pnt, start, end):
    m = math.sqrt((pnt[0]-start[0])*(pnt[0]-start[0])+(pnt[1]-start[1])*(pnt[1]-start[1])+(pnt[2]-start[2])*(pnt[2]-start[2]))
    n = math.sqrt((pnt[0]-end[0])*(pnt[0]-end[0])+(pnt[1]-end[1])*(pnt[1]-end[1])+(pnt[2]-end[2])*(pnt[2]-end[2]))
    if m > n:
        return n
    else:
        return m
    a = numpy.array(start)
    b = numpy.array(end)
    p = numpy.array(pnt)
    # normalized tangent vector
    d = numpy.divide(b - a, numpy.linalg.norm(b - a))
    # signed parallel distance components
    s = numpy.dot(a - p, d)
    t = numpy.dot(p - b, d)
    # clamped parallel distance
    h = numpy.maximum.reduce([s, t, 0])
    # perpendicular distance component
    c = numpy.cross(p - a, d)
    return numpy.hypot(h, numpy.linalg.norm(c))


def read_file_list(filename):
    fl = []
    with open(filename, "r") as fd:
        line = fd.readline()
        while line:
            sline = line.split(" ")
            if len(sline) == 2:
                fl.append([sline[0], sline[1][:-1]])
            else:
                fl.append([sline[0], sline[1], sline[2][:-1]])
            line = fd.readline()
    return fl


def generate_gnuplot_splot_input(input_filename, channel_name, new_file_name):
    with open(new_file_name, "w") as fw:
        with open(input_filename, "r") as fr:
            json_data = json.load(fr)
            for ch in json_data["Annotations"]:
                if ch['Name'] == channel_name.replace("+", " "):
                    ch_id = ch['Id']
            for ch_type in json_data['Channels']:
                for ch in json_data['Channels'][ch_type]:
                    if ch['Id'] == ch_id:
                        for j in ch['Profile']:
                            fw.write(f"{j['X']} {j['Y']} {j['Z']} \n")
    print(f"'{new_file_name}' lt rgb 'green'", end=", ")


def generate_gnuplot_splot_input2(input_filename, ch_id, new_file_name, near, domain, file_list_file):

    if near:
        # we have to test whether the channel is near enough <10A to
        # be shown
        with open(input_filename, "r") as fr, open(near, "r") as fa:
            json_data_chan = json.load(fr)
            json_data_sses = json.load(fa)
            if len(json_data_chan.keys()) >= 2:
              for ch_type in json_data_chan['Channels']:
                for ch in json_data_chan['Channels'][ch_type]:
                    if ch_id == None or ch['Id'] == ch_id: 
                        for j in ch['Profile']:
                            is_near = False
                            for sses in json_data_sses[domain]["secondary_structure_elements"]:
                                dist = pnt2line([j['X'], j['Y'], j['Z']], sses["start_vector"], sses["end_vector"])
                                if dist < 10:
                                    is_near = True
                                    break
                            if is_near:
                                break
                        if is_near:
                            if True:
                                with open(f"{new_file_name}_{ch['Id']}", "w") as fw:
                                    for j in ch['Profile']:
                                        fw.write(f"{j['X']} {j['Y']} {j['Z']} \n")
                                with open(file_list_file, "a") as fa:
                                    fa.write(f"'{new_file_name}_{ch['Id']}' lt rgb '{ch['Type']}'\n")
        return

    print("bb", input_filename, ch_id)
    with open(input_filename, "r") as fr:
        json_data = json.load(fr)
        print("jjj", json_data.keys())
        if len(json_data.keys()) >= 2:

            for ch_type in json_data['Channels']:
                print("ch_type", ch_type)
                for ch in json_data['Channels'][ch_type]:
                    print("ch", ch)
                    if ch_id == None or ch['Id'] == ch_id:
                        with open(f"{new_file_name}_{ch['Id']}", "w") as fw:
                            for j in ch['Profile']:
                                fw.write(f"{j['X']} {j['Y']} {j['Z']} \n")
                        print(f"'{new_file_name}_{ch['Id']}' lt rgb 'orange'")#, end=", ")


# read arguments
parser = argparse.ArgumentParser(
    description='generate gnuplot input data for input channel annotation files from input directory')
parser.add_argument('--dir', required=True, action="store",
                    help="channel files directory")
parser.add_argument('--unknown_file_list', action="store",
                    help="file containing protein names and chain numbers of hanis which would be shown")
parser.add_argument('--named_file_list', action="store",
                    help="file containing protein names and chain numbers of hanis which would be shown")
parser.add_argument('--near', action="store_true",
                    help="we want to have only channels which are near to sses of the given domain <10A")
parser.add_argument('--all', action="store_true",
                    help="we want to have only channels which are near to sses of the given domain <10A")
parser.add_argument('--af', action="store_true",
                    help="direct mole channel files (mole output dir)")
parser.add_argument('--name', required=True, action="store",
                    help="name of the set - probably family name")
parser.add_argument('--file_list', required=True, action="store",
                    help="set the location of output file list")
parser.add_argument('--whole_list', required=True, action="store",
                    help="set the location of gnuplot readable whole output list")
parser.add_argument('--annot_dir', action="store",
                    help="annotation files directory")


args = parser.parse_args()

timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
new_dir = f"gnuplot-{args.name}-{timestamp}"
os.mkdir(new_dir)

if args.named_file_list:
    file_list = read_file_list(args.named_file_list)
    for i in file_list:
        threed_layout_filename = f"{args.dir}/{i[0]}-channel.json"
        generate_gnuplot_splot_input(threed_layout_filename, i[2], f"gnuplot-{i[0][:4]}{i[0][6:]}{i[0][4:6]}-{i[1]}.txt")
    print()

if args.unknown_file_list:
    file_list = read_file_list(args.unknown_file_list)
    for i in file_list:
        if len(i) > 2:
            output_filename = f"gnuplot-{i[0][:4]}{i[0][6:]}{i[0][4:6]}-{i[1]}-{i[2]}.txt"
        else:
            output_filename = f"gnuplot-{i[0][:4]}{i[0][6:]}{i[0][4:6]}-{i[1]}.txt"
        if args.near:
            near = f"{args.dir}/{i[0][:4]}{i[0][6:]}{i[0][4:6]}-annotated.sses.json"
        else:
            near = False
        threed_layout_filename = f"{args.dir}/{i[0]}-channel.json"
        generate_gnuplot_splot_input2(threed_layout_filename, i[1], output_filename, near, f"{i[0][:4]}{i[0][6:]}{i[0][4:6]}")
    print()

if args.all:
    print("all")
    for file_name in os.listdir(args.dir):
         if file_name.endswith("-channel.json"):
             print("filename", file_name)
             sfn = file_name.split("-")
             pdc = sfn[0]
             if args.near:
                 near = f"{args.dir}/{pdc[0:4]}{pdc[6]}{pdc[4:6]}-annotated.sses.json"
             else:
                 near = False
             generate_gnuplot_splot_input2(f"{args.dir}/{file_name}",
                                           None,
                                           f"{new_dir}/gnuplot-{sfn[0]}.txt",
                                           near,
                                           f"{pdc[0:4]}{pdc[6]}{pdc[4:6]}",
                                           args.file_list)
    with open(args.whole_list, "w") as fw:
        fw.write(f"{args.file_list} None 0.5,0.5,1\n")


if args.af:
    for af_dir in os.listdir(args.dir):
        file_name=f"{args.dir}/{af_dir}/json/data.json"
        if os.path.isfile(file_name):
             sfn = file_name.split("-")
             pdc = sfn[0]
             near = False

             generate_gnuplot_splot_input2(f"{file_name}",
                                           None,
                                           f"{new_dir}/gnuplot-{af_dir}.txt",
                                           near,
                                           f"{pdc[0:4]}{pdc[6]}{pdc[4:6]}",
                                           args.file_list)
    with open(args.whole_list, "w") as fw:
        fw.write(f"{args.file_list} None 0.5,0.5,1\n")

