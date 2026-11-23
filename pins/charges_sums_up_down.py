#!/usr/bin/env python3
# tool which generate for given pin directory its standart graphs with up/down part charges sum

import argparse
import time
import os
import json
import svgwrite
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def up_part(ch_data, part1, part2, text):

    up = []
    down = []
    res_min = 10000
    res_max = 0
    for res in ch_data:
        if int(res) < res_min:
            res_min = int(res)
        if int(res) > res_max:
            res_max = int(res)
    res_middle = res_min + (res_max - res_min) / 2
    aux = []
    res = res_min
    while res < res_middle:
        aux.append(ch_data[str(res)])
        res += 1

    if part1 == "up":
        up = aux
    else:
        down = aux
    aux2 = []
    while res <= res_max:
        aux2.append(ch_data[str(res)])
        res += 1
    if part2 == "up":
        up += aux2
    else:
        down += aux2

    return up, down


def up_part2(name, data, sses_list, data_all_res, min_index, max_index, part, text):
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = 1
        for i in data[sses_list[min_index]][name]:
            if res_min is None or res_min < int(i):
                res_min = int(i)
        res_min += 1

        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max is None or res_max > int(i):
                res_max = int(i)

        ch_data = {}
        res = res_min
        while res < res_max:
            ch_data[str(res)] = data_all_res[res]
            res += 1
        return up_part(ch_data, part, part, text)
    else:
        return [], []


def up_part_start(name, data, sses_list, data_all_res, sindex, part, text):
    if name in data[sses_list[sindex]]:
        res_min = 1000
        for i in data[sses_list[sindex]][name]:
            if res_min is None or res_min > int(i):
                res_min = int(i)

        ch_data = {}
        res = 1
        while res < res_min:
            ch_data[str(res)] = data_all_res[res]
            res += 1
        return up_part(ch_data, part, part, text)
    else:
        return [], []


def up_part_end(name, data, sses_list, data_all_res, sindex, part, text):
    if name in data[sses_list[sindex]]:
        res_max = 1
        for i in data[sses_list[sindex]][name]:
            if res_max is None or res_max < int(i):
                res_max = int(i)

        ch_data = {}
        res = res_max+1
        while res in data_all_res:
            ch_data[str(res)] = data_all_res[res]
            res += 1
        return up_part(ch_data, part, part, text)
    else:
        return [], []


def dolu(x, y, res_list, sses, opacity, text="", mmax=0, avg=False, protein_name="", map_file="aux_map"):
    res_len = float(y[1]-y[0])/len(res_list)
    res_min = None
    res_max = None
    for i in res_list:
        if res_min is None or res_min > int(i):
            res_min = int(i)
        if res_max is None or res_max < int(i):
            res_max = int(i)

    for i in res_list:
        if res_list[i] < -40:
            color = "red"
            dwg.add(dwg.line(start=(x, y[0]+int((int(i)-res_min)*res_len)), end=(x, y[0]+int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
        elif res_list[i] > 40:
            color = "blue"
            dwg.add(dwg.line(start=(x, y[0]+int((int(i)-res_min)*res_len)), end=(x, y[0]+int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
        elif opacity == 1:
            color = "white"
            dwg.add(dwg.line(start=(x, y[0]+int((int(i)-res_min)*res_len)), end=(x, y[0]+int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
    dwg.add(dwg.text(f"{text}", insert=(x-3, mmax-3), color="white", font_size=4, style="text-anchor:center"))
    if avg == 2:
        dwg.add(dwg.text(res_min, insert=(x+3, y[0]+5), color="white", font_size=6, style="text-anchor:center"))
        dwg.add(dwg.text(res_max, insert=(x+3, y[1]), color="white", font_size=6, style="text-anchor:center"))
    n = int((x-20)/20*32.22)
    y1 = 15 + (y[0] - 10) * 1.7
    y2 = 15 + (y[1] - 10) * 1.7
    with open(map_file, 'a') as fp:
        fp.write(f'<area alt="" title="sse{text}" href="http://147.251.115.132/results?ID={protein_name}_7.2_4&residue_start={res_min}&residue_end={res_max}" shape="rect" coords="{25+n},{y1},{35+n},{y2}" />\n')


def nahoru(x, y, res_list, sses, opacity, text="", mmax=0, avg=False, protein_name="", map_file="aux_map"):
    res_len = float(y[1]-y[0])/len(res_list)
    res_min = None
    res_max = None
    for i in res_list:
        if res_min is None or res_min > int(i):
            res_min = int(i)
        if res_max is None or res_max < int(i):
            res_max = int(i)
    for i in res_list:
        if res_list[i] < -40:
            color = "red"
            dwg.add(dwg.line(start=(x, y[1]-int((int(i)-res_min)*res_len)), end=(x, y[1]-int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
        elif res_list[i] > 40:
            color = "blue"
            dwg.add(dwg.line(start=(x, y[1]-int((int(i)-res_min)*res_len)), end=(x, y[1]-int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
        elif opacity == 1:
            color = "white"
            dwg.add(dwg.line(start=(x, y[1]-int((int(i)-res_min)*res_len)), end=(x, y[1]-int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
    dwg.add(dwg.text(f"{text}", insert=(x-3, mmax-3), color="white", font_size=4, style="text-anchor:center"))
    if avg == 2:
        dwg.add(dwg.text(res_min, insert=(x+3, y[1]), color="white", font_size=6, style="text-anchor:center"))
        dwg.add(dwg.text(res_max, insert=(x+3, y[0]+5), color="white", font_size=6, style="text-anchor:center"))
    n = int((x-20)/20*32.22)
    y1 = 15+(y[0]-10)*1.7
    y2 = 15+(y[1]-10)*1.7
    with open(map_file, 'a') as fp:
        fp.write(f'<area alt="" title="sse{text}" href="http://147.251.115.132/results?ID={protein_name}_7.2_4&residue_start={res_min}&residue_end={res_max}" shape="rect" coords="{25+n},{y1},{35+n},{y2}" />\n')


def res_line(s, e, data_all_res, res_min, res_max, pom=""):
    r = res_min
    r_diff = res_max - res_min
    ss = s
    sdif = (e[0]-s[0], e[1]-s[1])
    while r < res_max:
          if data_all_res[r] < -40:
              color = "red"
          elif data_all_res[r] > 40:
              color = "blue"
          else:
              color = "white"
          dwg.add(dwg.line(start=(ss[0]+sdif[0]/r_diff*(r-res_min), ss[1]+sdif[1]/r_diff*(r-res_min)), end=(ss[0]+sdif[0]/r_diff*(r+1-res_min), ss[1]+sdif[1]/r_diff*(r+1-res_min)), stroke_width=1, stroke=color))
          r += 1


def connection(name, data, sses_list, data_all_res, min_index, max_index, start_coords, center_coords, end_coords):
    # draw a connection between start_coords and end_coords crossig start_coords point
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = 1
        for i in data[sses_list[min_index]][name]:
            if res_min is None or res_min < int(i):
                res_min = int(i)
        res_min += 1
        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max is None or res_max > int(i):
                res_max = int(i)

        if res_min + 1 == res_max:  # only one residuum, whole loop has one color
            res_line(start_coords, center_coords, data_all_res, res_min, res_max)
            res_line(center_coords, end_coords, data_all_res, res_min, res_max)
        elif (res_max - res_min) % 2 == 1:  # loop has even number of residue, the first part contains +1 res then second
            res_line(start_coords, center_coords, data_all_res, res_min, res_min + (res_max - res_min+1) / 2)
            res_line(center_coords, end_coords, data_all_res, int(res_min + (res_max - res_min+1) / 2), res_max)
        else:  # loop has odd residues - both parts have the same number
            res_line(start_coords, center_coords, data_all_res, res_min, res_min + (res_max - res_min) / 2)
            res_line(center_coords, end_coords, data_all_res, int(res_min + (res_max - res_min) / 2), res_max)
    else:
        # there are no start or end data - draw the lines only simple
        dwg.add(dwg.line(start=start_coords, end=center_coords, stroke_width=1, stroke="black"))
        dwg.add(dwg.line(start=center_coords, end=end_coords, stroke_width=1, stroke="black"))


def simple_connection(name, data, sses_list, data_all_res, min_index, max_index, start_coords, end_coords):
    # draw a connection between start_coords and end_coords crossig start_coords point
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = 1
        for i in data[sses_list[min_index]][name]:
            if res_min is None or res_min < int(i):
                res_min = int(i)
        res_min += 1

        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max is None or res_max > int(i):
                res_max = int(i)
        res_line(start_coords, end_coords, data_all_res, res_min, res_max, pom="nnn")
    else:
        dwg.add(dwg.line(start=start_coords, end=end_coords, stroke_width=1, stroke="black"))


def graph_full(data, aux_data, sses_order, fn_order, avg, first, data_all_res, aname, map_file):
    # generate diagram and record
    opacity = 1/avg*2
    c_sses = 0
    c_fn = 0
    cc = 80
    for sses in sses_order:
         c_fn = 0
         c_sses += 1
         for name in fn_order:
              if first:
                  mmin = 10+c_fn*0
                  mmax = 60+c_fn*0
                  dwg.add(dwg.line(start=(100-cc, mmax), end=(100-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(120-cc, mmax), end=(120-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(160-cc, mmax), end=(160-cc, mmax-20), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(160-cc, mmin), end=(160-cc, mmin+20), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(140-cc, mmax), end=(140-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(180-cc, mmax), end=(180-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(200-cc, mmax), end=(200-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(220-cc, mmax), end=(220-cc, mmax-20), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(220-cc, mmin), end=(220-cc, mmin+20), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(240-cc, mmax), end=(240-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(260-cc, mmax), end=(260-cc, mmin), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(280-cc, mmax), end=(280-cc, mmax-22), stroke_width=4, stroke="#C0C0C0"))
                  dwg.add(dwg.line(start=(280-cc, mmin), end=(280-cc, mmin+22), stroke_width=4, stroke="#C0C0C0"))

    mmin = 10
    mmax = 60
    sum_up = []
    sum_down = []
    for name in fn_order:
          with open(map_file, 'a') as fp:
              fp.write(f'<map name="{name}" id="Map">\n')
          sname = name.split("-")
          protein_name = sname[1]
          if not avg:
              mmin = 10 + c_fn * 70
              mmax = 60 + c_fn * 70
          c_fn += 1
#          sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"]
          sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"]   # BIG DATASET
#          sses_list = ["H0", "H1", "H3", "H4", "H6", "H7", "H15", "H16", "H17", "H18", "H19", "H21", "H22"]   # arabidopsis
          if name in data[sses_list[0]]:
                   r_up, r_down = up_part(data[sses_list[0]][name], "up", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   dolu(100-cc, (mmin, mmax), data[sses_list[0]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)               # 1

          if name in data[sses_list[1]]:
                   r_up, r_down = up_part(data[sses_list[1]][name], "down", "up", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(120-cc, (mmin, mmax), data[sses_list[1]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)             # 2
          if name in data[sses_list[2]]:
                   r_up, r_down = up_part(data[sses_list[2]][name], "up", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   dolu(140-cc, (mmin, mmax), data[sses_list[2]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)               # 3
          if name in data[sses_list[4]]:
                   r_up, r_down = up_part(data[sses_list[4]][name], "up", "up", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(160-cc, (mmin, mmin+20), data[sses_list[4]][name], sses, opacity, f"{s_up}", mmin, avg, protein_name, map_file)          # 4
          if name in data[sses_list[3]]:
                   r_up, r_down = up_part(data[sses_list[3]][name], "down", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(160-cc, (mmax-20, mmax), data[sses_list[3]][name], sses, opacity, f".... {s_down}", mmin, avg, protein_name, map_file)          # 4
          if name in data[sses_list[5]]:
                   r_up, r_down = up_part(data[sses_list[5]][name], "up", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   dolu(180-cc, (mmin, mmax), data[sses_list[5]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)               # 5
          if name in data[sses_list[7]]:
                   r_up, r_down = up_part(data[sses_list[7]][name], "up", "up", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(280-cc, (mmin, mmin+22), data[sses_list[7]][name], sses, opacity, f"{s_up} ", mmin, avg, protein_name, map_file)          # 6
          if name in data[sses_list[6]]:
                   r_up, r_down = up_part(data[sses_list[6]][name], "down", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(280-cc, (mmax-22, mmax), data[sses_list[6]][name], sses, opacity, f".... {s_down}", mmin, avg, protein_name, map_file)          # 6
          if name in data[sses_list[8]]:
                   r_up, r_down = up_part(data[sses_list[8]][name], "up", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   dolu(260-cc, (mmin, mmax), data[sses_list[8]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)                  # 7
          if name in data[sses_list[9]]:
                   r_up, r_down = up_part(data[sses_list[9]][name], "down", "up", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(240-cc, (mmin, mmax), data[sses_list[9]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)             # 8
          if name in data[sses_list[11]]:
                   r_up, r_down = up_part(data[sses_list[11]][name], "down", "down", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   dolu(220-cc, (mmax-20, mmax), data[sses_list[11]][name], sses, opacity, f".... {s_down}", mmin, avg, protein_name, map_file)           # 9
          if name in data[sses_list[10]]:
                   r_up, r_down = up_part(data[sses_list[10]][name], "up", "up", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   dolu(220-cc, (mmin, mmin+20), data[sses_list[10]][name], sses, opacity, f"{s_up}", mmin, avg, protein_name, map_file)           # 9
          if name in data[sses_list[12]]:
                   r_up, r_down = up_part(data[sses_list[12]][name], "down", "up", name)
                   sum_up += r_up
                   sum_down += r_down
                   s_up = int(sum(r_up))
                   s_down = int(sum(r_down))
                   nahoru(200-cc, (mmin, mmax), data[sses_list[12]][name], sses, opacity, f"{s_up} {s_down}", mmin, avg, protein_name, map_file)            # 10
          with open(map_file, 'a') as fp:
              fp.write("</map>\n")

    pom_up = sum(sum_up)
    pom_down = sum(sum_down)
    # in the order
    r_up, r_down = up_part_start(name, data, sses_list, data_all_res, 0, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 0, 1, (100-cc, mmax), (110-cc, mmax+4), (120-cc, mmax))   # 1 -> 2
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 0, 1, "down", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 1, 2, (120-cc, mmin), (130-cc, mmin-4), (140-cc, mmin))   # 2 -> 3
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 1, 2, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 2, 3, (140-cc, mmax), (150-cc, mmax+4), (160-cc, mmax))   # 3 -> 4
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 2, 3, "down", name)
    sum_up += r_up
    sum_down += r_down

    simple_connection(name, data, sses_list, data_all_res, 3, 4, (160-cc, mmin+20), (160-cc, mmax-20))  # 6 -> 6
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 3, 4, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 4, 5, (160-cc, mmin), (170-cc, mmin-4), (180-cc, mmin))   # 4 -> 5
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 4, 5, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 5, 6, (180-cc, mmax), (230-cc, mmax+12), (280-cc, mmax))  # 5 -> 6
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 5, 6, "down", name)
    sum_up += r_up
    sum_down += r_down

    simple_connection(name, data, sses_list, data_all_res, 6, 7, (280-cc, mmax-22), (280-cc, mmin+22))  # 6 -> 6
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 6, 7, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 7, 8, (260-cc, mmin), (270-cc, mmin-4), (280-cc, mmin))   # 6 -> 7
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 7, 8, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 8, 9, (240-cc, mmax), (250-cc, mmax+4), (260-cc, mmax))   # 7 -> 8
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 8, 9, "down", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 9, 10, (220-cc, mmin), (230-cc, mmin-4), (240-cc, mmin))  # 8 -> 9
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 9, 10, "up", name)
    sum_up += r_up
    sum_down += r_down

    simple_connection(name, data, sses_list, data_all_res, 10, 11, (220-cc, mmin+20), (220-cc, mmax-20))  # 6 -> 6
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 10, 11, "up", name)
    sum_up += r_up
    sum_down += r_down

    connection(name, data, sses_list, data_all_res, 11, 12, (200-cc, mmax), (210-cc, mmax+4), (220-cc, mmax))  # 9 -> 10
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 11, 12, "down", name)
    sum_up += r_up
    sum_down += r_down

    r_up, r_down = up_part_end(name, data, sses_list, data_all_res, 12, "up", name)
    sum_up += r_up
    sum_down += r_down

    # pom_up, pom_down is up/down parts of all parts
    s_up = sum(sum_up)
    s_down = sum(sum_down)

    # pom_up, pom_down is up/down parts of all loops  + start end only
    pom_up = sum(sum_up) - pom_up
    pom_down = sum(sum_down) - pom_down
    r_up, r_down = up_part2(name, data, sses_list, data_all_res, 5, 6, "down", name)

    dwg.add(dwg.text(f"({int(s_up)},{int(s_down)})", insert=(10, 70), color="white", font_size=8, style="text-anchor:center"))
    dwg.add(dwg.text(f"<{int(pom_up)},{int(pom_down)}>", insert=(10, 50), color="white", font_size=8, style="text-anchor:center"))
    # printloops only charges
    return s_up, s_down, pom_up, pom_down, sum(r_up), sum(r_down)


# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument('--charges_dir', required=True, action="store",
                    help="charges data dir")
parser.add_argument('--org_name_dir', required=False, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument("--name", help="name")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = os.path.basename(os.path.normpath(args.directory))
new_dir = f"charges-{name}-{timestamp}"
map_file = f'{new_dir}/map_file.html'
os.makedirs(new_dir, exist_ok=True)

# Creates a new map file
with open(map_file, 'w') as fp:
    pass

print("new_dir", new_dir)

dall = {}
d2all = []
dup = {}
ddown = {}
for ff in os.listdir(args.directory):
  dall[ff] = []
  dup[ff] = []
  ddown[ff] = []
  data = {}      # charges data
  aux_data = {}  # organism name
  fn_order = []
  data_all_res = {}
  for filename in os.listdir(f"{args.directory}/{ff}"):
    # go through all PINs in PIN directory and create html record and diagram for them

    # read PIN json and read residues charges in SSEs
    with open(f"{args.directory}/{ff}/{filename}", "r") as fd:
        jdata = json.load(fd)
        fn_order.append(filename)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in data:
                    data[sse] = {}
                data[sse][filename] = jdata["charges"][sse]

    # and (if org_name_dir is set) organism name as well
    if args.org_name_dir:
        sfilename = filename.split("-")
        with open(f"{args.org_name_dir}/{sfilename[1]}-F1.cif") as fd:
            line = fd.readline()
            while line:
#                if line.startswith("_ma_target_ref_db_details.ncbi_taxonomy_id"):
#                    pom = int(line[55:])
                if line.startswith("_ma_target_ref_db_details.organism_scientific"):
                    aux_data[filename] = line[55:]
                    pom = line[55:]
                    break
                line = fd.readline()
        aux_data[filename] = pom
    else:
#        continue
        aux_data[filename] = "--"
#    if pom == 3702:
#    if "Arabidop" in pom:
#        aux_data[filename] = pom
#        aux_data[filename] = "Arabidopsis thaliana"
#    else:
#        continue

    # read PIN json and read residues charges in parts outside SSEs
    sfilename = filename.split("-")
    ac_filename = f"{sfilename[1]}.alphacharges"
    with open(f"{args.charges_dir}/{ac_filename}", "rb") as fd:
        ac_res = {}
        r = fd.read()

        line = ""
        start = False

        for c in r:
            if chr(c) != "\n":
                # read separate character and put them to to lines
                line = line + chr(c)
            else:
                if "ATOM" in line:
                    # the first line contains the rest of old stuff at the beginning, need to be moved
                    pos = line.find("ATOM")
                    line = line[pos:]
                if line.startswith("ATOM") and start:
                    if int(line[21:26]) not in ac_res:
                        ac_res[int(line[21:26])] = []
                    ac_res[int(line[21:26])].append(float(line[55:62]))

                line = ""
            if line.startswith("TER"):
                # read only stuff between "TER" signs
                if not start:
                    start = True
                else:
                    break

    # data aggregation
    data_all_res[filename] = {}
    for r in ac_res:
        s = 0
        for ss in ac_res[r]:
            s += ss
        data_all_res[filename][r] = s*100
    for sse in data:
        if filename in data[sse]:
            for res in data[sse][filename]:
                data[sse][filename][res] = data_all_res[filename][int(res)]

  aux_data_order = sorted(aux_data.items(), key=lambda x: x[1])
  fn_order = []

  for i in aux_data_order:
      fn_order.append(i[0])
  sses_order = list(data.keys())

  for f in fn_order:

    # read filename and diagram size
    fnsplit = f.split("-")
    prot_name = fnsplit[1]
    if len(aux_data[f]) > 21:
        org_name = aux_data[f][1:20] + "."
    else:
        org_name = aux_data[f][1:-2]
    c = 0
    fsplit = f.split("-")
    name = fsplit[1]
    width = 220
    height = 75
    dwg = svgwrite.Drawing(f"{new_dir}/charges-{prot_name}-full-{prot_name}.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))

    # generate diagram and recors
    u1, d1, u2, d2, r1, r2 = graph_full(data, aux_data, sses_order, [f], 2, True, data_all_res[f], args.name, map_file)
    #        PINx      B09I2E8  charge of whole protein          up part        down part       only loops up   only loops down
    print(f"{ff}	{name}	{int(round((u1+d1)/100))}	{int(u1)}	{int(d1)}	{int(u2)}	{int(d2)}	{f}")
    dall[ff].append(int(round((u1+d1)/100)))
    d2all.append({"":ff, "average charge":int(round((u1+d1)/100)), "average charge of non central loop side": u1/100, "average charge of central loop side":d1/100, "big loop":r2/100})
    dup[ff].append(u1/100)
    ddown[ff].append(d1/100)

    dwg.save()
    c += 1

df = pd.DataFrame(d2all)
sns.violinplot(data=df, y="average charge", x="", cut=0, bw=.2)
plt.show()

sns.violinplot(data=df, y="average charge of non central loop side", x="", cut=0, bw=.2)
plt.show()

sns.violinplot(data=df, y="average charge of central loop side", x="", cut=0, bw=.2)
plt.show()

sns.violinplot(data=df, y="big loop", x="", cut=0, bw=.2)
plt.show()
