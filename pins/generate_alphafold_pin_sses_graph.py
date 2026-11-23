#!/usr/bin/env python3
# tool which generate for given alphafold family
# its diagrams

import argparse
import time
import os
import json
import math
import svgwrite


def dolu(x, y, res_list, sses, opacity, text="", mmax=0, avg=False, protein_name=""):
    res_len = float(y[1]-y[0])/len(res_list)
    res_min = None
    res_max = None
    for i in res_list:
        if res_min == None or res_min > int(i):
            res_min = int(i)
        if res_max == None or res_max < int(i):
            res_max = int(i)

    for i in res_list:
        if res_list[i] < -40:
            color = "red"
            dwg.add(dwg.line(start=(x,y[0]+int((int(i)-res_min)*res_len)), end=(x,y[0]+int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
        elif res_list[i] > 40:
            color = "blue"
            dwg.add(dwg.line(start=(x,y[0]+int((int(i)-res_min)*res_len)), end=(x,y[0]+int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
        elif opacity == 1:
            color = "white"
            dwg.add(dwg.line(start=(x,y[0]+int((int(i)-res_min)*res_len)), end=(x,y[0]+int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
    dwg.add(dwg.text(f"{text}", insert=(x-3, mmax-3), color="white", font_size=8, style="text-anchor:center"))
    if avg == 2:
        dwg.add(dwg.text(res_min, insert=(x+3, y[0]+5), color="white", font_size=6, style="text-anchor:center"))
#        dwg.add(dwg.text(f"({res_max-res_min})", insert=(x+3, y[0]+12), color="white", font_size=5, style="text-anchor:center"))
        dwg.add(dwg.text(res_max, insert=(x+3, y[1]), color="white", font_size=6, style="text-anchor:center"))
    n = int((x-20)/20*32.22)
    y1 = 15+(y[0]-10)*1.7
    y2 = 15+(y[1]-10)*1.7
#    print("min", y[0], "max", y[1])
    print(f'<area alt="" title="sse{text}" href="http://147.251.115.132/results?ID={protein_name}_7.2_4&residue_start={res_min}&residue_end={res_max}" shape="rect" coords="{25+n},{y1},{35+n},{y2}" />')


def nahoru(x,y,res_list, sses, opacity, text="", mmax=0, avg=False, protein_name=""):
    res_len = float(y[1]-y[0])/len(res_list)
    res_min = None
    res_max = None
    for i in res_list:
        if res_min == None or res_min > int(i):
            res_min = int(i)
        if res_max == None or res_max < int(i):
            res_max = int(i)
    for i in res_list:
        if res_list[i] < -40:
            color = "red"
            dwg.add(dwg.line(start=(x,y[1]-int((int(i)-res_min)*res_len)), end=(x,y[1]-int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color,  style=f"opacity:{opacity}"))
        elif res_list[i] > 40:
            color = "blue"
            dwg.add(dwg.line(start=(x,y[1]-int((int(i)-res_min)*res_len)), end=(x,y[1]-int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color,  style=f"opacity:{opacity}"))
        elif opacity == 1:
            color = "white"
            dwg.add(dwg.line(start=(x,y[1]-int((int(i)-res_min)*res_len)), end=(x,y[1]-int((int(i)-res_min+1)*res_len)), stroke_width=4, stroke=color,  style=f"opacity:{opacity}"))
    dwg.add(dwg.text(f"{text}", insert=(x-3, mmax-3), color="white", font_size=8, style="text-anchor:center"))
    if avg == 2:
        dwg.add(dwg.text(res_min, insert=(x+3, y[1]), color="white", font_size=6, style="text-anchor:center"))
#        dwg.add(dwg.text(f"({res_max-res_min})", insert=(x+3, y[0]+12), color="white", font_size=5, style="text-anchor:center"))
        dwg.add(dwg.text(res_max, insert=(x+3, y[0]+5), color="white", font_size=6, style="text-anchor:center"))
    n = int((x-20)/20*32.22)
    y1 = 15+(y[0]-10)*1.7
    y2 = 15+(y[1]-10)*1.7
#    print("min", y[0], "max", y[1])
    print(f'<area alt="" title="sse{text}" href="http://147.251.115.132/results?ID={protein_name}_7.2_4&residue_start={res_min}&residue_end={res_max}" shape="rect" coords="{25+n},{y1},{35+n},{y2}" />')


def res_line(s, e, data_all_res, res_min, res_max, opacity, pom = ""):
    r = res_min
    r_diff = res_max - res_min
    ss = s
    sdif = (e[0]-s[0], e[1]-s[1])
    while r < res_max:
#          print(f"{pom} krok - diff {r_diff} krok min {r-res_min} max{r+1-res_min}, e {e}, vypisu na {(ss[0]+sdif[0]/r_diff*(r+1-res_min))}" )
       if data_all_res != 0:
          if data_all_res[r] < -40:
              color = "red"
          elif data_all_res[r] > 40:
              color = "blue"
          else: #if opacity == 1:
              color = "white"
#          print("<", color,r, data_all_res,  opacity, ">", end=": ")
          dwg.add(dwg.line(start=(ss[0]+sdif[0]/r_diff*(r-res_min), ss[1]+sdif[1]/r_diff*(r-res_min)), end=(ss[0]+sdif[0]/r_diff*(r+1-res_min), ss[1]+sdif[1]/r_diff*(r+1-res_min)), stroke_width=1, stroke=color, style=f"opacity:{opacity}"))

#       else:
#           return
#           color = "black"
#          print("mmm", pom, r, data_all_res[r], color)
##       print("<", color,r, opacity, ">", end=": ")
##       dwg.add(dwg.line(start=(ss[0]+sdif[0]/r_diff*(r-res_min), ss[1]+sdif[1]/r_diff*(r-res_min)), end=(ss[0]+sdif[0]/r_diff*(r+1-res_min), ss[1]+sdif[1]/r_diff*(r+1-res_min)), stroke_width=1, stroke=color, style=f"opacity:{opacity}"))
       r += 1

def connection(name, data, sses_list, data_all_res, min_index, max_index, start_coords, center_coords, end_coords, opacity):
    print("opacity", name, opacity)
    # draw a connection between start_coords and end_coords crossig start_coords point
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = 1
        for i in data[sses_list[min_index]][name]:
            if res_min == None or res_min < int(i):
                res_min = int(i)
        res_min += 1
        
        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max == None or res_max > int(i):
                res_max = int(i)
        
        if res_min + 1 == res_max: #only one residuum, whole loop has one color
            res_line(start_coords, center_coords, data_all_res, res_min, res_max, opacity)
            res_line(center_coords, end_coords, data_all_res, res_min, res_max, opacity)
        elif (res_max - res_min) % 2 == 1: # loop has even number of residue, the first part contains +1 res then second 
            res_line(start_coords, center_coords, data_all_res, res_min, res_min +(res_max - res_min+1) / 2, opacity)
            res_line(center_coords, end_coords, data_all_res, int(res_min +(res_max - res_min+1) / 2), res_max, opacity)
        else: # loop has odd residues - both parts have the same number
            res_line(start_coords, center_coords, data_all_res, res_min, res_min +(res_max - res_min) / 2, opacity)
            res_line(center_coords, end_coords, data_all_res, int(res_min +(res_max - res_min) / 2), res_max, opacity)
    else:
        # there are no start or end data - draw the lines only simple
        dwg.add(dwg.line(start=start_coords, end=center_coords, stroke_width=1, stroke="black"))
        dwg.add(dwg.line(start=center_coords, end=end_coords, stroke_width=1, stroke="black"))


def simple_connection(name, data, sses_list, data_all_res, min_index, max_index, start_coords, end_coords, opacity):
    # draw a connection between start_coords and end_coords crossig start_coords point
    print("opacity", name, opacity)
    if name in data[sses_list[min_index]] and name in data[sses_list[max_index]]:
        res_min = 1
        for i in data[sses_list[min_index]][name]:
            if res_min == None or res_min < int(i):
                res_min = int(i)
        res_min += 1

        res_max = None
        for i in data[sses_list[max_index]][name]:
            if res_max == None or res_max > int(i):
                res_max = int(i)
#        print("nnn", min_index, max_index,",", data[sses_list[min_index]][name], data[sses_list[max_index]][name], "---", res_min, res_max)
        res_line(start_coords, end_coords, data_all_res, res_min, res_max, opacity, pom="nnn")
    else:
        dwg.add(dwg.line(start=start_coords, end=end_coords, stroke_width=1, stroke="black"))


def graph_full(data, aux_data, sses_order, fn_order, avg, first, data_all_res):
    opacity = 1/math.sqrt(avg)*2
    c_sses = 0
    c_fn = 0
    cc = 80
    for sses in sses_order:
         c_fn = 0
         c_sses += 1
         for name in fn_order:
              if first:
                  mmin=10 +c_fn*0
                  mmax=60 +c_fn*0
                  color = "#C0C0C0"
                  color = "white"
                  
                  dwg.add(dwg.line(start=(100-cc,mmax), end=(100-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(120-cc,mmax), end=(120-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(160-cc,mmax), end=(160-cc, mmax-20), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(160-cc,mmin), end=(160-cc, mmin+20), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(140-cc,mmax), end=(140-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(180-cc,mmax), end=(180-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(200-cc,mmax), end=(200-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(220-cc,mmax), end=(220-cc, mmax-20), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(220-cc,mmin), end=(220-cc, mmin+20), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(240-cc,mmax), end=(240-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(260-cc,mmax), end=(260-cc, mmin), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(280-cc,mmax), end=(280-cc, mmax-22), stroke_width=4, stroke=color))
                  dwg.add(dwg.line(start=(280-cc,mmin), end=(280-cc, mmin+22), stroke_width=4, stroke=color))

    mmin = 10
    mmax = 60
#    print("================", fn_order)
    for name in fn_order:
          print(f'<map name="{name}" id="Map">')
          sname = name.split("-")
          protein_name = sname[1]
          if not avg:
              mmin = 10 + c_fn * 70
              mmax = 60 + c_fn * 70
          c_fn += 1
#          sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"] #standart 45 set
          sses_list = ["H6", "H7", "H9", "H11", "H13", "H15", "H44", "H45", "H48", "H51", "H53", "H56", "H59"] # BIG DATASET

          if name in data[sses_list[0]]:
                   dolu(100-cc, (mmin,mmax), data[sses_list[0]][name], sses, opacity, "1", mmin, avg, protein_name)               # 1
#                   print('<area alt="" title="sse1" href="http://147.251.115.132/results?ID=A0A0H3V2K9_7.2_4&residue_start=4&residue_end=30" shape="rect" coords="25,15,35,100" />')
          connection(name, data, sses_list, data_all_res, 0, 1, (100-cc, mmax), (110-cc, mmax+4), (120-cc, mmax), opacity)   # 1 -> 2
          if name in data[sses_list[1]]:
                   nahoru(120-cc, (mmin,mmax), data[sses_list[1]][name], sses, opacity, "2", mmin, avg, protein_name)             # 2
          connection(name, data, sses_list, data_all_res, 1, 2, (120-cc, mmin), (130-cc, mmin-4), (140-cc, mmin), opacity)   # 2 -> 3
          if name in data[sses_list[2]]:
                   dolu(140-cc, (mmin,mmax), data[sses_list[2]][name], sses, opacity, "3", mmin, avg, protein_name)               # 3
          connection(name, data, sses_list, data_all_res, 2, 3, (140-cc, mmax), (150-cc, mmax+4), (160-cc, mmax), opacity)   # 3 -> 4
          if name in data[sses_list[4]]:
                   nahoru(160-cc, (mmin,mmin+20), data[sses_list[4]][name], sses, opacity, "4", mmin, avg, protein_name)          # 4
          simple_connection(name, data, sses_list, data_all_res, 3, 4, (160-cc, mmin+20), (160-cc, mmax-20), opacity)  # , opacity6 -> 6
#         dwg.add(dwg.line(start=(160-cc,mmin+20), end=(160-cc, mmax-20), stroke_width=1, stroke="black"))          # 4 -> 4
          if name in data[sses_list[3]]:
                   nahoru(160-cc, (mmax-20,mmax), data[sses_list[3]][name], sses, opacity, "4", mmin, avg, protein_name)          # 4
          connection(name, data, sses_list, data_all_res, 4, 5, (160-cc, mmin), (170-cc, mmin-4), (180-cc, mmin), opacity)   # 4 -> 5
          if name in data[sses_list[5]]:
                   dolu(180-cc, (mmin,mmax), data[sses_list[5]][name], sses, opacity, "5", mmin, avg, protein_name)               # 5
          connection(name, data, sses_list, data_all_res, 5, 6, (180-cc, mmax), (230-cc, mmax+12), (280-cc, mmax), opacity)  # 5 -> 6
          if name in data[sses_list[7]]:
                   nahoru(280-cc, (mmin,mmin+22), data[sses_list[7]][name], sses, opacity, "6", mmin, avg, protein_name)          # 6
          simple_connection(name, data, sses_list, data_all_res, 6, 7, (280-cc, mmax-22), (280-cc, mmin+22), opacity)  # 6 -> 6
#          dwg.add(dwg.line(start=(280-cc,mmin+22), end=(280-cc, mmax-22), stroke_width=1, stroke="black"))          # 6 -> 6
          if name in data[sses_list[6]]:
                   nahoru(280-cc, (mmax-22,mmax), data[sses_list[6]][name], sses, opacity, "6", mmin, avg, protein_name)          # 6
          connection(name, data, sses_list, data_all_res, 7, 8, (280-cc, mmin), (270-cc, mmin-4), (260-cc, mmin), opacity)   # 6 -> 7
          if name in data[sses_list[8]]:
                   dolu(260-cc, (mmin,mmax), data[sses_list[8]][name], sses, opacity, "7", mmin, avg, protein_name)                  # 7
          connection(name, data, sses_list, data_all_res, 8, 9, (260-cc, mmax), (250-cc, mmax+4), (240-cc, mmax), opacity)   # 7 -> 8
          if name in data[sses_list[9]]:
                   nahoru(240-cc, (mmin,mmax), data[sses_list[9]][name], sses, opacity, "8", mmin, avg, protein_name)             # 8
          connection(name, data, sses_list, data_all_res, 9, 10, (240-cc, mmin), (230-cc, mmin-4), (220-cc, mmin), opacity)  # 8 -> 9
          if name in data[sses_list[11]]:
                   dolu(220-cc, (mmax-20,mmax), data[sses_list[11]][name], sses, opacity, "9", mmin, avg, protein_name)           # 9
          simple_connection(name, data, sses_list, data_all_res, 10, 11, (220-cc, mmin+20), (220-cc, mmax-20), opacity)  # 6 -> 6
#          dwg.add(dwg.line(start=(220-cc,mmin+20), end=(220-cc, mmax-20), stroke_width=1, stroke="black"))          # 9 -> 9
          if name in data[sses_list[10]]:
                   dolu(220-cc, (mmin,mmin+20), data[sses_list[10]][name], sses, opacity, "9", mmin, avg, protein_name)           # 9
          connection(name, data, sses_list, data_all_res, 11, 12, (220-cc, mmax), (210-cc, mmax+4), (200-cc, mmax), opacity) # 9 -> 10
          if name in data[sses_list[12]]:
                   nahoru(200-cc, (mmin,mmax), data[sses_list[12]][name], sses, opacity,"10", mmin, avg, protein_name)            # 10
          print("</map>")


def graph(data, aux_data, sses_order, fn_order, avg, label):
    opacity = 1/math.sqrt(avg)

    if avg:
        dwg.add(dwg.text("multiple transparent - ", insert=(10, 70-20), color="#C0C0C0", font_size=9, style="text-anchor:center"))
    else:
        for name in fn_order:
            if len(aux_data[name])>21:
                org_name = aux_data[name][1:20]+"."
            else:
                org_name = aux_data[name][1:-2]
    for sses in sses_order:
         for name in fn_order:
              if name in data[sses]:
                  res_min = None
                  for i in data[sses][name]:
                      if res_min == None or res_min > int(i):
                          res_min = int(i)
#                      if res_max == None or res_max < int(i):
#                          res_max = int(i)
#                      dwg.add(dwg.text(res_min, insert=(5+c_sses*100+ len(data[sses][name])*1+10-cc-11, c_fn*20-3), color="white", font_size=6, style="text-anchor:center"))
#                      dwg.add(dwg.text(f"{text}", insert=(5+c_sses*100+ len(data[sses][name])*1+10-cc, c_fn*20-3), color="white", font_size=8, style="text-anchor:center"))
#                      dwg.add(dwg.text(res_max, insert=(5+c_sses*100+ len(data[sses][name])*1+10-cc+10, c_fn*20-3), color="white", font_size=6, style="text-anchor:center"))
##                  if label in ["1","3", "5", "7"]:
##                      dwg.add(dwg.text("Nc", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
#                      dwg.add(dwg.text("Cy", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                      dwg.add(dwg.text("Cy", insert=(2+20*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                  elif label in ["2", "8", "10"]:
##                      if (int(i)-res_min) > 5:
##                          dwg.add(dwg.text("Cy", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
#                          dwg.add(dwg.text("Nc", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                          dwg.add(dwg.text("Nc", insert=(2+20*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                      else:
##                          dwg.add(dwg.text("Cy", insert=(0, 8), color="white", font_size=8, style="text-anchor:center"))
#                          dwg.add(dwg.text("Nc", insert=(3+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                          dwg.add(dwg.text("Nc", insert=(3+(20)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                  elif label in ["4a", "6a"]:
##                      dwg.add(dwg.text("Cy", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
##                  elif label in ["4b"]:
#                      dwg.add(dwg.text("Nc", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                      dwg.add(dwg.text("Nc", insert=(2+(8)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                  elif label in ["6b"]:
#                      dwg.add(dwg.text("Nc", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                      dwg.add(dwg.text("Nc", insert=(2+(15)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                  elif label in ["9a"]:
##                      dwg.add(dwg.text("Nc", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
##                  elif label in ["9b"]:
#                      dwg.add(dwg.text("Cy", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
##                      dwg.add(dwg.text("Cy", insert=(2+(8)*2, 8), color="white", font_size=8, style="text-anchor:center"))
                  for i in data[sses][name]:
                      if data[sses][name][i] < -40:
                          color = "red"
                      elif data[sses][name][i] > 40:
                          color = "blue"
                      else: color = "white"
                      dwg.add(dwg.line(start=(5+(int(i)-res_min)*2, 12), end=(5+(int(i)-res_min+1)*2, 12), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))


# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument('--charges_dir', required=True, action="store",
                    help="charges data dir")
parser.add_argument('--org_name_dir', required=False, action="store",
                    help="name of the directory, which contains domain charges layout")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = os.path.basename(os.path.normpath(args.directory))

new_dir = f"charges-{name}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)

print("new_dir", new_dir)


data = {}      # charges data
aux_data = {}  # organism name
fn_order = []
data_all_res = {}
for filename in os.listdir(args.directory):
    print("filename", filename, "dir", args.directory)
    
    with open(f"{args.directory}/{filename}", "r") as fd:
        jdata = json.load(fd)
        fn_order.append(filename)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in data:
                    data[sse] = {}
                data[sse][filename] = jdata["charges"][sse]
                if args.org_name_dir:
                     sfilename=filename.split("-")
                     with open(f"{args.org_name_dir}/{sfilename[1]}-F1.cif") as fd:
                         line = fd.readline()
                         while line:
                             if line.startswith("_ma_target_ref_db_details.organism_scientific"):
                                  aux_data[filename] = line[55:]
                                  break
                             line = fd.readline()
                else:
                     aux_data[filename] = "--"
    sfilename = filename.split("-")
    ac_filename=f"{sfilename[1]}.alphacharges"

    
    with open(f"{args.charges_dir}/{ac_filename}", "rb") as fd:
        ac_res={}
        counter = 0
        line = fd.readline()
        start = False
        data_all_res[filename] = {}
        while line:
            counter += 1
            strline = str(line)[2:-2]
            if strline.startswith("TER"):
                if start == False:
                    start = True
                else:
                    break
            if start and strline.startswith("ATOM") and len(line)<80:
                if int(line[21:26]) not in ac_res:
                     ac_res[int(line[21:26])]=[]
                ac_res[int(line[21:26])].append(float(line[55:62]))
            line = fd.readline()

    for r in ac_res:
         s = 0
         for ss in ac_res[r]:
              s += ss
         data_all_res[filename][r] = int(s*100)

aux_data_order = sorted(aux_data.items(), key=lambda x:x[1])
fn_order = []

for i in aux_data_order:
    fn_order.append(i[0])
sses_order = list(data.keys())


print("=========================")


if False:
    # simple
    print('<tr><td></td><td></td><td colspan="2" bgcolor="#ADD8E6">Scaffold domain</td><td colspan="4" bgcolor="#ADD8E6">Transporter domain - left</td><td colspan="3"  bgcolor="#90EE90">Scaffold domain</td><td colspan="4" bgcolor="#90EE90">Transporter domain - right</td></tr>')
    print("<tr><td></td><td></td><td>1</td><td>2</td><td>3</td><td>4a</td><td>4b</td><td>5</td><td>6a</td><td>6b</td><td>7</td><td>8</td><td>9a</td><td>9b</td><td>10</td></tr>")
    for f in fn_order:
        fnsplit = f.split("-")
        prot_name = fnsplit[1]
        if len(aux_data[f])>21:
            org_name = aux_data[f][1:20]+"."
        else:
            org_name = aux_data[f][1:-2]
        print(f'<tr><td>{args.directory}</td><td>{prot_name} - {org_name}<br> <a href="https://alphacharges.ncbr.muni.cz/results?ID={prot_name}_7.2_4">AlphaCharges</a>, <a href="https://alphafold.ebi.ac.uk/entry/{prot_name}">AlphaFold</a></td>')
        c = 0
        for sse in sses_order:
            fsplit = f.split("-")
            name = fsplit[1]
            if f not in data[sse]:
                print("<td></td>")
                c+=1
                continue
            width = 2*(len(data[sse][f]))+10
            height = 20
            text = ["1", "2", "3", "4a", "4b", "5", "6a", "6b", "7", "8", "9a", "9b", "10"]
            dwg = svgwrite.Drawing(f"{new_dir}/charges-{prot_name}-{text[c]}.svg", size=(width, height))
            dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
            graph(data, aux_data, [sse], [f], 2, text[c])
            dwg.save()
            res_min = None
            res_max = None
            for i in data[sse][f]:
                if res_min == None or res_min > int(i):
                    res_min = int(i)
                if res_max == None or res_max < int(i):
                    res_max = int(i)
            print(f'<td bgcolor="#808080"><img src="{args.directory}/charges-{prot_name}-{text[c]}.svg" height=30 usemap="#{f}_{c}"></td>')
            print(f'<map name="{f}_{c}" id="Map"><area alt="" title="sse{len(data[sse][f])}" href="http://147.251.115.132/results?ID={prot_name}_7.2_4&residue_start={res_min}&residue_end={res_max}" shape="rect" coords="10,18,{10+len(data[sse][f])*5},30" /></map></td>')

            c+=1
        print("</tr>")

if True:
    # avg simple
    print("<tr><td></td><td></td><td>1</td><td>2</td><td>3</td><td>4a</td><td>4b</td><td>5</td><td>6a</td><td>6b</td><td>7</td><td>8</td><td>9a</td><td>9b</td><td>10</td></tr>")
    c = 0
    print(F"<td>{args.directory}</td><td>multi</td>", end=" ")
    for sse in sses_order:
        mmax = 0
        for f in fn_order:
            if f in data[sse] and len(data[sse][f]) > mmax:
                mmax = len(data[sse][f])
        width = 2*(mmax)+10
        height = 20
        text = ["1", "2", "3", "4a", "4b", "5", "6a", "6b", "7", "8", "9a", "9b", "10"]
        dwg = svgwrite.Drawing(f"{new_dir}/charges-avg-{text[c]}.svg", size=(width, height))
        dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
        for f in fn_order:
            fsplit = f.split("-")
            name = fsplit[1]
            if f not in data[sse]:
                continue
            graph(data, aux_data, [sse], [f], len(fn_order), text[c])
            dwg.save()
        print(f'<td bgcolor="#808080"><img src="{args.directory}/charges-avg-{text[c]}.svg" usemap="#{f}" height=40></td>')
        c+=1
    print("</tr>")


if False:
    # full per single
    for f in fn_order:
        fnsplit = f.split("-")
        prot_name = fnsplit[1]
        if len(aux_data[f])>21:
            org_name = aux_data[f][1:20]+"."
        else:
            org_name = aux_data[f][1:-2]
        print(f'<tr><td>{args.directory}</td><td>{prot_name} <br> {org_name}<br> <a href="http://147.251.115.132/results?ID={prot_name}_7.2_4">AlphaCharges</a><br> <a href="https://alphafold.ebi.ac.uk/entry/{prot_name}">AlphaFold</a></td>')
        c = 0
        fsplit = f.split("-")
        name = fsplit[1]
        width = 220
        height = 75
#        text = ["1", "2", "3", "4a", "4b", "5", "6a", "6b", "7", "8", "9a", "9b", "10"]
        dwg = svgwrite.Drawing(f"{new_dir}/charges-{prot_name}-full-{prot_name}.svg", size=(width, height))
        dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
        graph_full(data, aux_data, sses_order, [f], 2, True, data_all_res[f])
        dwg.save()
        print(f'<td bgcolor="#808080"><img src="{args.directory}/charges-{prot_name}-full-{prot_name}.svg" usemap="#{f}" height=120></td>')
        c+=1
        print("</tr>")


if False:
    # avg full
    print("<tr><td></td><td></td><td>1</td><td>2</td><td>3</td><td>4a</td><td>4b</td><td>5</td><td>6a</td><td>6b</td><td>7</td><td>8</td><td>9a</td><td>9b</td><td>10</td></tr>")
    c = 0
    print(F"<td>{args.directory}</td><td>multi</td>", end=" ")
    width = 220
    height = 75
#    text = ["1", "2", "3", "4a", "4b", "5", "6a", "6b", "7", "8", "9a", "9b", "10"]
    dwg = svgwrite.Drawing(f"{new_dir}/charges-{name}-full-avg.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
    first = True
    for f in fn_order:
        fsplit = f.split("-")
        name = fsplit[1]
        graph_full(data, aux_data, sses_order, [f], len(fn_order), first, data_all_res[f])
#   def graph_full(data, aux_data, sses_order, fn_order, avg, first, data_all_res):

        first = False
        dwg.save()
#    print(f'<td bgcolor="#808080"><img src="{args.directory}/charges-avg-{text[c]}.svg" height=120></td>')
    c+=1
    print("</tr>")

