#!/usr/bin/env python3
# tool which generate for given alphafold family
# its diagrams


import sys
import argparse
import time
import os
import subprocess
import shutil
import json

import utils_common
import svgwrite


def dolu(x,y,res_list, sses, opacity, text="", mmax=0, avg=False):
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
    if avg == False:
        dwg.add(dwg.text(res_min, insert=(x+3, y[0]+5), color="white", font_size=6, style="text-anchor:center"))
        dwg.add(dwg.text(f"({res_max-res_min})", insert=(x+3, y[0]+12), color="white", font_size=5, style="text-anchor:center"))
        dwg.add(dwg.text(res_max, insert=(x+3, y[1]), color="white", font_size=6, style="text-anchor:center"))


def nahoru(x,y,res_list, sses, opacity, text="", mmax=0, avg=False):
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
    if avg == False:
        dwg.add(dwg.text(res_min, insert=(x+3, y[1]), color="white", font_size=6, style="text-anchor:center"))
        dwg.add(dwg.text(f"({res_max-res_min})", insert=(x+3, y[0]+12), color="white", font_size=5, style="text-anchor:center"))
        dwg.add(dwg.text(res_max, insert=(x+3, y[0]+5), color="white", font_size=6, style="text-anchor:center"))


def graph_full(data, aux_data, sses_order, fn_order, avg):
#    print("sses_order", sses_order)
#    print("fn_order", fn_order)
    if avg:
        opacity = 1/len(fn_order)*2
    else:
        opacity = 1
    c_sses = 0
    c_fn = 0
    if avg:
        dwg.add(dwg.text("multiple transparent - ", insert=(10, 70-20), color="white", font_size=9, style="text-anchor:center"))
    else:
        for name in fn_order:
            c_fn += 1
            if len(aux_data[name])>21:
                org_name = aux_data[name][1:20]+"."
            else:
                org_name = aux_data[name][1:-2]
#            dwg.add(dwg.text(org_name, insert=(10, c_fn*70-30), color="white", font_size=9, style="text-anchor:center"))
#            dwg.add(dwg.text(name[7:-5], insert=(10, c_fn*70-20), color="white", font_size=9, style="text-anchor:center"))
##    fnsplit=name.split("-") 
##    fnn = fnsplit[1]
##    print(f'<tr><td>{args.directory}</td><td><a href="https://alphacharges.ncbr.muni.cz/results?ID={fnn}_7.2_4">{fnn}</a><br>{org_name}<br></td><td><img src="{args.directory}/charges-{fnn}_full_avg.svg" width=400></td></tr>')
##    cc = 90
    cc = 0
    for sses in sses_order:
         c_fn = 0
         c_sses += 1
         for name in fn_order:
              if avg:
                  mmin=10 +c_fn*0
                  mmax=60 +c_fn*0
                  dwg.add(dwg.line(start=(100-cc,mmax), end=(100, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(120-cc,mmax), end=(120, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(160-cc,mmax), end=(160, mmax-20), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(160-cc,mmin), end=(160, mmin+20), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(140-cc,mmax), end=(140, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(180-cc,mmax), end=(180, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(200-cc,mmax), end=(200, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(220-cc,mmax), end=(220, mmax-20), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(220-cc,mmin), end=(220, mmin+20), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(240-cc,mmax), end=(240, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(260-cc,mmax), end=(260, mmin), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(280-cc,mmax), end=(280, mmax-22), stroke_width=4, stroke="white"))
                  dwg.add(dwg.line(start=(280-cc,mmin), end=(280, mmin+22), stroke_width=4, stroke="white"))
    for sses in sses_order:
         c_fn = 0
         c_sses += 1
         for name in fn_order:
              if avg:
                  mmin=10 +c_fn*0
                  mmax=60 +c_fn*0
              else:
                  mmin=10 +c_fn*70
                  mmax=60 +c_fn*70
              c_fn += 1
              sses_list = ["H2", "H3", "H5", "H6", "H8", "H9", "H19", "H20", "H21", "H22", "H23", "H25", "H26"]
#              sses_list = ["H0", "H1", "H3", "H4", "H6", "H7", "H11", "H12", "H13", "H14", "H15", "H17", "H18"] # pin5_more
#              sses_list = ["H0", "H1", "H3", "H4", "H6", "H7", "H20", "H21", "H22", "H23", "H24", "H26", "H27"] # pin1_less
#              sses_list = ["H0", "H1", "H3", "H4", "H6", "H7", "H9", "H10", "H11", "H12", "H13", "H15", "H16"] # pin8_less
              if name in data[sses]:
                   if sses == sses_list[0]:
                       xx = 100 + 0-cc
                       dolu(xx, (mmin,mmax), data[sses][name], sses, opacity, "1", mmin, avg)
#                       dwg.add(dwg.text("1", insert=(10, 70-20), color="white", font_size=9, style="text-anchor:center"))
                   elif sses == sses_list[1]:
                       xx = 100 + 20-cc
                       nahoru(xx, (mmin,mmax), data[sses][name], sses, opacity, "2", mmin, avg)
                   elif sses == sses_list[2]:
                       xx = 100 + 40-cc
                       dolu(xx, (mmin,mmax), data[sses][name], sses, opacity, "3", mmin, avg)
                   elif sses== sses_list[4]:
                       xx = 100 + 60-cc
                       nahoru(xx, (mmin,mmin+20), data[sses][name], sses, opacity, "4", mmin, avg)
                   elif sses== sses_list[3]:
                       xx = 100 + 60-cc
                       nahoru(100+60-cc, (mmax-20,mmax), data[sses][name], sses, opacity, "4", mmin, avg)
                   elif sses== sses_list[5]:
                       dolu(100+80-cc, (mmin,mmax), data[sses][name], sses, opacity, "5", mmin, avg)
                   elif sses== sses_list[12]:
                       nahoru(100+100-cc, (mmin,mmax), data[sses][name], sses, opacity,"10", mmin, avg)
                   elif sses== sses_list[11]:
                       dolu(100+120-cc, (mmax-20,mmax), data[sses][name], sses, opacity, "9", mmin, avg)
                   elif sses== sses_list[10]:
                       dolu(100+120-cc, (mmin,mmin+20), data[sses][name], sses, opacity, "9", mmin, avg)
                   elif sses== sses_list[9]:
                       xx = 100 + 140-cc
                       nahoru(xx, (mmin,mmax), data[sses][name], sses, opacity, "8", mmin, avg)
                   elif sses== sses_list[8]:
                       xx = 100 + 160-cc
                       dolu(xx, (mmin,mmax), data[sses][name], sses, opacity, "7", mmin, avg)
                   elif sses== sses_list[7]:
                       xx = 100 + 180-cc
                       nahoru(xx, (mmin,mmin+22), data[sses][name], sses, opacity, "6", mmin, avg)
                   elif sses== sses_list[6]:
                       xx = 100 + 180-cc
                       nahoru(xx, (mmax-22,mmax), data[sses][name], sses, opacity, "6", mmin, avg)
              
              dwg.add(dwg.line(start=(100-cc,mmax), end=(110-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(120-cc,mmax), end=(110-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(140-cc,mmax), end=(150-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(160-cc,mmax), end=(150-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(200-cc,mmax), end=(210-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(220-cc,mmax), end=(210-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(240-cc,mmax), end=(250-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(260-cc,mmax), end=(250-cc, mmax+4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(120-cc,mmin), end=(130-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(140-cc,mmin), end=(130-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(160-cc,mmin), end=(170-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(180-cc,mmin), end=(170-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(220-cc,mmin), end=(230-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(240-cc,mmin), end=(230-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(260-cc,mmin), end=(270-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(280-cc,mmin), end=(270-cc, mmin-4), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(160-cc,mmin+20), end=(160-cc, mmax-20), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(220-cc,mmin+20), end=(220-cc, mmax-20), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(280-cc,mmin+22), end=(280-cc, mmax-22), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(180-cc,mmax), end=(230-cc, mmax+12), stroke_width=1, stroke="black"))
              dwg.add(dwg.line(start=(280-cc,mmax), end=(230-cc, mmax+12), stroke_width=1, stroke="black"))
              


def graph2(data, aux_data, sses_order, fn_order, avg, label, res_data):
    if avg:
        dwg.add(dwg.text("multiple transparent - ", insert=(10, 70-20), color="white", font_size=9, style="text-anchor:center"))
    else:
        for name in fn_order:
            if len(aux_data[name])>21:
                org_name = aux_data[name][1:20]+"."
            else:
                org_name = aux_data[name][1:-2]
    for sses in sses_order:
         for name in fn_order:
              opacity = 1/avg*2
              if name in data[sses]:
                  res_min = None
#                  res_max = None
                  for i in data[sses][name]:
                      if res_min == None or res_min > int(i):
                          res_min = int(i)
                 
#                      if res_max == None or res_max < int(i):
#                          res_max = int(i)
#                      dwg.add(dwg.text(res_min, insert=(5+c_sses*100+ len(data[sses][name])*1+10-cc-11, c_fn*20-3), color="white", font_size=6, style="text-anchor:center"))
#                      dwg.add(dwg.text(f"{text}", insert=(5+c_sses*100+ len(data[sses][name])*1+10-cc, c_fn*20-3), color="white", font_size=8, style="text-anchor:center"))
#                      dwg.add(dwg.text(res_max, insert=(5+c_sses*100+ len(data[sses][name])*1+10-cc+10, c_fn*20-3), color="white", font_size=6, style="text-anchor:center"))
                  if label in ["1","3", "5", "7"]:
                      dwg.add(dwg.text("N", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
                      dwg.add(dwg.text("C", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
                  elif label in ["2", "8", "10"]:
                      dwg.add(dwg.text("C", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
                      dwg.add(dwg.text("N", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
                  elif label in ["4a", "6a"]:
                      dwg.add(dwg.text("C", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
                  elif label in ["4b", "6b"]:
                      dwg.add(dwg.text("N", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
                  elif label in ["9a"]:
                      dwg.add(dwg.text("N", insert=(5, 8), color="white", font_size=8, style="text-anchor:center"))
                  elif label in ["9b"]:
                      dwg.add(dwg.text("C", insert=(2+(int(i)-res_min)*2, 8), color="white", font_size=8, style="text-anchor:center"))
                  for i in data[sses][name]:
                   if res_data[name][int(i)] in ["GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP" ]:
                          #nonpolar
                          color = "pink"
                   elif res_data[name][int(i)] in ["SER", "THR", "TYR", "ASN", "CYS", "GLN"]:
                          #polar
                          color = "red"
                   elif res_data[name][int(i)] in ["ASP", "GLU"]:
                          #negativly charged
                          color = "#0072CF"
                   elif res_data[name][int(i)] in ["LYS", "ARG", "HIS"]:
                       #positively charged
                          color = "green"
                   else:
                          print("neni uveden ", res_data[name][i])#, i , name)
                          color = "black"
                   dwg.add(dwg.line(start=(5+(int(i)-res_min)*2, 12), end=(5+(int(i)-res_min+1)*2, 12), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))

def graph(data, aux_data, sses_order, fn_order, avg, res_data):
#    print("sses_order", sses_order)
#    print("fn_order", fn_order)
    c_sses = 0
    c_fn = 0
#    for name in fn_order:
#        c_fn += 1
#        if not avg:
#            dwg.add(dwg.text(name[7:-5], insert=(10, c_fn*20+5), color="white", font_size=9, style="text-anchor:center"))
#        else:
#            dwg.add(dwg.text(f"multi-avg ({len(fn_order)}pcs)", insert=(10, 1*20+5), color="white", font_size=9, style="text-anchor:center"))
    counter = 1
    cc = 100
    name = fn_order[0]
    fnsplit=name.split("-") 
    fnn = fnsplit[1]
    if avg:
        dwg.add(dwg.text("multiple transparent - ", insert=(10, 70-20), color="white", font_size=9, style="text-anchor:center"))
    else:
        for name in fn_order:
            c_fn += 1
            if len(aux_data[name])>21:
                org_name = aux_data[name][1:20]+"."
            else:
                org_name = aux_data[name][1:-2]

    print(f'<tr><td>{args.directory}</td><td><a href="https://alphacharges.ncbr.muni.cz/results?ID={fnn}_7.2_4">{fnn}</a><br>{org_name}<br></td><td><img src="{args.directory}/charges-{fnn}.svg" width=1600></td></tr>')
##
    for sses in sses_order:
         c_fn = 0
         c_sses += 1
#         counter = 1
         for name in fn_order:
              if not avg:
                  c_fn += 1
                  opacity = 1
              else:
                  c_fn = 1
                  opacity = 1/len(fn_order)*2
              if name in data[sses]:
                  res_min = None
                  for i in data[sses][name]:
                      if res_min == None or res_min > int(i):
                          res_min = int(i)
                  if True:#not avg:
#                      dwg.add(dwg.text(sses, insert=(5+c_sses*100, c_fn*20+5), color="white", font_size=8, style="text-anchor:center"))
#                      dwg.add(dwg.text(res_min, insert=(5+c_sses*100, c_fn*20+5), color="white", font_size=8, style="text-anchor:center"))
                      if counter == 1:
                          text = "1"
                      elif counter == 2:
                          text = "2"
                      elif counter == 3:
                          text = "3"
                      elif counter == 4:
                          text = "4a"
                      elif counter == 5:
                          text = "4b"
                      elif counter == 6:
                          text = "5"
                      elif counter == 7:
                          text = "6a"
                      elif counter == 8:
                          text = "6b"
                      elif counter == 9:
                          text = "7"
                      elif counter == 10:
                          text = "8"
                      elif counter == 11:
                          text = "9a"
                      elif counter == 10:
                          text = "8"
                      elif counter == 11:
                          text = "9a"
                      elif counter == 12:
                          text = "9b"
                      elif counter == 13:
                          text = "10"
                      else:
                          text = str(counter)
                      dwg.add(dwg.text(text, insert=(5+c_sses*100-cc+20, c_fn*20-3), color="white", font_size=8, style="text-anchor:center"))
                  for i in data[sses][name]:
#                      print("rrr", res_data[name], name, i)
                      j = int(i)
                      
                      if res_data[name][j] in ["GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP" ]:
                          #nonpolar
                          color = "pink"
                      elif res_data[name][j] in ["SER", "THR", "TYR", "ASN", "CYS", "GLN"]:
                          #polar
                          color = "red"
                      elif res_data[name][j] in ["ASP", "GLU"]:
                          #negativly charged
                          color = "white"
                      elif res_data[name][j] in ["LYS", "ARG", "HIS"]:
                          #positively charged
                          color = "green"
                      else:
                          print("neni uveden ", res_data[name][j])#, i , name)
                          color = "black"

#                      print("---", res_min, int(i)-res_min, c_sses, c_fn, "---", (c_sses*80+(int(i)-res_min)*2, c_fn*20), (c_sses*80+(int(i)-res_min+1)*2,c_fn*20+2))
                      dwg.add(dwg.line(start=(20+c_sses*100+(int(i)-res_min)*2-cc, c_fn*20), end=(20+c_sses*100+(int(i)-res_min+1)*2-cc,c_fn*20+2), stroke_width=4, stroke=color, style=f"opacity:{opacity}"))
         counter += 1


def diff(data, sses_order, fn_order, aux):
    stats = {}
    for sses in sses_order:
#         print("stats", stats)
         if sses not in stats:
#             print("neni tam ", sses)
             stats[sses] = {}
#         print("stats", stats)
#         print("sses", sses)
         for name in data[sses]:
#              print("name", name)
##              if name in data[sses]:
              print(f"{aux}_{sses}_{name}", end=" ")
              res_min = None
              res_max = None
              for i in data[sses][name]:
                  if res_min == None or res_min > int(i):
                      res_min = int(i)
                  if res_max == None or res_max < int(i):
                      res_max = int(i)
#              print("sses", sses, res_max-res_min)
              for c in data[sses][name]:
#                       print("stats", stats[sses], data[sses])
#                       print("ssss",  data[sses][name])
                       if int(c)-res_min not in stats[sses]:
                           stats[sses][int(c)-res_min] = []
                       print(f"{data[sses][name][c]}", end=" ")
#                       print("c", c, int(c)-res_min, data[sses][name][c])
                       stats[sses][int(c)-res_min].append(data[sses][name][c])
              print()
    
#         print("===============")
#         print(stats, stats[sses])
#    print("==================")
    stats_avg = {}
    stats_d = {}
    for sses in sses_order:
        for c in stats[sses]:
            ssum = 0
            for l in stats[sses][c]:
                ssum += l
            if sses not in stats_avg:
                stats_avg[sses] = {}
                stats_d[sses] = {}
            stats_avg[sses][c] = int(ssum/len(stats[sses][c]))

            ssum = 0
            for l in stats[sses][c]:
#                print("ll ", ssum,  l, stats[sses][c])
                ssum += abs(l-stats_avg[sses][c])

            stats_d[sses][c] = int(ssum/len(stats[sses][c]))
    for sses in stats_avg:
        
#        print("avg", stats_avg[sses])
#        print("d  ", stats_d[sses])
        for c in stats_avg[sses]:
            print(f"{stats_avg[sses][c]}", end= " ")
        print()
    print("---")
    for sses in stats_avg:
        for c in stats_d[sses]:
            print(f"{stats_d[sses][c]}", end= " ")
        print()
    sys.exit()

# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--directory', required=True, action="store",
                    help="name of the directory, which contains domain charges layout")
parser.add_argument('--org_name_dir', required=False, action="store",
                    help="name of the directory, which contains domain charges layout")
#parser.add_argument('--full', action="store_true",
#                    help="short or full output, default short")
#parser.add_argument('--avg', action="store_true",
#                    help="short or full output, default short")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
if args.directory[-1] =="/":
    name = args.directory[:-1]
else:
    name = args.directory

new_dir = f"charges-{name}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)
print("new_dir", new_dir)


data = {}
aux_data = {}
res_data = {}
fn_order = []
for filename in os.listdir(args.directory):
    print("filename", filename)
    res_data[filename] = {}
    with open(f"{args.directory}/{filename}", "r") as fd:
        jdata = json.load(fd)
#        print("data", jdata["h_type"])
        fn_order.append(filename)
        for sse in jdata["h_type"]:
            if jdata["h_type"][sse] == 1 and (sse[0] != "E"):
                if sse not in data:
                    data[sse] = {}
                data[sse][filename] = jdata["charges"][sse]
                if args.org_name_dir:
                     sfilename=filename.split("-")
#                     print("name is")
                     with open(f"{args.org_name_dir}/{sfilename[1]}-F1.cif") as fd:
                         line = fd.readline()
                         while line:
                             if line.startswith("_ma_target_ref_db_details.organism_scientific"):
                                  aux_data[filename] = line[55:]
                                  break
                             line = fd.readline()
                else:
                     aux_data[filename] = "--"

                if args.org_name_dir:
                     sfilename=filename.split("-")
#                     print("name is")
                     with open(f"{args.org_name_dir}/{sfilename[1]}-F1.pdb") as fd:
                         line = fd.readline()
                         while line:
                             if line.startswith("ATOM"):
#                                  print("line", line)
#                                  print(f"line[22:26] >{line[22:26]}<")
#                                  print(f"line[17:20] >{line[17:20]}<")
                                  resnum = int(line[22:26])
                                  restype = line[17:20]
#                                  print("pom", res_data[filename], resnum, type(resnum))
#                                  print("pred")
                                  res_data[filename][resnum] = restype
                                 # break
                             line = fd.readline()
                else:
                     aux_data[filename] = "--"



fn_order.sort()


sses_order = list(data.keys())

print("sses order", sses_order)

#for sse in data:
#    print("---", sse, data[sse])
# we want to download/read relevant cif files

##diff(data, sses_order, fn_order, args.directory[-1])

if False:#True:
    print("1---------------------------")
    width = (len(sses_order)+1)*20+70
    height =2*70-50
    dwg = svgwrite.Drawing(f"{new_dir}/charges-{name}_full_avg.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
    graph_full(data, aux_data, sses_order, fn_order, True)
    dwg.save()
    print("2---------------------------")

    height =(len(fn_order)+1)*70-50
    dwg = svgwrite.Drawing(f"{new_dir}/charges-{name}_full.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
    graph_full(data, aux_data, sses_order, fn_order, False)
    dwg.save()

if False: #True:
    print("3---------------------------")
    width = (len(sses_order)+1)*100
    height =2*70-50-50
    print("name", f"{new_dir}/charges-{name}_avg.svg")
    dwg = svgwrite.Drawing(f"{new_dir}/charges-{name}_avg.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
    graph(data, aux_data, sses_order, fn_order, True)
    dwg.save()
    sys.exit()
    print("4---------------------------")
    height =(len(fn_order)+1)*20
    dwg = svgwrite.Drawing(f"{new_dir}/charges-{name}.svg", size=(width, height))
    dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
    graph(data, aux_data, sses_order, fn_order, False)

print("=======================")
if False:
    # single images
    for f in fn_order:
        fsplit = f.split("-")
        name = fsplit[1]
#        print(f"<tr><> </tr>"
        width = (len(sses_order)+1)*100-100
        height =2*70-50-50
        dwg = svgwrite.Drawing(f"{new_dir}/charges-{name}.svg", size=(width, height))
        dwg.add(dwg.rect(insert=(1, 1), size=(width, height), fill="grey"))
        graph(data, aux_data, sses_order, [f], False, res_data)
        dwg.save()
        

if True:
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
            graph2(data, aux_data, [sse], [f], 2, text[c], res_data)
            dwg.save()
            print(f'<td bgcolor="#808080"><img src="{args.directory}/charges-{prot_name}-{text[c]}.svg" height=40></td>')
            c+=1
        print("</tr>")


