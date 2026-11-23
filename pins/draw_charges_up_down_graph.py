#!/usr/bin/env python3
# tool which generates diagram like:
# http://147.251.21.23/pub/2023_12_05_piny/up_and_down/diagram-up2.png
import argparse
import time
import os
import svgwrite

const = 2

# read arguments
parser = argparse.ArgumentParser(
    description='generated charges special diagram - simple one')
parser.add_argument('--input', required=True, action="store",
                    help="name of the file containing recordst like: 'pin1	302	397	A0A0H3V2K9'")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
name = args.input

new_dir = f"draw-{name}-{timestamp}"
os.makedirs(new_dir, exist_ok=True)
print("new_dir", new_dir)
print("const", const)

with open(args.input, "r") as fd:
    line = fd.readline()
    width = 1200
    height = 1200
    dwg = svgwrite.Drawing(f"{new_dir}/diagram-up_down-{name}.svg", size=(width*const, height*const))
    dwg.add(dwg.rect(insert=(1, 1), size=(width*const, height*const), fill="white"))
    dwg.add(dwg.line(start=((width/2+500)*const, (height/2-500)*const), end=((width/2-500)*const, (height/2+500)*const), stroke="black", stroke_width=2))
    dwg.add(dwg.line(start=((width/2-500)*const, (height/2)*const), end=((width/2+500)*const, (height/2)*const), stroke="black", stroke_width=2))
    dwg.add(dwg.text(f"{0}", insert=((width/2+3)*const, (height/2-3)*const), color="white", font_size=12, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-500)*const, (height/2-200/2)*const), end=((width/2+500)*const, (height/2-200/2)*const), stroke="grey"))
    dwg.add(dwg.text(f"{2}", insert=((width/2+3)*const, (height/2-200/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-500)*const, (height/2-400/2)*const), end=((width/2+500)*const, (height/2-400/2)*const), stroke="grey"))
    dwg.add(dwg.text(f"{4}", insert=((width/2+3)*const, (height/2-400/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-500)*const, (height/2+200/2)*const), end=((width/2+500)*const, (height/2+200/2)*const), stroke="grey"))
    dwg.add(dwg.text(f"{-2}", insert=((width/2+3)*const, (height/2+200/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2-500)*const, (height/2+400/2)*const), end=((width/2+500)*const, (height/2+400/2)*const), stroke="grey"))
    dwg.add(dwg.text(f"{-4}", insert=((width/2+3)*const, (height/2+400/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2)*const, (height/2+500)*const), stroke="black", stroke_width=2))
    dwg.add(dwg.text("UP - part", insert=((width/2+3)*const, (height/2-300/2-3)*const), color="white", font_size=16, style="text-anchor:center"))

    dwg.add(dwg.line(start=((width/2-200/2)*const, (height/2-500)*const), end=((width/2-200/2)*const, (height/2+500)*const), stroke="grey"))
    dwg.add(dwg.text(f"{-2}", insert=((width/2-200/2+3)*const, (height/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2+200/2)*const, (height/2-500)*const), end=((width/2+200/2)*const, (height/2+500)*const), stroke="grey"))
    dwg.add(dwg.text(f"{2}", insert=((width/2+200/2+3)*const, (height/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2+400/2)*const, (height/2-500)*const), end=((width/2+400/2)*const, (height/2+500)*const), stroke="grey"))
    dwg.add(dwg.text(f"{4}", insert=((width/2+400/2+3)*const, (height/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2+600/2)*const, (height/2-500)*const), end=((width/2+600/2)*const, (height/2+500)*const), stroke="grey"))
    dwg.add(dwg.text(f"{6}", insert=((width/2+600/2+3)*const, (height/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.line(start=((width/2+800/2)*const, (height/2-500)*const), end=((width/2+800/2)*const, (height/2+500)*const), stroke="grey"))
    dwg.add(dwg.text(f"{8}", insert=((width/2+800/2+3)*const, (height/2-3)*const), color="white", font_size=20, style="text-anchor:center"))
    dwg.add(dwg.text("DOWN - part", insert=((width/2+900/2+3)*const, (height/2-3)*const), color="white", font_size=16, style="text-anchor:center"))
    pcRotate = 'rotate(45,%s, %s)' % ((width/2+600/2+3)*const, (height/2+200/2-3)*const)
    pcRotate = 'rotate(45,%s, %s)' % ((width/2+400/2+3)*const, (height/2+400/2-3)*const)

    dwg.add(dwg.line(start=((width/2-650)*const, (height/2-500)*const), end=((width/2+350)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-600)*const, (height/2-500)*const), end=((width/2+400)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-550)*const, (height/2-500)*const), end=((width/2+450)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-500)*const, (height/2-500)*const), end=((width/2+500)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-450)*const, (height/2-500)*const), end=((width/2+550)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-400)*const, (height/2-500)*const), end=((width/2+600)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-350)*const, (height/2-500)*const), end=((width/2+650)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-300)*const, (height/2-500)*const), end=((width/2+700)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-250)*const, (height/2-500)*const), end=((width/2+750)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-200)*const, (height/2-500)*const), end=((width/2+800)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-150)*const, (height/2-500)*const), end=((width/2+850)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-100)*const, (height/2-500)*const), end=((width/2+900)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2-50)*const, (height/2-500)*const), end=((width/2+950)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2+1000)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2+50)*const, (height/2-500)*const), end=((width/2+1050)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2+100)*const, (height/2-500)*const), end=((width/2+1100)*const, (height/2+500)*const), stroke="#DDDDDD", stroke_width=1))
    dwg.add(dwg.line(start=((width/2)*const, (height/2-500)*const), end=((width/2)*const, (height/2+500)*const), stroke="black", stroke_width=2))

    while (line):
        sline = line.split("\t")
        line = fd.readline()
        if sline[0] == "pin1":
            color = "red"
        elif sline[0] == "pin2":
            color = "green"
        elif sline[0] == "pin3":
            color = "blue"
        elif sline[0] == "pin4":
            color = "pink"
        elif sline[0] == "pin5":
            color = "magenta"
        elif sline[0] == "pin6":
            color = "cyan"
        elif sline[0] == "pin7":
            color = "brown"
        elif sline[0] == "pin8":
            color = "orange"
        dwg.add(dwg.line(start=((width/2+int(sline[2])/2-3)*const, (height/2-int(sline[1])/2-3)*const), end=((width/2+int(sline[2])/2+3)*const, (height/2-int(sline[1])/2+5)*const), stroke_width=4, stroke=color))
        dwg.add(dwg.line(start=((width/2+int(sline[2])/2+3)*const, (height/2-int(sline[1])/2-3)*const), end=((width/2+int(sline[2])/2-3)*const, (height/2-int(sline[1])/2+5)*const), stroke_width=4, stroke=color))
        print("sline", sline[3][:-1])
        if sline[3][:-1] == "Q673E8":
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2+3)*const+5), color="white", font_size=12, style="text-anchor:center"))
        elif sline[3][:-1] == "A0A1L2D3F5":
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2+3)*const+7), color="white", font_size=12, style="text-anchor:center"))
        elif sline[3][:-1] == "D7KRJ6":
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2+3)*const-5), color="white", font_size=12, style="text-anchor:center"))
        elif sline[3][:-1] == "D7LLR0":
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2+3)*const+5), color="white", font_size=12, style="text-anchor:center"))
        elif sline[3][:-1] == "Q8LKH1":
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2+3)*const-3), color="white", font_size=12, style="text-anchor:center"))
        elif sline[3][:-1] == "Q9SQH6":
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2-3)*const-3), color="white", font_size=12, style="text-anchor:center"))
        else:
            dwg.add(dwg.text(f"{sline[0][3]}-{sline[3]}", insert=((width/2+int(sline[2])/2+5)*const, (height/2-int(sline[1])/2+3)*const), color="white", font_size=12, style="text-anchor:center"))

dwg.save()


