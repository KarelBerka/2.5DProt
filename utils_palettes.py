#!/usr/bin/env python

def revert_cs(cs):
    r_cs = []
    for c in reversed(cs):
       new_c = f"#{c[5]}{c[6]}{c[1]}{c[2]}{c[3]}{c[4]}"
       new_c = f"#{c[3]}{c[4]}{c[5]}{c[6]}{c[1]}{c[2]}"
       r_cs.append(new_c)
    return r_cs


def domain_order(name):
    if (len(name) == 1):
        order=ord(name[0])-ord('A')
    elif (len(name) == 2):
        order = 26*(ord(name[0])-ord('A')+1)+(ord(name[1])-ord('A'))
    elif (len(name) == 3):
        order = 26*26*(ord(name[0])-ord('A')+1) + 26*(ord(name[0])-ord('A')+1) + (ord(name[1])-ord('A'))
    else:
        order = 26*26*(ord(name[0])-ord('A')+1) + 26*(ord(name[0])-ord('A')+1) + (ord(name[1])-ord('A'))
    return order


def set_domain_color(name):
   do_modulo = domain_order(name) % len(protein_domain_colors)
   return protein_domain_colors[do_modulo]

color_scheme = {}

color_scheme[10] = [
    "#006400",
    "#00008b",
    "#b03060",
    "#ff0000",
    "#ffff00",
    "#00ff00",
    "#00ffff",
    "#ff00ff",
    "#6495ed",
    "#ffdab9" ]

color_scheme[20] = [
    "#ffe4c4",
    "#2e8b57",
    "#800000",
    "#808000",
    "#00008b",
    "#ff0000",
    "#ff8c00",
    "#ffd700",
    "#ba55d3",
    "#00ff7f",
    "#0000ff",
    "#f08080",
    "#adff2f",
    "#ff00ff",
    "#1e90ff",
    "#dda0dd",
    "#87ceeb",
    "#ff1493",
    "#7fffd4",
    "#696969"]

color_scheme[30] = [
    "#ffb6c1",
    "#556b2f",
    "#8b4513",
    "#483d8b",
    "#008000",
    "#008b8b",
    "#4682b4",
    "#9acd32",
    "#00008b",
    "#daa520",
    "#8fbc8f",
    "#800080",
    "#b03060",
    "#ff0000",
    "#ffff00",
    "#7cfc00",
    "#8a2be2",
    "#00ff7f",
    "#dc143c",
    "#00ffff",
    "#0000ff", 
    "#da70d6",
    "#ff7f50",
    "#ff00ff",
    "#1e90ff",
    "#90ee90",
    "#add8e6",
    "#ff1493",
    "#ffdead",
    "#696969"]

color_scheme[40] = [
    "#7f0000",
    "#ffb6c1",
    "#556b2f",
    "#a0522d",
    "#808000",
    "#483d8b",
    "#008000",
    "#008080", 
    "#4682b4",
    "#000080",
    "#9acd32",
    "#32cd32",
    "#daa520",
    "#8fbc8f",
    "#8b008b",
    "#d2b48c",
    "#ff4500",
    "#00ced1",
    "#ff8c00",
    "#00ff00",
    "#00fa9a",
    "#8a2be2",
    "#dc143c",
    "#00bfff",
    "#f4a460",
    "#9370db",
    "#0000ff",
    "#adff2f",
    "#ff00ff",
    "#1e90ff",
    "#db7093",
    "#f0e68c",
    "#fa8072",
    "#ffff54",
    "#dda0dd",
    "#add8e6",
    "#ff1493",
    "#ee82ee",
    "#7fffd4",
    "#808080" ]

protein_domain_colors = [
"#1b9e77",
"#d95f02",
"#7570b3",
"#e7298a",
"#66a61e",
"#e6ab02",
"#a6761d",
"#666666",
"#e41a1c",
"#377eb8",
"#4daf4a",
"#984ea3",
"#ff7f00",
"#ffff33",
"#a65628",
"#f781bf",
"#999999",
"#66c2a5",
"#fc8d62",
"#8da0cb",
"#e78ac3",
"#a6d854",
"#ffd92f",
"#e5c494",
"#b3b3b3",]
