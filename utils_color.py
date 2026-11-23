#!/usr/bin/env python
import json
import utils_palettes
import utils_color_ligands
import random


def set_ligands_colors(ligands_data, input_ligands_file, args_multi=False, dominik_color_scheme=False):

    if input_ligands_file:
        # if we have an input file with ligand colors use it
        with open(input_ligands_file, "r") as fd:
            color_data = json.load(fd)
        for i in ligands_data:
            for lig in color_data:
                # set color to each ligand type
                if lig["edsc"] == i["edsc"] and lig["e_type"] == i["e_type"] and lig["ecmp"] == i["ecmp"]:
                    i["color"] = lig["color"]
                    i["ocolor"] = lig["ocolor"]
                    break

    else:
        # we do not have input file with ligand colors
        # there are three posibilities:
        if dominik_color_scheme:
            # use Dominik's color scheme use args.dominik
            for i in ligands_data:
                # use args.dominik "ch"
                if dominik_color_scheme == "ch":
                    i["color"] = utils_color_ligands.set_ligand_chemical_color(i["ecmp"])
                    i["ocolor"] = i["color"]

                # use args.dominik "f"
                elif dominik_color_scheme == "f":
                    i["color"] = utils_color_ligands.set_ligand_functional_color(i["ecmp"])
                    i["ocolor"] = i["color"]

        else:
            # there is no input file and we do not want to use dominik color scheme -
            # use palette from utils_palletes
            # choose the best palette, based on ligands number:
            if len(ligands_data) <= 10:
                palette = utils_palettes.color_scheme[10]
            elif len(ligands_data) <= 20:
                palette = utils_palettes.color_scheme[20]
            elif len(ligands_data) <= 30:
                palette = utils_palettes.color_scheme[30]
            else:
                palette = utils_palettes.color_scheme[40]

            ligand_type_counter = 0
            for i in ligands_data:
                # set color to each ligand type
                if ligand_type_counter < 40:
                    # first 40 ligand types uses chosen palette
                    i["color"] = palette[ligand_type_counter]
                    i["ocolor"] = palette[ligand_type_counter]
                else:
                    if args_multi:
                        # 40th and higher ligand types are shown by the same color
                        # - palette[39] (silver one) and labelled as others
                        # if we use diagram with legend (multiple diagram)
                        i["color"] = palette[39]
                        i["ocolor"] = palette[39]
                    else:
                        # otherwise - go on and use random colors
                        red = random.randint(0, 256)
                        green = random.randint(0, 256)
                        blue = random.randint(0, 256)
                        i["color"] = '#%02x%02x%02x' % (red, green, blue)
                        i["ocolor"] = "black"
                ligand_type_counter += 1


def set_channels_colors(channels_data, args_multi=False):

    color_counter = 0
    for i in channels_data:
        palette = utils_palettes.color_scheme[20]
        if i["text_names"] == "Tunnel":
            i["color"] = palette[10]
        elif i["text_names"] == "Wather channel" or i["text_names"] == "Wather channel ":
            i["color"] = palette[1]
        elif i["text_names"] == "Solvent channel" or i["text_names"] == "Solvent channel ":
            i["color"] = palette[2]
        elif i["text_names"] == "Channel 2b" or i["text_names"] == "Channnel 2b":
            i["color"] = palette[3]
        elif i["text_names"] == "Channel 2e":
            i["color"] = palette[4]
        elif i["text_names"] == "Channel 2c":
            i["color"] = palette[5]
        elif i["text_names"] == "Channel 2f":
            i["color"] = palette[6]
        elif i["text_names"] == "Channel 2a" or i["text_names"] == "Channnel 2a":
            i["color"] = palette[7]
        elif i["text_names"] == "Channel 3":
            i["color"] = palette[8]
        elif i["text_names"] == "Channel 2ac":
            i["color"] = palette[9]
        elif i["text_names"] == "Channel 5":
            i["color"] = palette[11]
        elif i["text_names"] == "Channel 2d" or i["text_names"] == "Chnnel 2d":
            i["color"] = palette[12]
        elif i["text_names"] == "Channel 4":
            i["color"] = palette[13]
        elif i["text_names"] == "Channel 2bs":
            i["color"] = palette[14]
        else:
            if (args_multi):
                # multiple diagrams use only one color for tunnels (with transparency)
                i["color"] = "#665646"
            else:
                # single domain color (not too dark to shade ligands)
                match color_counter % 7:
                    case 0:
                        i["color"] = "#88bb88"
                    case 1:
                        i["color"] = "#bb8888"
                    case 2:
                        i["color"] = "#8888bb"
                    case 3:
                        i["color"] = "#b070b0"
                    case 4:
                        i["color"] = "#70b0b0"
                    case 5:
                        i["color"] = "#b0b070"
                    case 6:
                        i["color"] = "#999999"
                color_counter += 1
