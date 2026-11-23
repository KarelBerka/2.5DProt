import math
import utils
import utils_vector
from utils_vector import add
from utils_vector import rm
import utils_sses_types
import utils_image
import utils_color_ligands
import utils_color
import utils_membrane
import loader
try:
    import pymol
except ImportError:
    pymol = None
import superimpose2d
import colorsys
import json
import os
import os.path
import logging
import sys
import numpy as np

MIN_CANVAS_DIMENSIONS = (30, 10)

def fix_conflict(layout, angles, size, h1, h2, beta_con):
    """
    This function found sses which cross each other and fix their positions
    the movement depends on sses sizes and types - if it is strand we have
    to move its whole sheet together
    This function modifies layout variable
    """

    # At first we choose which sse will be moved - we compute the
    # size (size1, size2) of the relevant helix / sheet and find
    # whole structures to which sses belonged (cl1, cl2)
    size1 = 0
    size2 = 0
    for cluster in beta_con:
        cl_sses = cluster.split(":")
        for sse in cl_sses:
            if sse == h1:
                cl1 = cl_sses
                for s in cl_sses:
                    size1 = size1 + size[s]
            if sse == h2:
                cl2 = cl_sses
                for s in cl_sses:
                    size2 = size2 + size[s]

    if (cl1 == cl2):
        # conflict can be fixed only if it is in the different structures
        return

    # we have to choose which structure sse have to be moved
    if utils_sses_types.is_sheet(h1) and utils_sses_types.is_helic(h2):
        # the first is sheet and the second is helic
        sse_to_move, cluster_to_move = h2, cl2
        sse_to_leave = h1
    elif utils_sses_types.is_helic(h1) and utils_sses_types.is_sheet(h2):
        # the first is helic and the second is sheet
        sse_to_move, cluster_to_move = h1, cl1
        sse_to_leave = h2
    else:
        # both have the same type (both are sheets /helices)
        # we move the smaller one
        if (size1 > size2):
            sse_to_move, cluster_to_move = h2, cl2
            sse_to_leave = h1
        else:
            sse_to_move, cluster_to_move = h1, cl1
            sse_to_leave = h2

    # now we find the direction and length of optimal movement
    dx = utils_vector.angle_target([0, 0], (angles[sse_to_leave] + 90) % 360, 1)

    vl = utils_vector.vector_length(dx) * 2
    l1 = [layout[sse_to_move][0] + dx[0] / vl,
          layout[sse_to_move][1] + dx[1] / vl]
    l2 = [layout[sse_to_move][0] - dx[0] / vl,
          layout[sse_to_move][1] - dx[1] / vl]
    if (math.dist(layout[sse_to_leave], l1) <
       math.dist(layout[sse_to_leave], l2)):
        dx = (-dx[0], -dx[1])

    # modify layout variable
    for i in cluster_to_move:
        if i in layout:
            layout[i] = [
                layout[i][0] + dx[0] / vl, layout[i][1] + dx[1] / vl]


def count_start_end(layout, angles, sizes, chosen_sses):
    """ count sses start and end points  """
    s = {}
    e = {}
    d = {}
    for sse in chosen_sses:
        if (utils_sses_types.is_helic(sse)):
            # compute endsize based on the type of sse (sheet symbol is longer)
            endsize = 0.5
        else:
            endsize = 1
        s[sse] = utils_vector.angle_target(layout[sse], angles[sse] + 180,
                                           sizes[sse] * utils_image.Image.CIRCLE_MULTIPLIER + 0.5)
        e[sse] = utils_vector.angle_target(layout[sse], angles[sse],
                                           sizes[sse] * utils_image.Image.CIRCLE_MULTIPLIER + endsize)
        d[sse] = utils_vector.angle_target([0, 0], angles[sse] + 90, 0.8)
    if len(chosen_sses) == 1:
        return s[sse], e[sse], d[sse]
    else:
        return s, e, d


def intersection(se1, se2, h1, h2):
    """ return whether two line segments intersect """
    p = utils_vector.vectors_intersection(se1[0], se1[1], se2[0], se2[1])
    if (p[0] is None) and (p[1] is None):
        return False
    eps = 1e-5
    in_x1 = min(se1[0][0], se1[1][0]) - eps <= p[0] <= max(se1[0][0], se1[1][0]) + eps
    in_y1 = min(se1[0][1], se1[1][1]) - eps <= p[1] <= max(se1[0][1], se1[1][1]) + eps
    in_x2 = min(se2[0][0], se2[1][0]) - eps <= p[0] <= max(se2[0][0], se2[1][0]) + eps
    in_y2 = min(se2[0][1], se2[1][1]) - eps <= p[1] <= max(se2[0][1], se2[1][1]) + eps
    return in_x1 and in_y1 and in_x2 and in_y2


def point_to_point_min_dist(p1, d1, p2, d2):
    dmin = math.dist(add(p1, d1), add(p2, d2))
    dd = math.dist(rm(p1, d1), add(p2, d2))
    if dd < dmin:
        dmin = dd
    dd = math.dist(add(p1, d1), rm(p2, d2))
    if dd < dmin:
        dmin = dd
    dd = math.dist(rm(p1, d1), rm(p2, d2))
    if dd < dmin:
        dmin = dd
    return dmin


def resolve_sses_conflicts(layout, angles, sizes, relevant_sses, beta_con):
    """
    If two SSEs are overlapping, move one of them aside
    """
    start, end, d = count_start_end(layout, angles, sizes, relevant_sses)
    conflict_fixed = True
    countc = 0                                   # solved conflicts counter
    if (len(layout) > 500):
        return layout
    # TODO: use countc maximum as an input parameter for this tool
    # fix only first 10 conflicts to prevent loop
    try:
      while (conflict_fixed == True) and (countc < 10):
        conflict_fixed = False
        countc = countc + 1                      # conflict counter increased

        for h1 in relevant_sses:
            for h2 in relevant_sses:
                if h1 == h2:
                    continue

                # for all different pairs of sses h1 and h2, test wheher their position is not in conflict
                fixed = 0
                countn = 0                       # counter of tries to fix conflict positions h1 x h2
                conflict = False
                while (fixed == 0) and (countn < 10):
                    countn = countn + 1          # one more try to fix

                    if abs((angles[h1] - angles[h2]) % 180) < 1e-4:
                        # both SSEs are parallel or antiparallel
                        dif = [end[h1][0] - start[h1][0], end[h1][1] - start[h1][1]]
                        if ((utils.point_to_line_dist2(start[h1], dif, start[h2]) < 1) and
                           (utils.point_to_line_min(start[h1], end[h1], start[h2]) == 0)) or \
                           (math.dist(start[h1], start[h2]) < 1) or (math.dist(start[h1], end[h2]) < 1) or \
                           (math.dist(end[h1], start[h2]) < 1) or (math.dist(end[h1], end[h2]) < 1):
                            # if there is conflict
                            conflict = True
                    else:
                        # both SSEs are not parallel parallel try to find their intersection
                        i1 = intersection([add(start[h1], d[h1]), add(end[h1], d[h1])],
                                          [add(start[h2], d[h2]), add(end[h2], d[h2])], h1, h2)
                        i2 = intersection([rm(start[h1], d[h1]), rm(end[h1], d[h1])],
                                          [add(start[h2], d[h2]), add(end[h2], d[h2])], h1, h2)
                        i3 = intersection([add(start[h1], d[h1]), add(end[h1], d[h1])],
                                          [rm(start[h2], d[h2]), rm(end[h2], d[h2])], h1, h2)
                        i4 = intersection([rm(start[h1], d[h1]), rm(end[h1], d[h1])],
                                          [rm(start[h2], d[h2]), rm(end[h2], d[h2])], h1, h2)
                        if (i1 or i2 or i3 or i4):
                            # intersection means there is a conflict
                            conflict = True
                        else:
                            s1s2 = point_to_point_min_dist(start[h1], d[h1], start[h2], d[h2])
                            s1e2 = point_to_point_min_dist(start[h1], d[h1], end[h2], d[h2])
                            e1s2 = point_to_point_min_dist(end[h1], d[h1], start[h2], d[h2])
                            e1e2 = point_to_point_min_dist(end[h1], d[h1], end[h2], d[h2])
                            if ((s1s2 < 0.7) or (s1e2 < 0.7) or (e1s2 < 0.7) or (e1e2 < 0.7)):
                                # SSEs are to near, classify it as conflict
                                conflict = True

                    if conflict == True:
                        # conflict was found - fix it and update layout
                        fix_conflict(layout, angles, sizes, h1, h2, beta_con)
                        [start[h1], end[h1], d[h1]] = count_start_end(layout, angles, sizes, [h1])
                        [start[h2], end[h2], d[h2]] = count_start_end(layout, angles, sizes, [h2])
                        conflict = False
                        conflict_fixed = True                    # we have fixed conflict - we have to test SSE pairs in next round
                    else:
                        fixed = 1
      return layout
    except:
      return layout


def write_ligands_legend(img, ligands_data, ligands, dimensions, args):
    """ draw ligands legend (only if ligands are shown and its description is needed in add_description tag)"""

    # ligand legend should be indented to the center of the diagram
    right_border = img.get_upper_right()[0]
    left_border = img.get_lower_left()[0]
    image_width = right_border - left_border

    # the number of ligands legend columns
    if len(ligands_data) % 5 > 0:
        columns = (len(ligands_data) // 5) + 1
    else:
        columns = len(ligands_data) // 5

    ligand_legend_width = min(8*11, columns*11)                         # 8 is max number of columns
    ligand_legend_indentation = (image_width - ligand_legend_width) / 2

    # the position of ligand in legend
    n = 2   # row
    m = 0   # column

    # for all versions except of dominik type of output - dominik ligand colors
    if len(ligands_data) > 5:
        const = 10
    else:
        const = 2 * len(ligands_data)

    # use dominik coloring
    if args.dominik == "f" or args.dominik == "ch":
        ligands_types = {}
        for l in ligands_data:
           if args.dominik == "ch":
                # for dominik ligand colors output chemical multi diagram with ligands
                ltype = utils_color_ligands.set_ligand_chemical_type(l["ecmp"])
                color = utils_color_ligands.set_ligand_chemical_color(l["ecmp"])
           else:
                # for dominik ligand colors output  functional  multi diagram with ligands
                ltype = utils_color_ligands.set_ligand_functional_type(l["ecmp"])
                color = utils_color_ligands.set_ligand_functional_color(l["ecmp"])

           if ltype not in ligands_types:
               ligands_types[ltype] = {"count":0, "color": color, "name": ltype}
               ligands_types[ltype]["count"] += l["count"]

        if len(ligands_types) > 5:  # do not overflow space - labels are longer
            const = 10
        else:
            const = 2 * len(ligands_types)

        for l in ligands_types:
            desc_from = utils_vector.add(img.get_upper_left(), (1.2 + m * 11 * 1.5 + ligand_legend_indentation, 2-img.dimension[1] -n + const))
            img.ligand(name=ligands_types[l]["count"], coord=desc_from, size=10,
                e_type=ligands_types[l]["name"], edsc=ligands_types[l]["name"], ecmp=ligands_types[l]["name"],
                add_descriptions=None, input_color=ligands_types[l]["color"], ocolor=ligands_types[l]["color"])

            desc_from = utils_vector.add(img.get_upper_left(), (1.6 + m * 11 *1.5 + ligand_legend_indentation, 2-img.dimension[1] -n + const))
            img.text('...%s (%s pcs)' % (ligands_types[l]["name"], ligands_types[l]["count"]), desc_from,
                 style='fill: black; font-family: Verdana; text-anchor: left;')
            if n == 10:    # new row
               n = 2
               m += 1
            else:          # new column
               n += 2
    else:
        # use random coloring
        for l in ligands_data:
            desc_from = utils_vector.add(img.get_upper_left(), (1.2 + m * 11 + ligand_legend_indentation, 2-img.dimension[1] -n + const))
            img.ligand(name=l["count"], coord=desc_from, size=l["size"],
                e_type=l["e_type"], edsc=l["edsc"], ecmp=l["ecmp"], 
                add_descriptions=None, input_color=l["color"], ocolor =l["ocolor"])

            r = img.ligand_diagram_size(l["size"])
            desc_from = utils_vector.add(img.get_upper_left(), (1.3 + m * 11 + ligand_legend_indentation+r, 2-img.dimension[1] -n + const))
            name = img.lig_text(l["edsc"], l["ecmp"])
            img.text('...%s (%s pcs)' % (name, l["count"]), desc_from,
                 style='fill: black; font-family: Verdana; text-anchor: left;')
            if n == 10:   # new row
               n = 2
               m += 1
            else:         # new column
               n += 2
            if m == 7 and n == 10 and len(ligands_data) > 40:
               # there are more then 39 ligand types, all legend other
               # members of ligand data should be put to the 40 place
               # to label "OTHERS"
               other_ligands_sum = 0
               for l in ligands_data[40:]:
                   other_ligands_sum += l["count"]
               desc_from = utils_vector.add(img.get_upper_left(), (1.2 + m * 11 + ligand_legend_indentation, 2-img.dimension[1] - n + const))
               img.ligand(name=l["count"], coord=desc_from, size=l["size"],
                   e_type=l["e_type"], edsc=l["edsc"], ecmp=l["ecmp"],
                   add_descriptions=None, input_color=l["color"], ocolor =l["ocolor"])

               desc_from = utils_vector.add(img.get_upper_left(), (1.6 + m * 11 + ligand_legend_indentation, 2-img.dimension[1] - n + const))
               img.text('...%s (%s pcs)' % ("OTHERS", other_ligands_sum), desc_from,
                 style='fill: black; font-family: Verdana; text-anchor: left;')
               break
               # the loop have to finish now, nothing more to show


def show_full_layout(img, layout, ind, helices_residues, angles, sizes,
                     pdb_file, chain_id, chain, domain, h_type, colors,
                     add_descriptions, opacity, img_data, spike, ligands,
                     ligands_data, img_data_ligands, args, threeddata, charges, am, membrane=None):
    """ Draw full layout -
        text info if it is necessary
    """
    #TODO: img_data_ligands should be removed
    chain_groups = {}
    chain_groups_links = {}
    chain_groups_sses = {}
    tm_helices = set(membrane.get("tm_helices", [])) if membrane else set()

    # decode color tag
    for sse in colors:
        if colors[sse][0] == 's':
            c = colors[sse]
            if pymol is not None:
                color_tuple = pymol.cmd.get_color_tuple(colors[sse])
            else:
                try:
                    val = float(c[1:]) / 999.0
                except ValueError:
                    val = 0.5
                # PyMOL spectrum 's000'..'s999' goes from blue (240 deg) through green to red (0 deg)
                hue = (1.0 - max(0.0, min(1.0, val))) * (240.0 / 360.0)
                color_tuple = colorsys.hsv_to_rgb(hue, 0.85, 0.90)
            color_tuple_p = (int(color_tuple[0] * 255),
                             int(color_tuple[1] * 255),
                             int(color_tuple[2] * 255))
            colors[sse] = '#%02x%02x%02x' % color_tuple_p

    # Filter only SSEs we are going to display and sort them based on res. numbers
    sses_to_display = []
    for h in layout.keys():
        res = utils_sses_types.sse_should_be_shown(h, h_type, sizes)
        if ( res == 0):
            sses_to_display.append(h)

    # Sort sses_to_display based on sse order in protein
    sses_to_display_unsorted = \
        dict([(i,sses_to_display[i]) for i in range(len(sses_to_display))])

    # at first sort to chains
    aux_sses_list = {}
    for h in sses_to_display_unsorted:
         if sses_to_display_unsorted[h] in chain_id:
             ch = chain_id[sses_to_display_unsorted[h]]
             if ch not in aux_sses_list:
                aux_sses_list[ch] = [sses_to_display_unsorted[h]]
             else:
                 aux_sses_list[ch].append(sses_to_display_unsorted[h])

    sses_to_display = {}
    for ch in aux_sses_list:
        sses_to_display_unsorted[ch] = \
            dict([(i,aux_sses_list[ch][i]) for i in range(len(aux_sses_list[ch]))])
        for h in sses_to_display_unsorted[ch]:
            for j in sses_to_display_unsorted[ch]:
               if helices_residues[sses_to_display_unsorted[ch][h]][0] < helices_residues[sses_to_display_unsorted[ch][j]][0]:
                   pom = sses_to_display_unsorted[ch][j]
                   sses_to_display_unsorted[ch][j] = sses_to_display_unsorted[ch][h]
                   sses_to_display_unsorted[ch][h] = pom

        sses_to_display[ch] = \
           [sses_to_display_unsorted[ch][h] for h in sorted(sses_to_display_unsorted[ch].keys())]

    # Draw connection lines
    for ch in sses_to_display:
        for i in range(1,len(sses_to_display[ch])):
            if chain_id[sses_to_display[ch][i]] not in chain_groups:
                chain_groups[chain_id[sses_to_display[ch][i]]] = img.add_group(chain_id[sses_to_display[ch][i]], opacity=opacity, color="grey", classname="chain")
            if chain_id[sses_to_display[ch][i]] not in chain_groups_links:
                chain_groups_links[chain_id[sses_to_display[ch][i]]] = img.add_group(chain_id[sses_to_display[ch][i]], opacity=opacity, color="#777777", classname="links", ggroup=chain_groups[chain_id[sses_to_display[ch][i]]], without_id=True)
            h1 = sses_to_display[ch][i - 1]
            h2 = sses_to_display[ch][i]
            if (pdb_file.startswith("protein")):
                if (chain_id[h1] != chain_id[h2]):
                    continue
            helic_type = h_type[h1] if h1 in h_type else 3
            if not spike:
                if (h1 in colors) and ((helic_type == 1) or ("protein" in add_descriptions)) or ("multiaf" in add_descriptions):
                    h_color = colors[h1]
                else:
                    h_color = "grey"
            else:
                if colors[h1][0] == "#":
                    i1 = int(colors[h1][1:3], base=16) // 2
                    i2 = int(colors[h1][3:5], base=16) // 2
                    i3 = int(colors[h1][5:7], base=16) // 2
                    h_color = '#%02x%02x%02x' % (i1, i2, i3)
                else:
                    h_color = colors[h1]
            if ('only' not in vars(args)) or ("sses" in args.only):
                h_color = "grey"
                img.helic_connector(
                    coord_from=(layout[h1][0] + ind[0], layout[h1][1] + ind[1]),
                    coord_to=(layout[h2][0] + ind[0], layout[h2][1] + ind[1]),
                    size_from=sizes[h1], size_to=sizes[h2],
                    angle_from=angles[h1], angle_to=angles[h2],
                    color=h_color, opacity=opacity, sse_group=chain_groups_links[chain_id[sses_to_display[ch][i]]])

    # Draw secondary structures
    for ch in sses_to_display:
        for helic in sses_to_display[ch]:
            if chain_id[helic] not in chain_groups:
                chain_groups[chain_id[helic]] = img.add_group(chain_id[helic], opacity=opacity, color=colors[helic], classname="chain")
            if chain_id[helic] not in chain_groups_sses:
                chain_groups_sses[chain_id[helic]] = img.add_group(chain_id[helic], opacity=opacity, color=colors[helic], classname="sses", ggroup=chain_groups[chain_id[helic]], without_id=True)

            coord = (layout[helic][0] + ind[0], layout[helic][1] + ind[1])
            helic_type = h_type[helic] if helic in h_type else 3

            if not spike:
                if (helic in colors) and ((helic_type == 1) or ("protein" in add_descriptions)) or ("multiaf" in add_descriptions):
                    h_color = colors[helic]
                elif (helic in colors) and (helic_type != 1):
                    h_color = "grey"
            else:
                if colors[helic][0] == "#":
                    i1 = int(colors[helic][1:3], base=16) // 2
                    i2 = int(colors[helic][3:5], base=16) // 2
                    i3 = int(colors[helic][5:7], base=16) // 2
                    h_color = '#%02x%02x%02x' % (i1, i2, i3)
                else:
                    h_color = colors[helic]

            if ('only' not in vars(args)) or (("sses" in args.only) or ("charges" in args.only) or ("am" in args.only)):
                if (charges != None):
                    # draw sses based on sse type
                    ccmin = None
                    for cc in charges[helic]:
                        if ccmin == None or (ccmin > int(cc)):
                            ccmin = int(cc)
                    cc = {}
                    # set colors for each res - count sses res from 0
                    for c in charges[helic]:
                        if charges[helic][c] < -40:
                            cc[int(c)-ccmin] = "red"
                        elif charges[helic][c] > 40:
                            cc[int(c)-ccmin] = "blue"
                        else:
                            cc[int(c)-ccmin] = "grey"
                    img.secstrperres(name=helic, coord=coord, size=sizes[helic],
                       angle=angles[helic], h_type=helic_type,
                       color=cc, opacity=opacity, threeddata = threeddata[helic])

                elif (am != None):
                    # draw sses based on sse type
                    ccmin = int(min(am[helic]))
                    cc = {}
                    # set colors for each res - count sses res from 0
                    for c in am[helic]:
                        if am[helic][c] < 35:
                            cc[int(c)-ccmin] = "blue"
                        elif am[helic][c] < 65:
                            cc[int(c)-ccmin] = "grey"
                        else:
                            cc[int(c)-ccmin] = "red"
                    img.secstrperres(name=helic, coord=coord, size=sizes[helic],
                       angle=angles[helic], h_type=helic_type,
                       color=cc, opacity=opacity, threeddata=threeddata[helic])
                else:
                    if ((args.groups is not None) and (args.groups != False)):
                        gr = chain_groups_sses[chain_id[helic]]
                    else:
                        gr = False
                    img.secstr(name=helic+"_"+str(helices_residues[helic][0])+"_"+str(helices_residues[helic][1]),
                        coord=coord, size=sizes[helic], angle=angles[helic], h_type=helic_type,
                        color=h_color, opacity=opacity, threeddata=threeddata[helic], sse_group=gr,
                        is_tm=(helic in tm_helices))
                img_data[helic] = {"size": sizes[helic] * 2 * img.CIRCLE_MULTIPLIER,
                               "color": h_color,
                               "layout": coord,
                               "angles": int(angles[helic]),
                               "residues": helices_residues[helic]}


def show_ligands(img, ind, ligands, opacity, add_descriptions, ligands_data, img_data_ligands, args):
    """ Draw full layout -
        text info if it is necessary
    """
    # Draw ligands
    for l in ligands:
        # if there is a color array found out ligand color
        for i in ligands_data:
            if (ligands[l]["e_type"] == i["e_type"]) and  \
               (ligands[l]["edsc"] == i["edsc"]) and  \
               (ligands[l]["ecmp"] == i["ecmp"]):
                color = i["color"]
                ocolor = i["ocolor"]
                break

        coord = (ligands[l]["layout"][0] + ind[0], ligands[l]["layout"][1] + ind[1])
        # with the given color draw the ligand
        if ligands[l]["ecmp"] == None:
            if "RNA" in ligands[l]["edsc"]:
                ligand_name = "RNA"
            elif "DNA" in ligands[l]["edsc"]:
                ligand_name = "DNA"
            else:
                ligand_name = ""
        else:
            ligand_name = ligands[l]["ecmp"]
        if ligands[l]["edsc"] == None:
            ligand_name_suffix= ""
        else:
            ligand_name_suffix = ligands[l]["edsc"].replace(" ","#").replace("(","[").replace(")","]").replace(",","-")

        if "res" in ligands[l]:
            ligand_group = img.add_group(f"{ligand_name}_{str(ligands[l]['res'])}_{ligands[l]['auth_chain']}_{ligand_name_suffix}", opacity=opacity, color=color, classname="ligand")
        else:
            ligand_group = img.add_group(f"{ligand_name}_x_x_{ligand_name_suffix}", opacity=opacity*5, color=color, classname="ligand")
        img.ligand(name=l, coord=coord, size=ligands[l]["atoms"],
                   e_type=ligands[l]["e_type"], edsc=ligands[l]["edsc"], ecmp=ligand_name,
                   add_descriptions=add_descriptions, input_color=color, ocolor=ocolor, ligand_group=ligand_group)
        # find out the name of ligands withit the digit part (thus either full name, or for DNA/RNA only part
        # without the residuum number - DNA/ RNA name)
        d = -1
        for i, c in enumerate(l):
            if c.isdigit():
                d = i
                break
        if d != -1:
            m_l = l[:d]
        else:
            m_l = l

        img_data_ligands[l] = {"size": img.ligand_diagram_size(ligands[l]["atoms"]),
                               "input_color": color,
                               "ocolor": ocolor,
                               "layout": coord,
                               "multi_ligand": m_l,
                               "name": img.lig_text(edsc=ligands[l]["edsc"], ecmp=ligand_name)}

        # if the ligand is a part of a ligand chain (DNA/RNA and its residues,
        # draw a line between the neighbour
        for h, n in enumerate(l):
            if n.isdigit():
                name = l[:h]
                counter = str(int(l[h:]) + 1)
                newname = name + counter
                if newname in ligands:
                    tocoord = (ligands[newname]["layout"][0] + ind[0],
                               ligands[newname]["layout"][1] + ind[1])
                    img.ligand_connector(coord_from=coord,
                                         coord_to=tocoord,
                                         color=color)
                break


def show_channels(img, ind, channels, opacity, add_descriptions, channels_data, args):
    """
    Draw protein tunnels/channels as smooth variable-width ribbons and return
    tunnel mouth label descriptors to be rendered in the foreground layer.
    """
    if not channels:
        return []

    tunnel_groups = []
    current_prefix = None
    current_keys = []
    for l in channels:
        d = -1
        for i, c in enumerate(l):
            if c.isalpha():
                d = i
                break
        s_l = l[:d] if d != -1 else l
        if current_prefix is None or s_l != current_prefix:
            if current_keys:
                tunnel_groups.append((current_prefix, current_keys))
            current_prefix = s_l
            current_keys = [l]
        else:
            current_keys.append(l)
    if current_keys:
        tunnel_groups.append((current_prefix, current_keys))

    label_items = []
    for idx_t, (s_l, t_keys) in enumerate(tunnel_groups):
        first_l = t_keys[0]
        ch_entry = channels[first_l]
        text_name = str(ch_entry.get("text_names") or "")
        color = ch_entry.get("color")

        for i in channels_data:
            if (ch_entry.get("id") == i.get("id")) and (ch_entry.get("filename") == i.get("filename")):
                if not text_name:
                    text_name = str(i.get("text_names") or "")
                if not color or str(color).lower() in ("none", "grey", "gray", "#665646", "#887766"):
                    color = i.get("color")
                break

        if not color or str(color).lower() in ("none", "grey", "gray", "#665646", "#887766"):
            color = utils_membrane.pick_tunnel_color(text_name, idx=idx_t)

        pts = [(channels[k]["layout"][0] + ind[0], channels[k]["layout"][1] + ind[1]) for k in t_keys]
        rads = [float(channels[k].get("radius", 1.5)) for k in t_keys]

        channel_group = img.add_group(s_l, opacity=opacity, color=color, classname="tunnel")
        img.tunnel_ribbon(pts, rads, color=color, opacity=opacity, channel_group=channel_group)

        if len(pts) >= 2 and text_name:
            p_end = pts[-1]
            p_prev = pts[max(0, len(pts) - 3)]
            dx = p_end[0] - p_prev[0]
            dy = p_end[1] - p_prev[1]
            dn = math.hypot(dx, dy)
            ux, uy = (dx / dn, dy / dn) if dn > 1e-4 else (0.0, 1.0)
            w_end = max(0.24, min(0.75, 0.19 * rads[-1]))
            lbl_pt = [p_end[0] + ux * (w_end + 1.05), p_end[1] + uy * (w_end + 1.05)]
            label_items.append({
                "text": text_name,
                "pos": lbl_pt,
                "anchor": p_end,
                "ux": ux,
                "uy": uy,
                "color": color,
            })

    return label_items


def show_channel_labels(img, label_items, layout, ind, angles, sizes, ligands=None):
    """Render de-overlapped tunnel mouth labels in a foreground SVG group."""
    if not label_items:
        return

    # Collect 2D SSE segments and ligand centers in shifted canvas coordinates
    sse_segs = []
    if layout:
        for s, coord in layout.items():
            ang = math.radians(angles.get(s, 0.0)) if angles else 0.0
            hl = 0.1 * float(sizes.get(s, 8)) if sizes else 0.8
            cx, cy = coord[0] + ind[0], coord[1] + ind[1]
            a = (cx + hl * math.sin(ang), cy + hl * math.cos(ang))
            b = (cx - hl * math.sin(ang), cy - hl * math.cos(ang))
            sse_segs.append((a, b))

    lig_pts = []
    if ligands:
        for lk, lval in ligands.items():
            if isinstance(lval, dict) and lval.get("layout"):
                lig_pts.append((lval["layout"][0] + ind[0], lval["layout"][1] + ind[1]))

    for _ in range(18):
        # 1. Clear SSE segments and ligands
        for item in label_items:
            px, py = item["pos"]
            for (ax, ay), (bx, by) in sse_segs:
                vx, vy = bx - ax, by - ay
                v2 = vx * vx + vy * vy
                t = max(0.0, min(1.0, ((px - ax) * vx + (py - ay) * vy) / (v2 + 1e-9)))
                qx, qy = ax + t * vx, ay + t * vy
                ddx, ddy = px - qx, py - qy
                d = math.hypot(ddx, ddy)
                if d < 1.25:
                    if d > 1e-4:
                        px += (ddx / d) * (1.25 - d) * 0.5
                        py += (ddy / d) * (1.25 - d) * 0.5
                    else:
                        px += item["ux"] * 0.6
                        py += item["uy"] * 0.6

            for lx, ly in lig_pts:
                ddx, ddy = px - lx, py - ly
                d = math.hypot(ddx, ddy)
                if d < 1.75:
                    if d > 1e-4:
                        px += (ddx / d) * (1.75 - d) * 0.5
                        py += (ddy / d) * (1.75 - d) * 0.5
            item["pos"] = [px, py]

        # 2. Pairwise label-label repulsion
        n_lbl = len(label_items)
        for i in range(n_lbl):
            for j in range(i + 1, n_lbl):
                p1 = label_items[i]["pos"]
                p2 = label_items[j]["pos"]
                ddx = p1[0] - p2[0]
                ddy = p1[1] - p2[1]
                if abs(ddx) < 3.6 and abs(ddy) < 1.25:
                    push_y = (1.25 - abs(ddy)) * 0.45 * (1.0 if ddy >= 0 else -1.0)
                    push_x = (3.6 - abs(ddx)) * 0.25 * (1.0 if ddx >= 0 else -1.0)
                    p1[0] += push_x
                    p1[1] += push_y
                    p2[0] -= push_x
                    p2[1] -= push_y

        # 3. Tether label near its tunnel mouth
        for item in label_items:
            ax, ay = item["anchor"]
            ddx = item["pos"][0] - ax
            ddy = item["pos"][1] - ay
            d = math.hypot(ddx, ddy)
            if d > 2.6:
                item["pos"][0] = ax + (ddx / d) * 2.6
                item["pos"][1] = ay + (ddy / d) * 2.6

    lbl_group = img.dwg.add(img.dwg.g(id="tunnel_labels", class_="tunnel_labels"))
    for item in label_items:
        img.tunnel_label(item["text"], item["pos"], color=item["color"], group=lbl_group)


def show_template_layout(img, layout, ind, h_type, sizes, colors, angle):

    # Filter only SSEs we are going to display and sort them based on res. numbers
    sses_to_display = []
    for h in layout.keys():
        if (utils_sses_types.sse_should_be_shown(h, h_type, sizes) == 0):
            sses_to_display.append(h)

    # show template symbol for all SSES that showld be visible
    for helic in sses_to_display:
        coord = (layout[helic][0] + ind[0], layout[helic][1] + ind[1])

        if helic not in colors:
            h_color = 'grey'
        else:
            if colors[helic][0] == "s":
                color_tuple = pymol.cmd.get_color_tuple(colors[helic])
                color_tuple_p = (255 - int(color_tuple[0] * 255),
                                 255 - int(color_tuple[1] * 255),
                                 255 - int(color_tuple[2] * 255))
                h_color = '#%02x%02x%02x' % color_tuple_p
            else:
                h_color = colors[helic]
        h_color = "black"

        img.templsecstr(helic, coord=coord, size=sizes[helic],
                        h_type=h_type[helic], color=h_color, angle=angle[helic])


def stretch_dimensions_by_layout(l, width_left_max, width_right_max,
                                 height_up_max, height_down_max, avg):
    """
    Given loader.Layout object instance and current values for canvas
    dimensions, stretch these dimensions if necessary and return them
    """
    l_w_h, l_c = l.get_canvas_size(utils_sses_types.sse_should_be_shown, avg)
    if (width_left_max == "p") and (width_right_max == "p") and \
       (height_up_max == "p") and (height_down_max == "p"):
        width_left_max = l_c[0]
        width_right_max = l_w_h[0] + l_c[0] + 5
        height_down_max = l_c[1]
        height_up_max = l_w_h[1] + l_c[1]
    else:
        if width_left_max > l_c[0]:
            width_left_max = l_c[0]
        if width_right_max < (l_w_h[0] + l_c[0] + 5):
            width_right_max = l_w_h[0] + l_c[0] + 5
        if height_down_max > l_c[1]:
            height_down_max = l_c[1]
        if height_up_max < (l_w_h[1] + l_c[1]):
            height_up_max = l_w_h[1] + l_c[1]

    return width_left_max, width_right_max, height_up_max, height_down_max


def stretch_dimensions_by_ligands(l, width_left_max, width_right_max,
                                  height_up_max, height_down_max, avg):
    """
    Given ligands object instance and current values for canvas
    dimensions, stretch these dimensions if necessary and return them
    """
    c = [3, 3]
    if (width_left_max == "p") and (width_right_max == "p") and \
       (height_up_max == "p") and (height_down_max == "p"):
        width_left_max = l[0] - c[0]
        width_right_max = l[0] + c[0]
        height_down_max = l[0] - c[1]
        height_up_max = l[1] + c[1]
    else:
        if width_left_max > l[0] - c[0]:
            width_left_max = l[0] - c[0]
        if width_right_max < (l[0] + c[0]):
            width_right_max = l[0] + c[0]
        if height_down_max > l[1] - c[1]:
            height_down_max = l[1] - c[1]
        if height_up_max < (l[1] + c[1]):
            height_up_max = l[1] + c[1]

    return width_left_max, width_right_max, height_up_max, height_down_max


def rotate_and_translate_layout(
        rot_trans, orotation, irotation, layout_call_lists, layout_name_lists, l_all_loaded, data):
    # Rotate and translate all layouts which should be shown

    i = 0
    while i < len(layout_call_lists):
        f_layout = {}
        if irotation:
            jname = layout_name_lists[i]
            jnames = jname.split("/")
            j = jnames[len(jnames) - 1]
        else:
            j = i

        r, t = rot_trans[j]
        t_layout = superimpose2d.rotate_and_translate(
            layout_call_lists[i], r, t)

        # rotate SSEs angles as well (not only their positions
        l_all_loaded[i].data['layout'] = f_layout
        if (rot_trans[j][0][0][0] > 1):
            rot_trans[j][0][0][0] = 1
        t_angle = math.acos(rot_trans[j][0][0][0]) * 180 / math.pi
        # The sign depends on [[0][1]] member of matrix. If it is <0 the angle is - as well
        if rot_trans[j][0][0][1] < 0:
            t_angle = - t_angle
        for h in l_all_loaded[i].data['angles']:
            l_all_loaded[i].data['angles'][h] = l_all_loaded[i].data['angles'][h] + t_angle
            if (abs(rot_trans[j][0][0][0] + rot_trans[j][0][1][1])) < 0.001:
                l_all_loaded[i].data['angles'][h] = 180 - \
                    l_all_loaded[i].data['angles'][h]
            if l_all_loaded[i].data['angles'][h] > 360:
                l_all_loaded[i].data['angles'][h] = l_all_loaded[i].data['angles'][h] - 360

        # rotate ligands
        if "ligands" in l_all_loaded[i].data:
            l_layout_old = {}
            for l in l_all_loaded[i].data["ligands"]:
                l_layout_old[l] = l_all_loaded[i].data["ligands"][l]["layout"]

            l_layout = superimpose2d.rotate_and_translate(
                l_layout_old, r, t)

            for l in l_all_loaded[i].data["ligands"]:
                l_all_loaded[i].data["ligands"][l]["layout"] = l_layout[l]

        # rotate channels
        if "channels" in l_all_loaded[i].data:
            l_layout_old = {}
            for l in l_all_loaded[i].data["channels"]:
                l_layout_old[l] = l_all_loaded[i].data["channels"][l]["layout"]

            l_layout = superimpose2d.rotate_and_translate(
                l_layout_old, r, t)

            for l in l_all_loaded[i].data["channels"]:
                l_all_loaded[i].data["channels"][l]["layout"] = l_layout[l]

        # add rotation record to data structure - to be stored later
        if orotation:
            l_name = layout_name_lists[i].split("/")
            data[l_name[len(l_name) - 1]] = ([list(r[0]),
                                              list(r[1])], (list(t[0]), list(t[1])))

        for nn in t_layout:
            name_number = nn[:nn.rfind("_")]
            if (name_number not in f_layout):
                f_layout[name_number] = t_layout[nn]
        i = i + 1

    return data


def stretch_dimensions(l_all_loaded, data_ligands, ligands_set, ligands_data):
    # Set the appropriate canvas size based on l_all_loaded layouts
    width_left_max, width_right_max, height_up_max, height_down_max = "p", "p", "p", "p"

    for l in l_all_loaded:
        width_left_max, width_right_max, height_up_max, height_down_max = \
            stretch_dimensions_by_layout(l, width_left_max, width_right_max,
                                         height_up_max, height_down_max, avg=1)
        if "channels" in l.data:
            for ch_k in l.data["channels"]:
                ch_pt = l.data["channels"][ch_k].get("layout")
                if ch_pt:
                    width_left_max, width_right_max, height_up_max, height_down_max = \
                        stretch_dimensions_by_ligands(ch_pt, width_left_max, width_right_max,
                                                      height_up_max, height_down_max, avg=1)
        if "ligands" in l.data and l.data["ligands"]:
            for lg_k in l.data["ligands"]:
                lg_pt = l.data["ligands"][lg_k].get("layout")
                if lg_pt:
                    width_left_max, width_right_max, height_up_max, height_down_max = \
                        stretch_dimensions_by_ligands(lg_pt, width_left_max, width_right_max,
                                                      height_up_max, height_down_max, avg=1)
        if "ghost_sses" in l.data and l.data["ghost_sses"]:
            for gh in l.data["ghost_sses"]:
                gh_pt = gh.get("coord")
                if gh_pt:
                    width_left_max, width_right_max, height_up_max, height_down_max = \
                        stretch_dimensions_by_ligands(gh_pt, width_left_max, width_right_max,
                                                      height_up_max, height_down_max, avg=1)

    # Set the appropriate canvas size based on ligands layouts
    for l in data_ligands:
        width_left_max, width_right_max, height_up_max, height_down_max = \
            stretch_dimensions_by_ligands(data_ligands[l]["layout"], width_left_max, width_right_max,
                                          height_up_max, height_down_max, avg=1)


    # TODO should not be used ?? v
#    for i in ligands_set:
#        for l in i.data["layout"]:
#            width_left_max, width_right_max, height_up_max, height_down_max = \
#                stretch_dimensions_by_ligands(i.data["layout"][l], width_left_max, width_right_max,
#                                              height_up_max, height_down_max, avg=1)

    ind = [0, 0]
    if height_up_max == "p":
        # nothing to show
        dimensions = ((20, 20), (-10, -10))
    else:
        # something should be shown
        height_up_max = height_up_max + 6
        width = width_right_max - width_left_max
        height = height_up_max - height_down_max

        # if we draw ligand_legend the image have to be big enough for it
        ligand_legend_width = min( 8 * 11, (len(ligands_data) // 5) * 11) + 3
        if width < ligand_legend_width:
           width = ligand_legend_width
           width_left_max -= (ligand_legend_width -width) / 2

        if (width < MIN_CANVAS_DIMENSIONS[0]):
            ind[0] = (MIN_CANVAS_DIMENSIONS[0] - width) / 2
        if (height < MIN_CANVAS_DIMENSIONS[1]):
            ind[1] = (MIN_CANVAS_DIMENSIONS[1] - height) / 2
        dimensions = ((width, height), (-width_left_max, -height_down_max))

    return dimensions, ind


def set_opacity(transparency, number):
    # set opacity value
    if transparency == True:
        if number == 0:
            return 1
        else:
            if number > 10000:
                return math.sqrt(number) * 2 * 2.6 * math.sqrt(number/10000) / number
            else:
#                return math.sqrt(number) * 2 * 2.6 / number
                return math.sqrt(number) *2 / number
    else:
        return 1


def draw_nonempty_layouts_set(l_all_loaded, otemplate, opacity, img, template,
                              number_template, ind, start, ratio,
                              add_descriptions, img_data, spike,
                              ligands, ligands_data, img_data_ligands,
                              channels_data, args):

    c = 0
    for l in l_all_loaded:
        if ((c == 0) and (start is False)) or ((c >= ratio * len(l_all_loaded)) and (otemplate == 0)):
            # The first layout should be shown only once
            # This should not be shown and preprocessed to template
            logging.debug("%s should be ommited %s %s %s %s" % (l.filename, c, len(l_all_loaded), c % len(l_all_loaded), ratio))
            c = c + 1
            continue
        # Filter only SSEs we are going to display and sort them based on res. numbers
        sses_to_display = []
        for h in l.data['layout'].keys():
            if (utils_sses_types.sse_should_be_shown(
                    h, l.data['h_type'], l.data['size']) == 0):
                sses_to_display.append(h)
        sses_first_level = []
        for h in l.data['layout'].keys():
            if (utils_sses_types.sse_first_level(
                    h, l.data['h_type'], l.data['size']) == 0):
                sses_first_level.append(h)

        # Snap SSE angles to multiples of angle_step degrees (if enabled) and fix layout conflicts
        _mem_cfg = l.data.get('membrane')
        _step_cfg = _mem_cfg.get('angle_step', utils_membrane.ANGLE_STEP) if isinstance(_mem_cfg, dict) else utils_membrane.ANGLE_STEP
        l.data['angles'] = utils_membrane.snap_angles_dict(
            l.data['angles'], beta_connectivity=l.data.get('beta_connectivity'), step=_step_cfg
        )
        layout = resolve_sses_conflicts(l.data['layout'], l.data['angles'], l.data['size'],
                                        sses_to_display, l.data['beta_connectivity'])
        ghosts = []
        if len(l_all_loaded) == 1 and not getattr(args, "multi", False) and opacity == 1 and "3dlayout" in l.data:
            ghosts = utils_membrane.compute_ghost_sses(
                layout_2d={h: layout[h] for h in sses_to_display if h in layout},
                angles_2d=l.data["angles"],
                sizes=l.data["size"],
                vectors_3d=l.data["3dlayout"],
                h_type=l.data.get("h_type"),
                beta_connectivity=l.data.get("beta_connectivity"),
                mem_info=l.data.get("membrane"),
                helices_residues=l.data.get("helices_residues"),
                ligands=l.data.get("ligands") if ligands else None,
            )
        if ligands and "ligands" in l.data:
            layout_with_ghosts = {h: layout[h] for h in sses_to_display if h in layout}
            angles_with_ghosts = dict(l.data["angles"])
            sizes_with_ghosts = dict(l.data["size"])
            for idx_g, gh in enumerate(ghosts):
                gk = f"__ghost_{idx_g}_{gh['source_sse']}"
                layout_with_ghosts[gk] = gh["coord"]
                angles_with_ghosts[gk] = gh["angle"]
                sizes_with_ghosts[gk] = gh["size"]
            utils_membrane.deoverlap_ligands_from_sses(
                l.data["ligands"],
                layout_with_ghosts,
                angles_with_ghosts,
                sizes_with_ghosts,
                l.data.get("membrane"),
            )
        if "channels" in l.data and l.data["channels"]:
            utils_membrane.localize_channels_2d(
                l.data["channels"],
                layout_2d={h: layout[h] for h in sses_to_display if h in layout},
                angles_2d=l.data["angles"],
                sizes=l.data["size"],
                vectors_3d=l.data.get("3dlayout"),
                h_type=l.data.get("h_type"),
                ligands=l.data.get("ligands"),
                mem_info=l.data.get("membrane"),
            )
        if (c < ratio * len(l_all_loaded)):
            try:
                tunnel_labels = []
                if c == 0 and isinstance(l.data.get("membrane"), dict) and l.data["membrane"].get("is_transmembrane", True) and (('only' not in vars(args)) or ("sses" in args.only)):
                    mem = l.data["membrane"]
                    img.membrane_band(
                        y_center=mem.get("y_center_2d", 0.0) + ind[1],
                        half_height=mem.get("band_half_height_2d", 3.0),
                    )
                if ("channels" in l.data) and (('only' not in vars(args))  or ("channels" in args.only)):
                    tunnel_labels = show_channels(img, ind, l.data["channels"], opacity, add_descriptions, channels_data, args)
                logging.debug("Drawing full %s" % l.filename)
                if 'chain' in l.data.keys():
                    chain = l.data['chain']
                else:
                    chain = ""
                if ('only' not in vars(args)) or ("sses" in args.only) or ("template" in args.only):
                    if 'charges' in l.data:
                        charges = l.data['charges']
                    else:
                        charges = None

                    if 'am' in args and 'am' in l.data:
                        am = l.data['am']
                    else:
                        am = None

                    show_full_layout(img, layout, ind,
                                 l.data['helices_residues'], l.data['angles'], l.data['size'],
                                 l.data['protein'], l.data['chain_id'], chain, l.data['domain'],
                                 l.data['h_type'], l.data['color'], add_descriptions, opacity,
                                 img_data, spike, ligands, ligands_data, img_data_ligands, args,
                                 l.data['3dlayout'], charges, am, membrane=l.data.get("membrane"))

                    if ghosts:
                        for gh in ghosts:
                            src_s = gh["source_sse"]
                            part_s = gh.get("partner_sse")
                            g_color = l.data["color"].get(src_s, "#6b7682")
                            g_coord = (gh["coord"][0] + ind[0], gh["coord"][1] + ind[1])
                            p_coord = (layout[part_s][0] + ind[0], layout[part_s][1] + ind[1]) if part_s in layout else None
                            p_size = l.data["size"].get(part_s) if part_s else None
                            p_angle = l.data["angles"].get(part_s) if part_s else None
                            orig_c = gh.get("origin_coord")
                            orig_coord = (orig_c[0] + ind[0], orig_c[1] + ind[1]) if orig_c else None
                            img.ghost_secstr(
                                name=src_s,
                                coord=g_coord,
                                size=gh["size"],
                                angle=gh["angle"],
                                color=g_color,
                                partner_coord=p_coord,
                                partner_size=p_size,
                                partner_angle=p_angle,
                                origin_coord=orig_coord,
                                show_arc=bool(gh.get("show_arc", False)),
                            )

                if ligands and "ligands" in l.data and (('only' not in vars(args)) or ("ligands" in args.only)):
                    show_ligands(img, ind, l.data["ligands"], opacity, add_descriptions, ligands_data, img_data_ligands, args)

                if tunnel_labels and (not hasattr(args, "multi") or not args.multi):
                    show_channel_labels(
                        img,
                        tunnel_labels,
                        {h: layout[h] for h in sses_to_display if h in layout},
                        ind,
                        l.data["angles"],
                        l.data["size"],
                        l.data.get("ligands"),
                    )


            except Exception:
                raise
        if (otemplate != None) and (otemplate != "None") and (otemplate != 0) and (('only' not in vars(args)) or ("template" in args.only)):
            # create template layout
            for i in (layout):
                if i not in template['layout']:
                    template['layout'][i] = (layout[i][0], layout[i][1])
                    if i not in l.data['h_type']:
                        template['h_type'][i] = 3
                    else:
                        template['h_type'][i] = l.data['h_type'][i]
                    template['size'][i] = l.data['size'][i]
                    if (template['h_type'][i] == 1):
                        if (i in l.data['color']):
                            template['color'][i] = l.data['color'][i]
                        else:
                            template['h_type'][i] = 3
                            template['color'][i] = 'grey'
                    if i[0] == "E":
                        # 1. value of angle should be between 0 .. 360 and its the model of angle
                        # 2. value is the sum of all differences of strands angles in the same direction as the first layout
                        # 3. value is the number of all strands in the same direction as the first layout
                        # 4. second is the sum of all differences of strands angles in the opposit direction as the first layout
                        # 5. is the number of all strands in the opposit direction as the first layout
                        # find out whether there already is not a strand fromthe same sheet which direction is set
                        # (prevent nonparalell directions in one sheet if a strand is not annotated properly in the first protein
                        beta_connectivity = l.data['beta_connectivity']
                        found = False
                        for bc in beta_connectivity:
                            # if there is a strand from the same sheet and it has set an angle, 
                            # set the same angle to this strand
                            sts = bc.split(":")
                            if i in sts:
                                for j in sts:
                                    if i != j and j in template['angle']:
                                        if abs((l.data['angles'][i] -  l.data['angles'][j]) % 360) < 10:
                                            template['angle'][i] = (
                                                template['angle'][j][0], 0, 1, 0, 0)
                                        else:
                                            template['angle'][i] = (
                                                (template['angle'][j][0]+180) % 360, 0, 1, 0, 0)
                                        found = True
                                        break
                        if not found:
                            # there is no strand in sheet with set angle, thus set this as first
                            template['angle'][i] = (
                                l.data['angles'][i] % 360, 0, 1, 0, 0)
                    else:
                        template['angle'][i] = (
                            l.data['angles'][i], 0, 1, 0, 0)
                    number_template[i] = 1
                else:
                    template['layout'][i] = (layout[i][0] + template['layout'][i][0],
                                             layout[i][1] + template['layout'][i][1])
                    template['size'][i] = template['size'][i] + \
                        l.data['size'][i]
                    if i[0] == "E":
                        c1 = template['angle'][i][0]
                        c2 = l.data['angles'][i]
                        if (((c2 - c1) % 360) < 20):
                            # in the same direction
                            template['angle'][i] = (template['angle'][i][0],
                                                    template['angle'][i][1] + ((c2 - c1) % 360), template['angle'][i][2] + 1,
                                                    template['angle'][i][3], template['angle'][i][4])
                        elif (((c1 - c2) % 360) < 20):
                            # in the same direction
                            template['angle'][i] = (template['angle'][i][0],
                                                    template['angle'][i][1] - ((c1 - c2) % 360), template['angle'][i][2] + 1,
                                                    template['angle'][i][3], template['angle'][i][4])
                        elif (((c2 - c1 + 180) % 360) < 20):
                            # in the opposite direction
                            template['angle'][i] = (template['angle'][i][0],
                                                    template['angle'][i][1], template['angle'][i][2],
                                                    template['angle'][i][3] + ((c2 - c1 + 180) % 360), template['angle'][i][4] + 1)
                        elif (((c1 - c2 + 180) % 360) < 20):
                            # in the opposite direction
                            template['angle'][i] = (template['angle'][i][0],
                                                    template['angle'][i][1], template['angle'][i][2],
                                                    template['angle'][i][3] - ((c1 - c2 + 180) % 360), template['angle'][i][4] + 1)
                    else:
                        # the SSE is helic - jsut copy the angle
                        template['angle'][i] = (
                            template['angle'][i][0],
                            template['angle'][i][1] + l.data['angles'][i] - template['angle'][i][0],
                            template['angle'][i][2] + 1,
                            template['angle'][i][3],
                            template['angle'][i][4])
                    number_template[i] = number_template[i] + 1
        c = c + 1
    return template


def draw_empty_layouts_sets(l, img, l_empty_loaded, output, l_all_loaded, args, ligands, opacity):
    # draw diagram which does not contain any chain

    # we have to set diagram dimensions
    if 'ligands' in l.data:
        lig = l.data['ligands']
    else:
        lig = []

    dimensions, ind = stretch_dimensions(l_all_loaded, lig, {}, {})

    if "fill" in args:
        fill = args.fill
    else:
        fill = True

    img = l.get_canvas(
        filename_out=output,
        min_dimensions=MIN_CANVAS_DIMENSIONS,
        dimensions=dimensions, fill=fill)

    # if ligannd should be shown, show them
    if args.ligands and "ligands" in l.data:
        if ligands:
            show_ligands(img, ind, l.data["ligands"], opacity, args.add_descriptions, ligands, {}, args)

    return img, ind


def rot_trans_sync(rot_trans, start_set):

    # Find necessary rotation
    new_inv = np.linalg.inv(np.array(rot_trans[0][0]))
    transf_m = np.dot(start_set[0], new_inv)

    # Transform all input layouts
    c = 0
    for i in rot_trans:
        # rti = rot_trans[i]
        trans = [i[1][0][0], i[1][1][0]]
        rt = (np.dot(transf_m, np.array(i[0])),
              np.dot(transf_m, [[trans[0]], [trans[1]]]))
        rot_trans[c] = rt
        # t = np.dot(transf_m,np.array(i[1]))
        # rot_trans[c]= [np.dot(transf_m, np.array(i[0])), [0,0]]
        c = c + 1
    return rot_trans


def process_layout_list(l, layout_all_lists, layout_call_lists,
                        layout_name_lists, start_set, l_all_loaded,
                        otemplate, opacity, img, template, number_template,
                        ind, data, ratio, img_data, args,
                        ligands_data, img_data_ligands, channels_data):
    # set rotation based on either input file or on multi_superimpose fce
    # TODO -  rot_trand is read here from irotaion, but it is still set to t_a_t_l function above
    if (args.irotation):
        if (start_set == 0):
            with open(args.irotation, "r") as fd:
                rot_trans = json.load(fd)
    else:
        layout_short_list = []
        if "spike_protein" in args:
            for i in layout_all_lists:
                l1 = {}
                for j in i.keys():
                    if int(j[2]) < 4:
                        l1[j] = i[j]
                layout_short_list.append(l1)
            rot_trans = superimpose2d.multi_superimpose(
                layout_short_list, allow_mirror=True)
        elif any(ld.data.get("membrane") for ld in l_all_loaded):
            # Keep membrane-aligned orientation horizontal while superimposing family members
            rot_trans = utils_membrane.superimpose_membrane_layouts(layout_all_lists)
        else:
            rot_trans = superimpose2d.multi_superimpose(
                layout_all_lists, allow_mirror=True)

        if (start_set != 0):
            rot_trans = rot_trans_sync(rot_trans, start_set)
            dif = (start_set[1][0] - rot_trans[0][1][0],
                   start_set[1][1] - rot_trans[0][1][1])
            for i in rot_trans:
                i[1][0] = i[1][0] + dif[0]
                i[1][1] = i[1][1] + dif[1]

    # rotate input layout
    data = rotate_and_translate_layout(rot_trans, args.orotation, args.irotation,
                                       layout_call_lists, layout_name_lists,
                                       l_all_loaded, data)

    if args.ligands:
        ligands = True
    else:
        ligands = False


    dimensions = None
    if (start_set == 0):
        # Pre-localize channels and ghost SSEs on the rotated layout so canvas bounds include exterior elements
        for l_pre in l_all_loaded[:1]:
            if ("channels" in l_pre.data and l_pre.data["channels"]) or (len(l_all_loaded) == 1 and not getattr(args, "multi", False) and "3dlayout" in l_pre.data):
                vis_pre = [
                    h for h in l_pre.data["layout"].keys()
                    if utils_sses_types.sse_should_be_shown(h, l_pre.data["h_type"], l_pre.data["size"]) == 0
                ]
                _mem_pre = l_pre.data.get("membrane")
                _step_pre = _mem_pre.get("angle_step", utils_membrane.ANGLE_STEP) if isinstance(_mem_pre, dict) else utils_membrane.ANGLE_STEP
                ang_pre = utils_membrane.snap_angles_dict(
                    dict(l_pre.data["angles"]), beta_connectivity=l_pre.data.get("beta_connectivity"), step=_step_pre
                )
                lay_pre = resolve_sses_conflicts(
                    dict(l_pre.data["layout"]), ang_pre, l_pre.data["size"], vis_pre, l_pre.data["beta_connectivity"]
                )
                if ligands and "ligands" in l_pre.data:
                    utils_membrane.deoverlap_ligands_from_sses(
                        l_pre.data["ligands"],
                        {h: lay_pre[h] for h in vis_pre if h in lay_pre},
                        ang_pre,
                        l_pre.data["size"],
                        l_pre.data.get("membrane"),
                    )
                if "channels" in l_pre.data and l_pre.data["channels"]:
                    utils_membrane.localize_channels_2d(
                        l_pre.data["channels"],
                        layout_2d={h: lay_pre[h] for h in vis_pre if h in lay_pre},
                        angles_2d=ang_pre,
                        sizes=l_pre.data["size"],
                        vectors_3d=l_pre.data.get("3dlayout"),
                        h_type=l_pre.data.get("h_type"),
                        ligands=l_pre.data.get("ligands"),
                        mem_info=l_pre.data.get("membrane"),
                    )
                if len(l_all_loaded) == 1 and not getattr(args, "multi", False) and "3dlayout" in l_pre.data:
                    l_pre.data["ghost_sses"] = utils_membrane.compute_ghost_sses(
                        layout_2d={h: lay_pre[h] for h in vis_pre if h in lay_pre},
                        angles_2d=ang_pre,
                        sizes=l_pre.data["size"],
                        vectors_3d=l_pre.data["3dlayout"],
                        h_type=l_pre.data.get("h_type"),
                        beta_connectivity=l_pre.data.get("beta_connectivity"),
                        mem_info=l_pre.data.get("membrane"),
                        helices_residues=l_pre.data.get("helices_residues"),
                        ligands=l_pre.data.get("ligands") if ligands else None,
                    )
        dimensions, ind = stretch_dimensions(l_all_loaded, {}, l_all_loaded, ligands_data)
        if ligands and ("ligands" in args.add_descriptions):
            if len(ligands_data) > 5:
                dimensions = ((dimensions[0][0], dimensions[0][1] + 5 * 2),
                              (dimensions[1][0], dimensions[1][1] + 5 * 2))
            else:
                dimensions = ((dimensions[0][0], dimensions[0][1] + len(ligands_data) * 2),
                              (dimensions[1][0], dimensions[1][1] + len(ligands_data) * 2))
        if "fill" in args:
            fill = args.fill
        else:
            fill = True

        img = l.get_canvas(filename_out=args.output,
                           min_dimensions=MIN_CANVAS_DIMENSIONS,
                           dimensions=dimensions, fill=fill)
        if not (args.irotation):
            start_set = rot_trans[0]
        start = True
    else:
        start = False


    if "spike_protein" in args:
        spike = True
    else:
        spike = False

    if dimensions and ligands and ("ligands" in args.add_descriptions):
        # if we have some ligands and we should show legend
        write_ligands_legend(img, ligands_data, ligands, dimensions, args)
    template = draw_nonempty_layouts_set(l_all_loaded, otemplate, opacity, img,
                                         template, number_template, ind, start,
                                         ratio, args.add_descriptions, img_data,
                                         spike, ligands, ligands_data, img_data_ligands,
                                         channels_data, args)

    return img, template, ind, start_set, data


def output_template_layout(template, number_template, otemplate):
        for i in template['layout'].keys():
            template['layout'][i] = (template['layout'][i][0] / number_template[i],
                                     template['layout'][i][1] / number_template[i])
            template['size'][i] = template['size'][i] / number_template[i]
            if (template['angle'][i][2] > template['angle'][i][4]):
                template['angle'][i] = template['angle'][i][0] + \
                    (template['angle'][i][1] / template['angle'][i][2])
            else:
                template['angle'][i] = template['angle'][i][0] + \
                    (template['angle'][i][3] / template['angle'][i][4])
        filetempl = open(otemplate, "w")
        json.dump(template, filetempl, indent=4, separators=(',', ':'))
        filetempl.close()


def output_ligands(ligands_data, oligands):
    fd = open(oligands, "w")
    json.dump(ligands_data, fd, indent=4, separators=(',', ':'))
    fd.close()


def compute(args):
    ''' The main function, parse arguments read imput layouts and initializate
        variables
    '''
    ###
    # input arguments:
    # mandatory: --------------------------
    # args.output <file_name> - name of the returned svg diagram (e.g. 1tqn.svg)
    # args.layouts <file_name|file_dir>
    #                         - input layout file or the directory with layout files
    #                           test whether filename or filedir is set is based on
    #                           multi set in description TODO
    #
    # optional: ---------------------------
    # diagram style
    # args.multi <True|False>    - False (default) - means we want single diagram
    #                            - True means we want multi diagram
    #
    # arguments that defined result diagram
    # args.fill <True|False>     - True (default) means the backgroud of diagram is white False means
    #                            the backgroud is transparent
    # args.transparency <True|False>
    #                            - False (default) means no transparency standart for single images
    #                            if True SSEs and tunnels (not ligands) will be transparent (useful for
    #                            multiple diagrams)
    # args.dominik <f|ch>        - use dominik's color mapping of ligands - either functional or chemical
    #                            version default is False - not used
    # args.iligands <file_name>  - have to be used in single diagrams with ligands to have the same color
    #                            scheme for ligands in the whole family
    # args.irotation <file_name> - load rotation and translation of single domains in the diagram (to
    #                            createte comparable single and multi diagrams with given multi one)
    # args.add_descriptions <string>
    #                            - ligand|multi|chain|protein
    #                            - deals with text labels
    #                              - ligand  -
    #                              - multiaf -
    #                              - protein -
    #
    # arguments that defined auxiliary output files
    # args.data <file_name>      - returns data regarding the output diagam (positions of SSEs, ligands
    #                            and metadata (upper_right and down_left points) - used in
    #                            https://sestra.ncbr.muni.cz/ and https://moleonline.biodata.ceitec.cz.
    #                            Second version of 2dprot svg is used directly
    # args.oligands <file_name>  - save ligand colors, output file can be read by args.iligands
    #                            have to be done in multi picture and use in single picture
    # args.orotation <file_name> - save rotation and translation of single domains in the diagram (to
    #                            createte comparable single and multi diagrams with given multi one)
    #
    # obsolete: ---------------------------------
    # args.itemplate <file_name> - if args.only is not set or args,.only is "template"
    #                              output the template which is defined by input file - TODO remove
    # args.otemplate <file_name> - save template to otemplate file - TODO remove
    ###
    #
    # args.groups           - 2
    # args.ligands          - 4
    # args.only             - 4"channels", "ligands", "template", "sse"

    # set defaut arguments
    if "dominik" not in args:
        args.dominik = False
    else:
        if (args.dominik != "f") and (args.dominik != "ch") and (args.dominik is not False):
            print("argument dominik have to be either 'f' or 'ch', it is:", args.dominik)
            print("args", args)
            sys.exit()

    if "multi" not in args:
        args.multi = False
    else:
        if args.multi is not True and args.multi is not False:
            print("argument dominik have to be either 'True or 'False', it is:", args.multi)
            print("args", args)
            sys.exit()

    logging.debug("Generating image from 2D layout file(s) %s into %s"
                  % (','.join(args.layouts), args.output or '[default]'))

    if args.orotation:
        logging.debug("Used output rotation file %s" % args.orotation)
    if args.irotation:
        logging.debug("Set input rotation file %s" % args.irotation)

    # If it is multi image, set only the directory (otherwise there can be too much arguments)
    if (args.multi):
        layout_list = []
        for dfile in os.listdir(args.layouts[0]):
            if dfile.startswith("layout-") and dfile.endswith(".json") and ("_" not in dfile):
                layout_list.append(args.layouts[0] + "/" + dfile)
    else:
        layout_list = args.layouts

    opacity = set_opacity(args.transparency, len(layout_list))
    # Load layouts(s) and create canvas with proper size
    l_all_loaded = []
    layout_all_lists = []
    layout_call_lists = []
    layout_name_lists = []
    l_empty_loaded = []

    if args.otemplate:
        template = {'layout': {}, 'h_type': {}, 'size': {},
                    'color': {}, 'angle': {}}
        otemplate = args.otemplate
    else:
        template = {}
        otemplate = 0
    number_template = {}

    # Read all input layouts - will be used for superimpose and shown as well
    start_set = 0
    data = {}
    img_data = {}
    img_data_ligands = {}
    img = None
    ind = 0
    ligands_data_pre = []
    channels_data_pre = []

    if (len(layout_list) > 100):
        ratio = 1000 / len(layout_list)
    else:
        ratio = 1000

    for filename in layout_list:
        # Read all layout files
        logging.debug("Loading full layout %s" % filename)
        l = loader.Layout(filename=filename)

        # Automatically detect 3D transmembrane bundle & align to horizontal membrane if 3dlayout is present
        if "3dlayout" in l.data and "membrane" not in l.data:
            l.data["layout"], l.data["angles"], mem_info = utils_membrane.align_layout_to_membrane(
                l.data["layout"],
                l.data["angles"],
                l.data.get("size", {}),
                l.data.get("beta_connectivity", []),
                l.data["3dlayout"],
                helices_residues=l.data.get("helices_residues"),
                h_type=l.data.get("h_type"),
                ligands=l.data.get("ligands"),
                channels=l.data.get("channels"),
            )
            if mem_info:
                l.data["membrane"] = mem_info
        else:
            _mem_l = l.data.get("membrane")
            _step_l = _mem_l.get("angle_step", utils_membrane.ANGLE_STEP) if isinstance(_mem_l, dict) else utils_membrane.ANGLE_STEP
            l.data["angles"] = utils_membrane.snap_angles_dict(
                l.data["angles"], beta_connectivity=l.data.get("beta_connectivity"), step=_step_l
            )

        l_list = {}
        cl_list = {}
        prime = 0

        if (args.multi):
            # TODO: is it necesarry
            for i in (l.data['h_type'].keys()):
                if (l.data['h_type'][i] == 1):
                    prime = prime + 1

        # reads per ligands stuff, if it is necessary
        if 'ligands' in l.data:
            for li in l.data['ligands']:
                ecmp = l.data['ligands'][li]["ecmp"]
                edsc = l.data['ligands'][li]["edsc"]
                etype = l.data['ligands'][li]["e_type"]
                if "color" in l.data['ligands'][li]:
                    color = l.data['ligands'][li]["color"]
                else:
                    color = None
                found = False
                for lid in ligands_data_pre:
                    if lid["edsc"] == edsc and lid["e_type"] == etype and lid["ecmp"] == ecmp:
                        lid["count"] += 1
                        found = True
                        break

                if not found:
                    if "res" in l.data['ligands'][li]:
                        res = l.data['ligands'][li]["res"]
                    else:
                        res = -1
                    ligand_d = {"edsc": edsc,
                                "e_type": etype,
                                "ecmp": ecmp,
                                "size": l.data['ligands'][li]["atoms"],
                                "res": res,
                                "count": 1,
                                "color": color,
                                "ocolor": color}
                    ligands_data_pre.append(ligand_d)

        # reads per chanels stuff
        if 'channels' in l.data:
            for cha in l.data['channels']:
                # we have to add filename tag to channel points
                # to prevent mixing of channel points from the different
                # protein/domain files
                l.data["channels"][cha]["filename"] = os.path.basename(filename)

                ch_id = l.data['channels'][cha]["id"]
                found = False
                for chan in channels_data_pre:
                    if chan["id"] == ch_id:
                        chan["count"] += 1
                        found = True
                        break

                if not found:
                    if "color" in l.data['channels'][cha]:
                        color = l.data['channels'][cha]["color"]
                    else:
                        color = "grey"
                    channel_d = {"id": ch_id,
                                 "count": 1,
                                 "text_names": l.data['channels'][cha]["text_names"],
                                 "filename": os.path.basename(filename),
                                 "color": color}
                    channels_data_pre.append(channel_d)

        # reads sses stuff
        for i in (l.data['layout'].keys()):
            if i not in l.data['size']:
                continue
            # TODO - is this ok
            if (prime == 3) and (l.data['h_type'][i] == 3):
                continue
            counter = 0
            if len(layout_list) > 12000:
                start_name = "s" + str(i) + "_" + str(counter)
                end_name = "e" + str(i) + "_" + str(counter)
                name = str(i) + "_" + str(counter)
                start = utils_vector.angle_target(l.data['layout'][i],
                                                  l.data['angles'][i],
                                                  l.data['size'][i] * 0.1)
                end = utils_vector.angle_target(l.data['layout'][i],
                                                (l.data['angles'][i] + 180) % 360,
                                                l.data['size'][i] * 0.1)
                l_list[start_name] = start
                l_list[end_name] = end
                cl_list[name] = l.data['layout'][i]
            else:
                while counter < l.data['size'][i]:
                    start_name = "s" + str(i) + "_" + str(counter)
                    end_name = "e" + str(i) + "_" + str(counter)
                    name = str(i) + "_" + str(counter)
                    start = utils_vector.angle_target(l.data['layout'][i],
                                                      l.data['angles'][i],
                                                      l.data['size'][i] * 0.1)
                    end = utils_vector.angle_target(l.data['layout'][i],
                                                    (l.data['angles'][i] + 180) % 360,
                                                    l.data['size'][i] * 0.1)
                    l_list[start_name] = start
                    l_list[end_name] = end
                    cl_list[name] = l.data['layout'][i]
                    counter = counter + 1

        # add per sses stuff to layout per file lists
        if (l_list != {}):
            layout_all_lists.append(l_list)
            layout_call_lists.append(cl_list)
            layout_name_lists.append(filename)
            l_all_loaded.append(l)
        else:
            l_empty_loaded.append(l)

    ligands_data = sorted(ligands_data_pre, key=lambda x: x["count"], reverse=True)
    channels_data = sorted(channels_data_pre, key=lambda x: x["count"], reverse=True)

    # set ligand colors
    utils_color.set_ligands_colors(ligands_data, args.iligands, args.multi, args.dominik)
    # set channels colors
    utils_color.set_channels_colors(channels_data, args.multi)

    if layout_all_lists != []:
        # TODO: should be fixed is necessary to show only subset?
        img, template, ind, start_set, data = process_layout_list(
            l, layout_all_lists, layout_call_lists, layout_name_lists,
            start_set, l_all_loaded, otemplate, opacity, img, template,
            number_template, ind, data, ratio, img_data, args, ligands_data,
            img_data_ligands, channels_data)
    else:
        img, ind = draw_empty_layouts_sets(
            l, img, l_empty_loaded, args.output, l_all_loaded, args, ligands_data, opacity)

    # Save all wanted data
    # Save metadata for sestra web pages
    if args.data:
        fileor = open(args.data, "w")
        img_coords = {"upper_right": img.get_upper_right(),
                      "lower_left": img.get_lower_left()}
        json.dump({"metadata": img_coords, "sses": img_data, "ligands": img_data_ligands}, fileor, indent=4, separators=(',', ':'))
        fileor.close()

    # Output rotation file
    if args.orotation:
        fileor = open(args.orotation, "w")
        json.dump(data, fileor, indent=4, separators=(',', ':'))
        fileor.close()

    # Output template layout to a file
    if args.otemplate:
        output_template_layout(template, number_template, args.otemplate)

    # Output ligand colors file
    if args.oligands:
        output_ligands(ligands_data, args.oligands)

    if args.itemplate and (('only' not in vars(args)) or ("template" in args.only)):
        with open(args.itemplate, "r") as fd:
            data = json.load(fd)
            show_template_layout(img, data['layout'], ind, data['h_type'], data['size'],
                                 data['color'], data['angle'])

    # Save
    output = img.save()
    logging.info("INFO: Image saved into %s" % output)
    return output
