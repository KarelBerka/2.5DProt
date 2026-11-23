#!/usr/bin/env python3
# generates for given layout directory family
# or for family which have the layouts in 2DProts dbs
# multi and single images


import logging
import sys
import argparse
import time
import os
import utils_common
import color_me
sys.path.append('..')
import utils_generate_svg

TWODPROTS = "https://2dprots.ncbr.muni.cz/"


def start_logger(debug, info, directory, protein):
    root = logging.getLogger()
    if debug:
        root.setLevel(logging.DEBUG)
    else:
        root.setLevel(logging.INFO)
    logfile = f"{directory}/main-{protein}.log"
    fh = logging.FileHandler(logfile)
    if debug:
        fh.setLevel(level=logging.DEBUG)
    else:
        fh.setLevel(level=logging.INFO)
    root.addHandler(fh)

    return root


def read_layout_list(family, directory):
    # download list of all domains
    domain_list_url = f"{TWODPROTS}/files/family/{family}/domain_list.txt"
    domain_list_filename = f"{directory}/domain_list.txt"
    utils_common.get_url_to_file(domain_list_url, domain_list_filename)
    with open(domain_list_filename, "r") as fd:
        # download each domain from the list
        line = fd.readline()
        while line:
            dom_split = line.split(",")
            layout_url = f"{TWODPROTS}/files/domain/{dom_split[0][2:-1]}{dom_split[1][2:-1]}{dom_split[2][2:-3]}/layout.json"
            layout_file = f"{directory}/layout-{dom_split[0][2:-1]}{dom_split[1][2:-1]}{dom_split[2][2:-3]}.json"
            utils_common.get_url_to_file(layout_url, layout_file)
            print(f"get layout-{dom_split[0][2:-1]}{dom_split[1][2:-1]}{dom_split[2][2:-3]}.json")
            line = fd.readline()


def draw_multiple_images(work_dir):
    # draw multiple images
    svg_args = argparse.Namespace(
        add_descriptions="multi,ligands",
        data=False,
        debug=False,
        irotation=None,
        iligands=None,
        itemplate=None,
        layouts=[f"{work_dir}"],
        ligands=True,
        oligands=f"{work_dir}/ligands.json",
        only="sses",
        opac='1',
        orotation=f"{work_dir}/rotation.json",
        otemplate=f"{work_dir}/ctemplate.json",
        output=f"{work_dir}/only_sses.svg",
    )
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output,
                                work_dir)

    svg_args.irotation = f"{work_dir}/rotation.json"
    svg_args.itemplate = f"{work_dir}/ctemplate.json"
    svg_args.iligands = f"{work_dir}/ligands.json"
    svg_args.orotation = None
    svg_args.output = f"{work_dir}/only_ligands.svg"
    svg_args.oligands = None
    svg_args.only = "ligands"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output,
                                work_dir)

    svg_args.output = f"{work_dir}/only_template.svg"
    svg_args.only = "template"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output,
                                work_dir)

    svg_args.output = f"{work_dir}/only_sses_opacity.svg"
    svg_args.opac = "0"
    svg_args.only = "sses"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output,
                                work_dir)

    svg_args.output = f"{work_dir}/only_channels.svg"
    svg_args.only = "channels"
    svg_args.channels = True
    svg_args.fill = "white"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output,
                                work_dir)

    svg_args.output = f"{work_dir}/only_all.svg"
    svg_args.only = "channels, sses, templates, ligands"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output,
                                work_dir)

    for fn in os.listdir(work_dir):
        if fn.startswith("layout-") and fn.endswith(".json"):
            svg_args.layouts = [f"{work_dir}/{fn}"]
            svg_args.output = f"{work_dir}/{fn[:-5]}.svg"
            svg_args.only = "channels, sses, ligands, template"
            svg_args.itemplate = None
            svg_args.add_descriptions = "1"


# read arguments
parser = argparse.ArgumentParser(
    description='generate family diagrams for a set of layout')
parser.add_argument('--family', action="store",
                    help="family name which layout should be drawn, (layouts are taken from 2DProts dbs")
parser.add_argument('-g', action="store",
                    help="layouts directory")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
new_dir = f"draw-{args.family}-{timestamp}"
os.makedirs(new_dir)
print("new diagrams directory", new_dir)

if args.family:
    if args.g:
        layouts_dir = args.g
    else:
        layouts_dir=f"{new_dir}/layout"
        os.makedirs(layouts_dir)
        l_list = read_layout_list(args.family, layouts_dir)

# start logger
logger = start_logger(args.debug, args.info, new_dir, "aux")

# create stats.json file
link_variance_data = f"https://overprot.ncbr.muni.cz/data/db/family/consensus_sses//consensus-{args.family}.sses.json"
consensus_file = f"{new_dir}/consensus-{args.family}.sses.json"
ret = utils_common.wget_file(link_variance_data, new_dir, consensus_file, option="--no-check-certificate")
if not ret:
    sys.exit()

layouts2_dir= f"{new_dir}/layout2"
os.makedirs(layouts2_dir)
for file_name in os.listdir(layouts_dir):
    if file_name.startswith("layout-"):
        color_me.switch_color_to_variance(f"{layouts_dir}/{file_name}", f"{layouts2_dir}/{file_name}",  consensus_file)

print("draw_multiple_images", f"{new_dir}")
for fn in os.listdir(f"{layouts_dir}/"):
    if fn.startswith("layout-"):
        os.symlink(f"../{layouts_dir}/{fn}", f"{new_dir}/{fn}")

draw_multiple_images(f"{new_dir}")

sys.exit()

for dn in os.listdir(layouts_dir):
    if dn.startswith(args.family) and not dn.endswith(".json") and not dn.endswith(".svg") and not dn.endswith(".png"):
        print("dirs", f"{layouts_dir}/{dn}")
        os.makedirs(f"{new_dir}/{dn}")
        for fn in os.listdir(f"{layouts_dir}/{dn}"):
            if fn.startswith("layout-"):
                os.symlink(f"../../{layouts_dir}/{dn}/{fn}", f"{new_dir}/{dn}/{fn}")
        draw_multiple_images(f"{new_dir}/{dn}")

draw_multiple_images(f"{new_dir}")
sys.exit()

# draw single images -
svg_args = argparse.Namespace(
    add_descriptions="multi,ligands",
    data=False,
    debug=False,
    irotation=None,
    itemplate=None,
    iligands=None,
    layouts=[f"{new_dir}"],
    opac='1',
    orotation=f"{new_dir}/rotation.json",
    otemplate=f"{new_dir}/ctemplate.json",
    oligands=f"{new_dir}/ligands.json",
    output=f"{new_dir}/only_sses.svg",
    ligands=True,
    only="sses")

for fn in os.listdir(old_dir):
    if fn.startswith("layout-") and fn.endswith(".json"):
        svg_args.layouts = [f"{old_dir}/{fn}"]
        svg_args.output = f"{new_dir}/{fn[:-5]}.svg"
        svg_args.itemplate = None
        svg_args.add_descriptions = "1"
        utils_generate_svg.compute(svg_args)

print(f"output directory: {new_dir}")
