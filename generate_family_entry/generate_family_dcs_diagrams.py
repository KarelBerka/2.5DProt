#!/usr/bin/env python3
# tool which generate for given family name its directory with
# * single images
# * multi images per whole family
# * multi images for its s35 subfamilies

import logging
import logging.handlers
import sys
import argparse
import time
import os
import shutil
import subprocess
import json
import numpy

import utils_common
import utils_stats
# import utils_find_univ_start
import utils_synchronize_rotation_files
import utils_generate_svg
import utils_sort_to_prefamilies
import utils_pdbe_api

sys.path.append('..')
import utils_links


def start_logger(debug, info, directory, family):
    root = logging.getLogger()
    if info:
        root.setLevel(logging.INFO)
    if debug:
        root.setLevel(logging.DEBUG)

    logfile = f"{directory}/main-{family}.log"
    fh = logging.FileHandler(logfile)
    if info:
        fh.setLevel(level=logging.INFO)
    if debug:
        fh.setLevel(level=logging.DEBUG)
    root.addHandler(fh)
    return root


def download_annotated_files(annotation_dir, family, diagram_file):
    #  download annotation files and diagram file for given family
    #  overprot_link = "https://overprot.ncbr.muni.cz/data/db-dev/"
    overprot_link = "https://overprot.ncbr.muni.cz/data/db/"

    link_annot_zip = f"{overprot_link}/family/annotation/annotation-{family}.zip"
    annot_zip_file = f"{annotation_dir}/annotation-{family}.zip"
    ret = utils_common.curl_file(link_annot_zip, new_dir, annot_zip_file, option="--insecure")
    if not ret:
        sys.exit()

    ret = utils_common.unzip_archive(annot_zip_file, annotation_dir)
    if ret != 0:
        sys.exit()

    link_diagram = f"{overprot_link}/family/diagram/diagram-{family}.json"
    ret = utils_common.curl_file(link_diagram, new_dir, diagram_file, option="--insecure")
    if not ret:
        sys.exit()


def copy_ligand_json_to_dir(args, annotation_dir, new_dir):
    # we want to download read relevant ligand cif file
    ligand_filename = f"{annotation_dir}/ligands_data-{args.protein}.json"
    if (args.ligands_json_dir):
        logging.info("INFO: Copy ligand json data")
        shutil.copyfile(f"{args.ligands_json_dir}/{args.protein[1:3]}/{args.protein}.json", ligand_filename)
    else:
        logging.info("INFO: Download ligand json data")
        ligand_url = utils_links.link_ligands(args.protein)
        if not utils_common.curl_file(ligand_url, new_dir, ligand_filename, option="--insecure"):
            utils_common.error_report(new_dir, "copy ligand json problem - see curl fail", None, "", f"{new_dir}/missing_ligand_file.txt")
            logging.error("generate_domain_layout_from_scratch.py failed - report save to missing_ligand_file.txt")
            sys.exit()
    return ligand_filename


def copy_channel_json_to_dir(protein, annotation_dir, new_dir):
    # we want to download read relevant ligand cif file
    channel_filename = f"{annotation_dir}/channel-{protein}.json"
    logging.info("INFO: Download channel json data")
    channel_url = utils_links.link_channels(protein)
    if not utils_common.curl_file(channel_url, new_dir, channel_filename, option="--insecure", err_report=False):
        shutil.copyfile("empty_channel.json", channel_filename)
    return channel_filename


def rotate_ligand(ligand_filename, ligand_rotation, aux_ligand_filename):
    with open(ligand_filename, "r") as lf:
        l_data = json.load(lf)

    with open(ligand_rotation, "r") as rf:
        r_data = json.load(rf)

    for layout in l_data:
        l_data[layout]["center"] = list(numpy.array(numpy.dot(r_data["rotation"], numpy.array(l_data[layout]["center"])) + numpy.array([r_data["translation"][0][0], r_data["translation"][1][0], r_data["translation"][2][0]])))

    with open(aux_ligand_filename, "w") as lf:
        json.dump(l_data, lf, indent=4, separators=(',', ':'))


def rotate_channel(channel_filename, channel_rotation, aux_channel_filename):
    with open(channel_filename, "r") as lf:
        l_data = json.load(lf)

    with open(channel_rotation, "r") as rf:
        r_data = json.load(rf)

    aux = list(l_data.keys())
    if "Error" not in aux:
        if "Channels" in l_data:
          for channel_types in l_data["Channels"]:
            for ch in l_data["Channels"][channel_types]:
                for ch_inst in ch["Profile"]:
                    orig = [ch_inst['X'], ch_inst['Y'], ch_inst['Z']]
                    new = list(numpy.array(numpy.dot(r_data["rotation"], numpy.array(orig)) + numpy.array([r_data["translation"][0][0], r_data["translation"][1][0], r_data["translation"][2][0]])))
                    ch_inst['X'] = new[0]
                    ch_inst['Y'] = new[1]
                    ch_inst['Z'] = new[2]

    with open(aux_channel_filename, "w") as lf:
        json.dump(l_data, lf, indent=4, separators=(',', ':'))


def generate_layout_list(annotation_dir, pcd_list, color_layout_filename,
                         start_layout_filename, ostats, detection_dir, not_ligands):
    ret_layouts = []
    counter = 0
    for protein, chain, domain in pcd_list:
        pcd_name = f"{protein}{chain}{domain}"
        logging.debug(f"Generate layout for {pcd_name}")
        annotation_file = f"{annotation_dir}/{pcd_name}-annotated.sses.json"
        layout_filename = f"{detection_dir}/layout-{pcd_name}.json"
        stats_json = f"{detection_dir}/sses-stats.json"
        if not_ligands:
            command = f"./generate-layout.py {annotation_file} {stats_json} {chain} {domain} {layout_filename} {start_layout_filename} {ostats}"
        else:
            link_ligands = utils_links.link_ligands(protein)
            ligand_filename = f"{annotation_dir}/ligands_data-{protein}.json"
            ret = utils_common.curl_file(link_ligands, new_dir, ligand_filename, option="--insecure")
            if not ret:
                sys.exit()
            link_ligands_rotation = utils_links.link_ligands_rotation(protein, chain, domain)
            ligand_rotation = f"{annotation_dir}/ligands_data-{protein}{chain}{domain}-rotation.json"
            ret = utils_common.curl_file(link_ligands_rotation, new_dir, ligand_rotation, option="--insecure")
            if not ret:
                sys.exit()
            aux_ligand_filename = f"{annotation_dir}/{protein}{chain}{domain}-ligands.json"
            # change ligand rotation to sesctrannotator one
            rotate_ligand(ligand_filename, ligand_rotation, aux_ligand_filename)
            # rotate channel files
            channel_filename = f"{annotation_dir}/channel-{protein}.json"
            aux_channel_filename = f"{annotation_dir}/{protein}{chain}{domain}-channel.json"
            rotate_channel(channel_filename, ligand_rotation, aux_channel_filename)
            command = f"./generate-layout-ligands.py --annotation_file {annotation_file} --stats_file {stats_json} --domain {domain} --chain {chain} --output {layout_filename} --start_layout {start_layout_filename} --ostats {ostats} --ligands {aux_ligand_filename} --channels {aux_channel_filename}"

        logging.info(f"INFO: generating layout for {pcd_name}, start layout: {start_layout_filename != None}")
        counter += 1
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()
        if (r != 0) or (not os.path.isfile(layout_filename)):
            log_fn = f"{detection_dir}/generate_layout-fail-{time.time()}.txt"
            utils_common.error_report(detection_dir, command, None, err, log_fn)
            logging.error(f"generate_layout failed - report save to {log_fn}")
        else:
            ret_layouts.append((protein, chain, domain))
    return ret_layouts


def multiple_diagrams(detection_dir, layout_dir, first_svg):
    svg_args = argparse.Namespace(
        add_descriptions="",
        multi=True,
        layouts=[f"{layout_dir}"],
        ligands=False,
        transparency=False,
        debug=False,
        irotation=f"{detection_dir}/rotation.json",
        itemplate=None,
        iligands=f"{detection_dir}/ligands.json",
        orotation=None,
        otemplate=f"{detection_dir}/{args.family}.json",
        output=f"{layout_dir}/imageCircle-multiple.svg",
        oligands=None,
        data=None,
        groups=False,
    )
    if first_svg is False:
        utils_generate_svg.compute(svg_args)
    svg_args.otemplate = None
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")
    svg_args.output = f"{layout_dir}/imageCircle-multiple_3.svg"
    svg_args.itemplate = f"{detection_dir}/{args.family}.json"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")
    svg_args.output = f"{layout_dir}/imageCircle-multiple_2.svg"
    svg_args.transparency = True
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")
    svg_args.ligands = True
    svg_args.iligands = f"{detection_dir}/ligands.json"
    svg_args.add_descriptions = "ligands"
    svg_args.output = f"{layout_dir}/imageCircle-multiple_4.svg"
    if detection_dir != layout_dir:
        svg_args.oligands = f"{layout_dir}/ligands.json"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")
    svg_args.ligands = True
    svg_args.iligands = False
    svg_args.add_descriptions = "ligands"
    svg_args.output = f"{layout_dir}/{args.family}_5.svg"
    svg_args.dominik = "ch"
    if detection_dir != layout_dir:
        svg_args.oligands = f"{layout_dir}/ligands.json"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")
    svg_args.ligands = True
    svg_args.iligands = False
    svg_args.add_descriptions = "ligands"
    svg_args.output = f"{layout_dir}/{args.family}_6.svg"
    svg_args.dominik = "f"
    if detection_dir != layout_dir:
        svg_args.oligands = f"{layout_dir}/ligands.json"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")


def new_multiple_diagrams(detection_dir, layout_dir, first_svg, multiple_name, ligands_done=False):
    svg_args = argparse.Namespace(
        add_descriptions="ligands",
        multi=True,
        debug=False,
        irotation=f"{detection_dir}/rotation.json",
        itemplate=None,
        iligands=None,
        layouts=[f"{layout_dir}"],
        ligands=True,
        only="ligands",
        transparency=False,
        orotation=None,
        otemplate=f"{layout_dir}/template.json",
        output=f"{layout_dir}/{multiple_name}_ligands.svg",
        oligands=f"{layout_dir}/ligands.json",
        data=None,
        fill=False,
        groups=False,
    )
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")

    svg_args.add_descriptions = ""
    svg_args.only = "sses"
    svg_args.output = f"{layout_dir}/{multiple_name}_sses.svg"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")

    svg_args.only = "sses"
    svg_args.transparency = True
    svg_args.output = f"{layout_dir}/{multiple_name}_sses2.svg"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")

    svg_args.only = "template"
    svg_args.itemplate = f"{layout_dir}/template.json"
    svg_args.output = f"{layout_dir}/{multiple_name}_template.svg"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")

    svg_args.only = "channels"
    svg_args.output = f"{layout_dir}/{multiple_name}_channels.svg"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")

    svg_args.add_descriptions = "ligands"
    svg_args.iligands = f"{layout_dir}/ligands.json"
    svg_args.only = "ligands, sses, template, channels"
    svg_args.output = f"{layout_dir}/{multiple_name}_all.svg"
    svg_args.fill = True
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")


def rename_layout_fn(fn, label_old, label_new):
    with open(fn, "r") as fd:
        data_r = json.load(fd)

    data_w = {}
    for layout in data_r.keys():
        if layout == label_old:
            data_w[label_new] = data_r[label_old]
        else:
            data_w[layout] = data_r[layout]
    with open(f"{fn}", "w") as fd:
        json.dump(data_w, fd)


def sync_diagrams(family, directory, layout_list, not_ligands):
    logging.debug("Try to sync with the old version")
    svg_args = argparse.Namespace(
        add_descriptions="multi",
        layouts=[f"{directory}"],
        multi=True,
        ligands=False,
        transparency=False,
        debug=False,
        irotation=None,
        itemplate=None,
        iligands=None,
        orotation=f"{directory}/rotation.json",
        otemplate=f"{directory}/{args.family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{directory}/{args.family}.svg",
        data=None,
        groups=False,
    )
    if not_ligands:
        svg_args.ligands = False
    link = f"https://2dprots.ncbr.muni.cz/files/family/{family}/rotation.json"
    os.mkdir(f"{directory}/2dprots")
    ret = utils_common.curl_file(link, directory, f"{directory}/2dprots/rotation.json", err_report=False)
    logging.debug(f"Old version of rotation file exists: {ret}")
    if ret is False:
        logging.debug("no diagrams to synchronize with, create diagrams from scratch")
        # we have to sync the new versions of diagrams with the previous one
        utils_generate_svg.compute(svg_args)
        utils_common.convert_svg_to_png(svg_args.output, f"{directory}")
    else:
        logging.debug("Try to synchronize with old version")
        # try to find common layout (in new version and old one)
        found = None
        for layout in layout_list:
            # go through all layouts and find the one which is in old version
            pcd = f"{layout[0]}{layout[1]}{layout[2]}"
            fn = f"layout-{pcd}.json"
            link = f"https://2dprots.ncbr.muni.cz/files/domain/{pcd}/layout.json"
            ret = utils_common.curl_file(link, directory, f"{directory}/2dprots/{fn}", err_report=False)
            if ret:
                found = fn
                break
        if found is None:
            logging.info("INFO: there is no common layout")
            command = f"curl from layout list {layout_list} (trying to sync layout rotation with old version)"
            utils_common.error_report(directory, command, None, ret, f"{directory}/curl-warning-{time.time()}.txt")
            utils_generate_svg.compute(svg_args)
            utils_common.convert_svg_to_png(svg_args.output, f"{directory}")
        else:
            logging.info(f"common layout was found {found}")
            # sync can be done
            aux_layout_fn = "layout-aaaaA00.json"
            aux_layout_full_fn = f"{directory}/{aux_layout_fn}"
            shutil.move(f"{directory}/{fn}", aux_layout_full_fn)
            logging.info(f"move {directory}/{fn} to {aux_layout_full_fn}")
            # synchronize names in old and new layout
            ret = utils_synchronize_rotation_files.sync_layout_n(directory, found, aux_layout_full_fn)
            if ret is False:
                sys.exit()
            new_rotation = f"{directory}/rotation_new.json"
            svg_args.orotation = new_rotation
            # create rotation file for the new one
            utils_generate_svg.compute(svg_args)
            # synchronize rotation files
            old_rotation = f"{directory}/2dprots/rotation.json"
            aux_rotation = f"{directory}/rotation_aux.json"
            ret = utils_synchronize_rotation_files.synchronize_rotation_files(old_rotation, new_rotation, fn, aux_rotation)
            # aux_rotation file should contain all layout in new_rotation with the names synced the old_rotation file
            if ret == -1:
                # probably the synchronize layout has no sses (TODO2: - not necessary - try to skip trivial layouts)
                command = f"utils_synchronize_rotation_files.synchronize_rotation_files({old_rotation}, {new_rotation}, {fn}, {aux_rotation})"
                log_fn = f"{directory}/synchronize_rotation_files_warning.txt"
                err = "Rotation file does not contain record about domain" + str(fn)
                utils_common.error_report(directory, command, None, err, log_fn)
                logging.error(f"synchronize_rotation_files failed - report save to {log_fn}")
                shutil.move(aux_layout_full_fn, f"{directory}/{fn}")
                svg_args.orotation = f"{directory}/rotation.json"
                utils_generate_svg.compute(svg_args)
                utils_common.convert_svg_to_png(svg_args.output, f"{directory}")
            else:
                # cleanup_rotation_files
                svg_args.irotation = f"{directory}/rotation.json"
                utils_synchronize_rotation_files.cleanup_rotation_file(aux_rotation, fn, aux_layout_fn, f"{directory}/rotation.json", f"{directory}/rotation2.json")
                logging.debug("synchronize_rotation_files was succesfull")
                os.remove(f"{directory}/{fn}")
                shutil.move(aux_layout_full_fn, f"{directory}/{fn}")
                svg_args.orotation = None
                utils_generate_svg.compute(svg_args)


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagrams for a family')
parser.add_argument('--family', required=True, action="store",
                    help="generate diagram database for given family")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')
parser.add_argument('-n', '--not-ligands', action='store_true',
                    help='result diagrams does not contain ligands')
parser.add_argument('--new', action='store_true',
                    help='generate gnuplot2 output diagrams')

args = parser.parse_args()


# create new family diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
new_dir = "generated-" + args.family + "-" + timestamp
os.mkdir(new_dir)

# start logger
logger = start_logger(args.debug, args.info, new_dir, args.family)
g1, g2 = utils_common.git_version()
logging.error(f"2DProts version: {g1} {g2}")
last_time = utils_common.time_log("0 (Befor process starts)", last_time)

# we want to find family domains
logging.info("INFO: find family pcd set")
add_list = utils_pdbe_api.get_family_domains(args.family)
#add_list= {("4zwb", "A", "00"), ("4zw9", "A", "00"), ("4zwc", "A", "00"), ("4zwc", "B", "00"), ("5c65", "A", "00"), ("5c65", "B", "00"), ("7crz", "A", "00"), ("7sps", "A", "00"), ("7sps", "B", "00"), ("7spt", "A", "00")}
#print("add_list", add_list)

if len(add_list) == 0:
    with open(f"{new_dir}/domain_list.txt", "w") as fdl:
        json.dump([], fdl, indent=4)

    # can't find domain list - can't continue
    sys.exit()

# test whether we should not finished
if add_list == set():
    logging.info("INFO: No domains to add")
    sys.exit()
else:
    logging.info(f"INFO: Add {len(add_list)} domains")

# download annotated file and family statistic file diagram-{args.family}.json
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)
diagram_file = f"{annotation_dir}/diagram-{args.family}.json"
download_annotated_files(annotation_dir, args.family, diagram_file)
# get relevant ligand json and cif file
# ligand_filename = copy_ligand_json_to_dir(args, annotation_dir, new_dir)

for dom in add_list:
    channel_filename = copy_channel_json_to_dir(dom[0], annotation_dir, new_dir)
last_time = utils_common.time_log(f"1 ({len(add_list)} chain annotation downloaded)", last_time)

# generate sses_matches svg diagram - not necessary now
# ./show-sses-matches-for-family.py new_dir new_dir/sses-matches.svg

# copy statistics form diagram-{family}. json file to sses-stats.json statistics,
logging.info("INFO: Starting save statistics and counting average protein")
stats_file = new_dir + "/sses-stats.json"
ret = utils_stats.copy_statistic(annotation_dir, diagram_file, stats_file)
if not ret:
    command = f"utils_stats.copy_statistics({annotation_dir}, {diagram_file}, {stats_file})"
    log_fn = f"{new_dir}/copy_statistics_fail.txt"
    utils_common.error_report(new_dir, command, None, ret, log_fn)
    logging.error(f"save_statistics failed - report save to {log_fn}")
    sys.exit()

# find universal start layout
logging.info("INFO: Find universal start layout")
for pcd in add_list:
    universal_start = pcd          # utils_find_univ_start.find_universal_start(add_list, annotation_dir, 7)
    break
# TODO0: rewrite
logging.info(f"INFO: Universal_start {universal_start}")
last_time = utils_common.time_log("3 (stats, average protein, universal start layout counted)", last_time)
# generate layout for universal start domain
logging.debug("generate universal start {universal_start} layout")
res = generate_layout_list(annotation_dir, [universal_start], None, None, None, new_dir, args.not_ligands)
if res == []:
    logging.critical("can't create start layout")
    sys.exit()
start_layout_filename = f"{new_dir}/start-layout.json"
shutil.move(f"{new_dir}/layout-{universal_start[0]}{universal_start[1]}{universal_start[2]}.json", start_layout_filename)
# generate layouts for all annotated_list

logging.debug("generate all layouts")
color_layout_filename = f"{annotation_dir}/consensus-template.sses.json"
layout_list = generate_layout_list(annotation_dir, add_list, color_layout_filename, start_layout_filename, None, new_dir, args.not_ligands)
last_time = utils_common.time_log("4 (all layouts computed)", last_time)

link_variance_data = f"https://overprot.ncbr.muni.cz/data/db/family/consensus_sses//consensus-{args.family}.sses.json"
consensus_file = f"{new_dir}/consensus-{args.family}.sses.json"
ret = utils_common.wget_file(link_variance_data, new_dir, consensus_file)
if not ret:
    sys.exit()

# for file_name in os.listdir(new_dir):
#    if file_name.startswith("layout-"):
#        color_me.switch_color_to_variance(f"{new_dir}/{file_name}", f"{new_dir}/{file_name}", consensus_file)
####

# generate images ------------------------------------------
# if it is not update wget rotation.json and files from previous version of generating this family (to preserve the rotation)
# generate multiple images
sync_diagrams(args.family, new_dir, layout_list, args.not_ligands)
svg_list = []
if not args.new:
    multiple_diagrams(new_dir, new_dir, True)
    # generate per s35 cluster diagrams
    utils_sort_to_prefamilies.sort_to_s35clusters(new_dir, args.family)
    for dirname in os.listdir(new_dir):
        if dirname.startswith(args.family) and (not dirname.endswith(".svg")) and (not dirname.endswith(".png")) and (not dirname.endswith(".json")):
            logging.info(f"INFO: draw cluster {dirname} multi diagrams")
            # find weather the s35 cluster is not equivalent to the whole family
            s35_layout_list = []
            s35_dir = new_dir + "/" + dirname
            for filename in os.listdir(s35_dir):
                s35_layout_list.append(filename)
            if len(s35_layout_list) == len(layout_list):
                logging.info("INFO: s35 cluster is the same as whole family, just copy them")
                # the s35 group is the same as the whole group
                os.symlink(f"../{args.family}.svg", f"{s35_dir}/{args.family}.svg")
                os.symlink(f"../{args.family}.png", f"{s35_dir}/{args.family}.png")
                os.symlink(f"../{args.family}_2.svg", f"{s35_dir}/{args.family}_2.svg")
                os.symlink(f"../{args.family}_2.png", f"{s35_dir}/{args.family}_2.png")
                os.symlink(f"../{args.family}_3.svg", f"{s35_dir}/{args.family}_3.svg")
                os.symlink(f"../{args.family}_3.png", f"{s35_dir}/{args.family}_3.png")
                os.symlink(f"../{args.family}_4.svg", f"{s35_dir}/{args.family}_4.svg")
                os.symlink(f"../{args.family}_4.png", f"{s35_dir}/{args.family}_4.png")
                os.symlink("../ligands.json", f"{s35_dir}/ligands.json")
                os.symlink("../domain_list.txt", f"{s35_dir}/domain_list.txt")
            elif len(os.listdir(f"{new_dir}/{dirname}")) != 0:    # the directory is notempty
                multiple_diagrams(new_dir, f"{new_dir}/{dirname}", False)
            else:
                os.rmdir(f"{new_dir}/{dirname}")
            if len(s35_layout_list) != 0:                         # if list is empty, the directory was not created at all
                with open(f"{s35_dir}/domain_list.txt", "w") as fd:
                    json.dump(svg_list, fd, indent=4)

    # generate single images
    svg_args = argparse.Namespace(
        add_descriptions="chain",
        layouts=[f"{new_dir}"],
        ligands=True,
        transparency=True,
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        itemplate=None,
        iligands=f"{new_dir}/ligands.json",
        orotation=None,
        otemplate=None,
        oligands=None,
        output="",
        data=False,
        groups=False,
    )

    c = 0
    for i in layout_list:
        prot = i[0]
        chain = i[1]
        domain = i[2]
        pcd = f"{i[0]}{i[1]}{i[2]}"
        svg_args.layouts = [f"{new_dir}/layout-{pcd}.json"]
        svg_args.data = f"{new_dir}/image-{prot}_{chain}{domain}.json"
        svg_args.ligands = False
        svg_args.iligands = None
        svg_args.output = f"{new_dir}/image-{prot}_{chain}{domain}.svg"
        svg_args.channels = False
        r = utils_generate_svg.compute(svg_args)
        if r is not None:
            c += 1
        svg_list.append(list(i))
        svg_args.output = f"{new_dir}/image-{prot}_{chain}{domain}l.svg"
        svg_args.ligands = True
        svg_args.iligands = f"{new_dir}/ligands.json"
        r = utils_generate_svg.compute(svg_args)

else:
    new_multiple_diagrams(new_dir, new_dir, True, args.family)
    utils_sort_to_prefamilies.sort_to_s35clusters(new_dir, args.family)
    for dirname in os.listdir(new_dir):
        if dirname.startswith(args.family) and (not dirname.endswith(".svg")) and (not dirname.endswith(".png")) and (not dirname.endswith(".json")):
            logging.info(f"INFO: draw cluster {dirname} multi diagrams")
            # find weather the s35 cluster is not equivalent to the whole family
            s35_layout_list = []
            s35_dir = new_dir + "/" + dirname
            cluster_name = dirname
            for filename in os.listdir(s35_dir):
                protein=filename[7:11]
                chain=filename[11:-7]
                domain=filename[-7:-5]
                s35_layout_list.append([protein, chain, domain])
            if len(s35_layout_list) == len(layout_list):
                logging.info("INFO: s35 cluster is the same as whole family, just copy them")
                # the s35 group is the same as the whole group
                os.symlink(f"../{args.family}_all.svg", f"{s35_dir}/{cluster_name}_all.svg")
                os.symlink(f"../{args.family}_all.png", f"{s35_dir}/{cluster_name}_all.png")
                os.symlink(f"../{args.family}_sses.svg", f"{s35_dir}/{cluster_name}_sses.svg")
                os.symlink(f"../{args.family}_sses.png", f"{s35_dir}/{cluster_name}_sses.png")
                os.symlink(f"../{args.family}_sses2.svg", f"{s35_dir}/{cluster_name}_sses2.svg")
                os.symlink(f"../{args.family}_sses2.png", f"{s35_dir}/{cluster_name}_sses2.png")
                os.symlink(f"../{args.family}_ligands.svg", f"{s35_dir}/{cluster_name}_ligands.svg")
                os.symlink(f"../{args.family}_ligands.png", f"{s35_dir}/{cluster_name}_ligands.png")
                os.symlink(f"../{args.family}_channels.svg", f"{s35_dir}/{cluster_name}_channels.svg")
                os.symlink(f"../{args.family}_channels.png", f"{s35_dir}/{cluster_name}_channels.png")
                os.symlink(f"../{args.family}_template.svg", f"{s35_dir}/{cluster_name}_template.svg")
                os.symlink(f"../{args.family}_template.png", f"{s35_dir}/{cluster_name}_template.png")
                os.symlink("../ligands.json", f"{s35_dir}/ligands.json")
                os.symlink("../domain_list.txt", f"{s35_dir}/domain_list.txt")
            elif len(os.listdir(f"{new_dir}/{dirname}")) != 0:    # the directory is notempty
                new_multiple_diagrams(new_dir, f"{new_dir}/{dirname}", False, cluster_name, f"{new_dir}/ligands.json")
            else:
                os.rmdir(f"{new_dir}/{dirname}")
            if len(s35_layout_list) != 0:   # if the list is empty, the directory was not created at all
                with open(f"{s35_dir}/domain_list.txt", "w") as fd:
                    for domain in s35_layout_list:
                        fd.write(f"('{domain[7:11]}', '{domain[11:-7]}', '{domain[-7:-5]}')\n")
#                    json.dump(s35_layout_list, fd, indent=4)

    # generate single images
    svg_args = argparse.Namespace(
        add_descriptions="chain",
        layouts=[f"{new_dir}"],
        ligands=True,
        transparency=True,
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        itemplate=None,
        iligands=f"{new_dir}/ligands.json",
        orotation=None,
        otemplate=None,
        oligands=None,
        output="",
        data=False,
        fill=False,
        groups=False,
    )

    c = 0
    svg_list = []
    for i in layout_list:
        prot = i[0]
        chain = i[1]
        domain = i[2]
        pcd = f"{i[0]}{i[1]}{i[2]}"
        svg_args.layouts = [f"{new_dir}/layout-{pcd}.json"]
        svg_args.data = None
        svg_args.multi = False
        svg_args.iligands = None
        svg_args.output = f"{new_dir}/{prot}_{chain}{domain}.svg"
        svg_args.data = f"{new_dir}/{prot}_{chain}{domain}.json"
        svg_args.channels = True
        svg_args.ligands = True
        svg_args.only = "ligands, sses, template, channels"
        svg_args.fill = True
        svg_args.groups = True
        r = utils_generate_svg.compute(svg_args)
        svg_list.append(list(i))
        if r is not None:
            c += 1

# final cleanup
now_time = time.time()
if (not args.debug) and (len(add_list) == c):
    shutil.rmtree(f"{annotation_dir}")
    shutil.rmtree(f"{new_dir}/2dprots")
    if os.path.isfile(f"{new_dir}/rotation2.json"):
        os.remove(f"{new_dir}/rotation2.json")
    if os.path.isfile(f"{new_dir}/rotation_new.json"):
        os.remove(f"{new_dir}/rotation_new.json")
    if os.path.isfile(f"{new_dir}/rotation_aux.json"):
        os.remove(f"{new_dir}/rotation_aux.json")
    if os.path.isfile(f"{new_dir}/sses-stats.json"):
        os.remove(f"{new_dir}/sses-stats.json")

    last_time = utils_common.time_log("5 OK all diagrams created", last_time)
    last_time = utils_common.time_log("6 OK sum", start_time)
else:
    last_time = utils_common.time_log("5 BUG all diagrams created", last_time)
    last_time = utils_common.time_log("6 BUG sum", start_time)

if (len(svg_list) == c):
    with open(f"{new_dir}/domain_list.txt", "w") as fdl:
        for i in svg_list:
            fdl.write(f"{i}\n")

print("generating succesfully complete\n")
