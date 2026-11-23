#!/usr/bin/env python3
# tool which generate for given protein its directory with
# * protein images
# and image per chain


import logging
import sys
import argparse
import time
import os
import subprocess
import shutil
import copy
import json

import utils_stats
import utils_common
import utils_annotation
import utils_transform_layout_sses_names
import utils_color_protein
import utils_convert_c_to_p
import utils_synchronize_rotation_files
import convert_protein_method_chain_to_std

sys.path.append('..')
import utils_generate_svg


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

def generate_whole_protein_from_scratch(protein, working_dir, annotation_file, sses_stats, start_l):
    # run three level whoole protein generation - first two levels are primar and secondary trimer part, the third is the rest
    layout_file = f"{working_dir}/layout-{protein}.json"

    command = f"./generate_covid_layout.py {annotation_file} ALL ALL {layout_file} {sses_stats} {start_l}"
    logging.debug(f"DEBUG: Generating layout for {protein}")
    logging.info(f"INFO: Command {command}")
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0) or (not os.path.isfile(layout_file)):
        log_fn = f"{working_dir}/generate_protein_layout_fs-fail.txt"
        utils_common.error_report(working_dir, command, None, err, log_fn)
        logging.error(f"generate_domain_layout_from_scratch.py failed - report save to {log_fn}")
        sys.exit()
        return False
    else:
        return True


def generate_cl(cif_file, protein, chain, working_dir, annotation_dir, layout_dir):
    # generate separate layouts for potein chain
    layout_file = f"{layout_dir}/layout-{protein}_{chain}00.json"
    annotation_file = f"{annotation_dir}/{protein}{chain}00-annotated.sses.json"
    command = f"./generate_domain_layout_from_scratch.py {annotation_file} {chain} {layout_file} {protein}"
    logging.debug(f"DEBUG: Generating layout for {protein} {chain} re - the rest of chain which is not in any domain computed yet")
    logging.info(f"INFO: Command {command}")
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0) or (not os.path.isfile(layout_file)):
        log_fn = f"{working_dir}/generate_protein_layout_fs-fail.txt"
        utils_common.error_report(working_dir, command, None, err, log_fn)
        logging.error(f"generate_domain_layout_from_scratch.py failed - report save to {log_fn}")
        sys.exit()
        return False
    else:
        return True


def generate_chain_layout(cif_file, chain_layout_file, chain_domain_layout_files, working_dir, annot_file, prot):
    command = f"./generate_chain_layout.py {annot_file} {chain_layout_file} {prot} {chain_domain_layout_files}"
    logging.debug(f"DEBUG: generating chain layout {chain_layout_file}")
    logging.debug(f"DEBUG: command {command}")
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0) or (not os.path.isfile(chain_layout_file)):
        log_fn = f"{working_dir}/generate_chain_layout_fs-fail.txt"
        utils_common.error_report(working_dir, command, None, err, log_fn)
        logging.error(f"generate_domain_layout.py failed - report save to {log_fn}")
        sys.exit()
        return False
    else:
        return True


def multiple_diagrams(detection_dir, layout_dir, first_svg):
    svg_args = argparse.Namespace(
        add_descriptions="0",
        data=None,
        debug=False,
#        irotation=f"{detection_dir}/rotation.json",
        irotation=f"",
        otemplate=f"{detection_dir}/ctemplate.json",
        layouts=[f"{layout_dir}"],
        opac='1',
        orotation=f"{layout_dir}/rotation.json",
        itemplate=None,
        output=f"{layout_dir}/imageCircle-multiple.svg",
        spike_protein=True,
        ligands=False
    )
    if first_svg is False:
        utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}")
    svg_args.output = f"{layout_dir}/imageCircle-multiple3.svg"
    svg_args.orotation = None
    svg_args.irotation = f"{layout_dir}/rotation.json"
    svg_args.itemplate = f"{detection_dir}/ctemplate.json"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}")
    svg_args.output = f"{layout_dir}/imageCircle-multiple2.svg"
    svg_args.opac = "0"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}")


def sync_diagrams(directory, rot_dir):
    logging.debug("Try to sync with the unified version")
    svg_args = argparse.Namespace(
        add_descriptions="0",
        data=None,
        debug=False,
        irotation=None,
        otemplate=f"{directory}/ctemplate.json",
        layouts=[f"{directory}"],
        opac='1',
        orotation=f"{directory}/rotation.json",
        itemplate=None,
        output=f"{directory}/imageCircle-multiple.svg",
        spike_protein=True,
        ligands=False
    )
    link = rot_dir
    os.mkdir(f"{directory}/old")
    try:
        shutil.copyfile(f"{rot_dir}/rotation.json", f"{directory}/old/rotation.json")
    except:
        logging.debug("no diagrams to synchronize with, create diagrams from scratch")
        # we have to sync the new versions of diagrams with the previous one
        utils_generate_svg.compute(svg_args)
        utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{directory}")
        multiple_diagrams(directory, directory, True)
    else:
        logging.debug(f"Old version of rotation file exists")
        # try to find common layout (in new version and old one)
        found = None
        for layout in os.listdir(directory):
            if layout.startswith("layout-") and layout.endswith(".json") and not layout.endswith("00.json"):
                # go through all layouts and find the one which is in old version
                try:
                    shutil.copyfile(f"{rot_dir}/{layout}", f"{directory}/old/{layout}")
                except:
                    continue
                else:
                    found = layout
                    break
        if found is None:
            logging.info(f"INFO: there is no common layout")
            command = f"curl from layout list {layout_list}"
            utils_common.error_report(directory, command, None, ret, f"{directory}/curl-warning-{time.time()}.txt")
            utils_generate_svg.compute(svg_args)
            utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{directory}")
        else:
            logging.info(f"common layout was found {found}")
            # sync can be done
            aux_layout_fn = f"layout-aaaaA00.json"
            aux_layout_full_fn = f"{directory}/{aux_layout_fn}"
            shutil.move(f"{directory}/{found}", aux_layout_full_fn)
            logging.info(f"move {directory}/{found} to {aux_layout_full_fn}")
            # synchronize names in old and new layout
            ret = utils_synchronize_rotation_files.sync_layout_n(directory, found, aux_layout_full_fn)
            if ret is False:
                sys.exit
            new_rotation = f"{directory}/rotation_new.json"
            svg_args.orotation = new_rotation
            # create rotation file for the new one
            utils_generate_svg.compute(svg_args)
            shutil.copyfile(new_rotation, f"{directory}/rotation_new222.json")
            # synchronize rotation files
            old_rotation = f"{directory}/old/rotation.json"
            aux_rotation = f"{directory}/rotation_aux.json"
            ret = utils_synchronize_rotation_files.synchronize_rotation_files(old_rotation, new_rotation, found, aux_rotation)
            # aux_rotation file should contain all layout in new_rotation with the names synced the old_rotation file
            if ret == -1:
                # probably the synchronize layout has no sses (TODO - not necessary - try to skip trivial layouts)
                command = f"utils_synchronize_rotation_files.synchronize_rotation_files({old_rotation}, {new_rotation}, {found}, {aux_rotation})"
                log_fn = f"{directory}/synchronize_rotation_files_warning.txt"
                err = "Rotation file does not contain record about domain" + str(found)
                utils_common.error_report(directory, command, None, err, log_fn)
                logging.error(f"synchronize_rotation_files failed - report save to {log_fn}")
                shutil.move(aux_layout_full_fn, f"{directory}/{found}")
                svg_args.orotation = f"{directory}/rotation.json"
                utils_generate_svg.compute(svg_args)
                utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{directory}")
            else:
                # cleanup_rotation_files
                svg_args.irotation = f"{directory}/rotation.json"
                utils_synchronize_rotation_files.cleanup_rotation_file(aux_rotation, found, aux_layout_fn, f"{directory}/rotation.json", f"{directory}/rotation2.json")
                logging.debug(f"synchronize_rotation_files was succesfull")
                os.remove(f"{directory}/{found}")
                shutil.move(aux_layout_full_fn, f"{directory}/{found}")
                svg_args.orotation = None
                utils_generate_svg.compute(svg_args)
    multiple_diagrams(directory, directory, True)



# read arguments
parser = argparse.ArgumentParser(
    description='generate single and multiple diagrams for a list of spike protein chains')
parser.add_argument('--p_list', required=True, action="store",
                    help="generate diagram database for given protein-chain entry")
parser.add_argument('--t1', required=True, action="store",
                    help="t1 dir")
parser.add_argument('--t2', required=True, action="store",
                    help="t2 dir")
parser.add_argument('--t3', required=True, action="store",
                    help="t3 dir")
parser.add_argument('--t4', action="store",
                    help="t4 dir - only for q9byf1")
parser.add_argument('--t5', action="store",
                    help="t5 dir - only for q9byf1")
parser.add_argument('--t6', action="store",
                    help="t6 dir - only for q9byf1")
parser.add_argument('--tz', action="store",
                    help="t zbytek dir")
parser.add_argument('--method', action="store",
                    help="non trimer part (tz) will be computed per chains -> 'chain', by default it will be cmputed as one bunch together -> 'None'")
parser.add_argument('--start', action="store",
                    help="starts layout (optional)")
parser.add_argument('--orientation', action="store",
                    help="starts layout (optional)")
parser.add_argument('--input_cif', required=True, action="store",
                    help="which will be put the result alphafold directory")
#parser.add_argument('--start_layout', action="store",
#                    help="define start layout")
#parser.add_argument('--consensus', required=True, action="store",
#                    help="donsensus file location")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')


args = parser.parse_args()


# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
aux = args.p_list.split("/")
p_list_name = aux[-1]
new_dir = f"generated-{p_list_name}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)

# start logger
logger = start_logger(args.debug, args.info, new_dir, p_list_name)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# we want to download read relevant cif files
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)
p_list = []
with open(args.p_list, "r") as fd:
    line = fd.readline()
    while line:
        s = line.split()
        p_list.append(s[0])
        line = fd.readline()

# we want to download read relevant cif files
cif_dir = f"{new_dir}/cif"
os.mkdir(cif_dir)
for p in p_list:
    cif_file = f"{cif_dir}/{p}.cif"
    if not (os.path.isfile(cif_file)):
        shutil.copyfile(f"{args.input_cif}/{p}.cif", f"{cif_dir}/{p}.cif")
        shutil.copyfile(f"{args.input_cif}/{p}.pdb", f"{cif_dir}/{p}.pdb")

# convert diagram and annotation files
t=[]
t.append(args.t1)
t.append(args.t2)
t.append(args.t3)
if args.t4 and args.t5 and args.t6:
    t.append(args.t4)
    t.append(args.t5)
    t.append(args.t6)
utils_convert_c_to_p.convert_diagram_c_to_p(t, len(t), args.p_list, annotation_dir)

# prepare annotation files
for p in p_list:
    utils_annotation.detect_list(annotation_dir, cif_dir, [p], new_dir, whole_protein=True)
    if args.tz and args.method != "chain":
        tz = args.tz
    else:
        tz = None
    utils_convert_c_to_p.convert_annotation_c_to_p(t, len(t), tz, p, annotation_dir)
    utils_color_protein.switch_color_to_chain(len(t),f"{annotation_dir}/{p}-annotation.sses.json")

diagram_file = f"{annotation_dir}/diagram.json"


logging.info("INFO: Starting save statistics and counting average protein")
stats_file = new_dir + "/sses-stats.json"
avg_protein = new_dir + "/avg-protein.json"

# use overprot statistics to find prime/sec sses
ret = utils_stats.copy_statistic(annotation_dir, diagram_file, stats_file)
if ret != []:
    command = f"utils_stats.save_statistics({annotation_dir}, {add_list}, 100, {stats_file}, {avg_protein})"
    log_fn = f"{new_dir}/save_statistics_fail.txt"
    utils_common.error_report(new_dir, command, None, None, log_fn)
    if ret == 0:
        logging.critical(f"save_statistics critically failed - report save to {log_fn}")
        sys.exit()
    else:
        logging.error(f"save_statistics failed - report save to {log_fn}")

# generate trimer layout
pstart=f"{new_dir}/start.json"
if not args.start:
    if len(p_list) < 3:
       n = 0
    else:
       n = 3
    generate_whole_protein_from_scratch(p_list[n], new_dir, f"{annotation_dir}/{p_list[n]}-annotation.sses.json", stats_file, "")
    start_layout = f"{new_dir}/layout-{p_list[n]}.json"
else:
    start_layout = args.start

shutil.copy(start_layout, pstart)
svg_args = argparse.Namespace(#
    add_descriptions="1",
    data=False,
    debug=False,
    irotation=None,
    otemplate=None,
    layouts=[f"{new_dir}/start.json"],
    opac='0',
    orotation=None,
    itemplate=None,
    output=f"{new_dir}/start.svg",
    ligands=False
    
)
utils_generate_svg.compute(svg_args)
for p in p_list:
    generate_whole_protein_from_scratch(p, new_dir, f"{annotation_dir}/{p}-annotation.sses.json", stats_file, pstart)

if (args.method == "chain") and (args.tz):
    for p in p_list:
         # generate "not trimer" layouts
        shutil.move(f"{new_dir}/layout-{p}.json", f"{new_dir}/layout-{p}_Z00.json")
        protein_layout_files = f"{new_dir}/layout-{p}_Z00.json"
#        protein_layout_files =""
        cif_file=f"{cif_dir}/{p}.cif"
        for tz_file in os.listdir(args.tz):
            if tz_file.startswith(p) and tz_file.endswith("annotated.sses.json"):
                print("i", p, tz_file[4], tz_file, args.tz)
                generate_cl(cif_file, p, tz_file[4], new_dir, args.tz, new_dir)
                protein_layout_files += f" {new_dir}/layout-{p}_{tz_file[4]}00.json"
        # aggregate trimer and non-trimer layouts to whole protein layout
        full_layout_file=f"{new_dir}/layout-{p}.json"
        print(f'generate_chain_layout({cif_file}, {full_layout_file}, {protein_layout_files}, {new_dir}, f"{annotation_dir}/{p}-annotated.sses.json", {p})')
        generate_chain_layout(cif_file, full_layout_file, protein_layout_files, new_dir, f"{annotation_dir}/{p}-annotation.sses.json", p)
        convert_protein_method_chain_to_std.rename_layout_sses(full_layout_file)


last_time = utils_common.time_log("2 layout ready ", start_time)

if args.orientation:
    sync_diagrams(new_dir, args.orientation)
else:
    multiple_diagrams(new_dir, new_dir, False)

# draw single image per chain and whole together
svg_args = argparse.Namespace(#
    add_descriptions="1",
    data=False,
    debug=False,
    irotation=f"{new_dir}/rotation.json",
    otemplate=None,
    layouts=[f""],
    opac='0',
    orotation=None,
    itemplate=None,
    output=f"",
    ligands=False,
)


for p in p_list:
    svg_args.layouts = [f"{new_dir}/layout-{p}.json"]
    svg_args.output = f"{new_dir}/image-{p}.svg"
    utils_generate_svg.compute(svg_args)

