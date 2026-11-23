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
import utils_common
import utils_annotation
import utils_aggregate_annotations
import utils_pf_stats
import utils_generate_svg
import utils_pdbe_api
import utils_sses_types


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


def generate_layout_list(annotation_dir, pcd_list, color_layout_filename,
                         start_layout_filename, ostats, detection_dir, not_ligands):
    ret_layouts = []
    for protein, domain, chain in pcd_list:
        pcd_name = f"{protein}{domain}{chain}"
        logging.debug(f"Generate layout for {pcd_name}")
        annotation_file = f"{annotation_dir}/{protein}-annotated.sses.json"
        layout_filename = f"{detection_dir}/layout-{pcd_name}.json"
        stats_json = f"{detection_dir}/sses-stats.json"
        command = f"../generate-layout-ligands.py {annotation_file} {stats_json} {domain} {chain} {layout_filename} {start_layout_filename} {ostats}"
        logging.info(f"INFO: generating layout for {pcd_name}, start layout: {start_layout_filename != None}")
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()
        if (r != 0) or (not os.path.isfile(layout_filename)):
            log_fn = f"{detection_dir}/generate_layout-fail-{time.time()}.txt"
            utils_common.error_report(detection_dir, command, None, err, log_fn)
            logging.error(f"generate_layout failed - report save to {log_fn}")
        else:
            ret_layouts.append((protein, domain, chain))
    return ret_layouts


def generate_chain_layout(cif_file, chain_layout_file, chain_domain_layout_files, working_dir, annot_file, prot, not_ligands, annot_dir):
    # aggregate several layouts to one multi-layout - either domains to chain or chain to protein ones
    if not_ligands:
        command = f"./generate_chain_layout.py {annot_file} {chain_layout_file} {prot} None {chain_domain_layout_files}"
    else:
        ligand_url = utils_ligand_links.link_ligands(prot)
        ligand_filename = f"{annot_dir}/ligands_data-{prot}.json"
        ret = utils_common.curl_file(ligand_url, new_dir, ligand_filename, option="--insecure")
        if not ret:
            sys.exit
        command = f"./generate_chain_layout.py {annot_file} {chain_layout_file} {prot} {ligand_filename} {chain_domain_layout_files}"
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



def put_to_dir(to_dir, annot_file, chain, prot):
    with open(annot_file, "r") as fd:
        data = json.load(fd)
        data_chain = copy.deepcopy(data)
        data_chain[p]["secondary_structure_elements"] = []
        for sse in data[p]["secondary_structure_elements"]:
            if sse["chain_id"] == chain:
                data_chain[p]["secondary_structure_elements"].append(sse)
    to_file = f"{to_dir}/{p}{chain}00-annotated.sses.json"
    with open(to_file, "w") as fd:
        json.dump(data_chain, fd, indent=4, separators=(",", ":"))


def append_to_overprot_input_file(trimer_file, protein, chain):
    with open(trimer_file, "r") as fd:
        data = json.load(fd)
    rec = {}
    rec["domain"] = str(protein) + str(name) + "00"
    rec["pdb"] = str(protein)
    rec["chain_id"] = chain[protein]
    rec["ranges"] = ":"
    rec["auth_chain_id"] = rec["chain_id"]
    rec["auth_ranges"] = rec["ranges"]
    data.append(rec)

    with (open(trimer_file, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ": "))


def append_to_overprot_whole_input_file(output_json, protein, chain, name):
    with open(output_json, "r") as fd:
        data = json.load(fd)
    rec = {}
    rec["domain"] = str(protein) + str(name) + "00"
    rec["pdb"] = str(protein)
    rec["chain_id"] = chain[protein]
    rec["ranges"] = ":"
    rec["auth_chain_id"] = rec["chain_id"]
    rec["auth_ranges"] = rec["ranges"]
    data.append(rec)

    with (open(output_json, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ": "))


def set_ubertemplate(out_dir, d_file, log_dir):

    # use overprot tool to download set ubertemplate for whole family

    # we have to call overprot for its directory,
    # and read activate file before we call overprot itself
    cwd = os.getcwd()
    os.chdir("../overprot/OverProtCore")

    activ_co = "source venv/bin/activate"
    overprot_co = f"python3 overprot.py 1.1.1.1 {cwd}/{out_dir} --domains {cwd}/{d_file}"
    command = f"bash -c '{activ_co}; {overprot_co}'"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = f"{log_dir}/overprot-fail-{time.time()}.txt"
        utils_common.error_report(log_dir, command, out, err, log_fn)
        logging.error(f"overprot on 1.1.1.1 failed, report save to {log_fn}")
    os.chdir(cwd)
    return r


def all_multiple_diagrams(l, directory, layout_list, not_ligands):
    list = (l.split("/"))[-1]
    logging.debug("Try to sync with the old version")
    svg_args = argparse.Namespace(
        add_descriptions="multi",
        layouts=[f"{directory}"],
        ligands=False,
        opac='1',
        debug=False,
        irotation=None,
        itemplate=None,
        iligands=None,
        orotation=f"{directory}/rotation.json",
        otemplate=f"{directory}/{list}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{directory}/{list}.svg",
        data=None,
    )
    if not_ligands:
        svg_args.ligands = False
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{directory}")
    multiple_diagrams(directory, directory, True, not_ligands, name)


def multiple_diagrams(detection_dir, layout_dir, first_svg, not_ligands, name):
    svg_args = argparse.Namespace(
        add_descriptions="multi",
        layouts=[f"{layout_dir}"],
        ligands=False,
        opac='1',
        debug=False,
        irotation=f"{detection_dir}/rotation.json",
        itemplate=None,
        iligands=f"{detection_dir}/ligands.json",
        orotation=None,
        otemplate=f"{detection_dir}/{name}.json",
        output=f"{layout_dir}/{name}.svg",
        oligands=None,
        data=None,
    )
    if not_ligands:
        svg_args.ligands = False
    if first_svg is False:
        utils_generate_svg.compute(svg_args)
    svg_args.otemplate = None
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}", "")
    svg_args.output = f"{layout_dir}/{name}_3.svg"
    svg_args.itemplate = f"{detection_dir}/{name}.json"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}", "")
    svg_args.output = f"{layout_dir}/{name}_2.svg"
    svg_args.opac = "0"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}", "")
    svg_args.output = f"{layout_dir}/{name}_4.svg"
    svg_args.add_descriptions = "multi,ligands"
    svg_args.ligands = True
    svg_args.oligands = f"{layout_dir}/ligands.json"
    svg_args.opac = "0"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, svg_args.output[:-4]+".png", f"{detection_dir}", "")


# read arguments
parser = argparse.ArgumentParser(
    description='generate single and multiple diagrams for a list of spike protein chains')
parser.add_argument('--protein_list', required=True, action="store",
                    help="which will be put the result alphafold directory")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('-n', '--not-ligands', action='store_true',
                    help='result diagrams does not contain ligands')
parser.add_argument('--family_dbs', action='store', required=True,
                    help='"dev" - build against development version, or "2dprots" - build against actual 2dprots dbs')
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')


args = parser.parse_args()


# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
pl_split = args.protein_list.split("/")
name = pl_split[-1]
new_dir = f"generated-{name}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)

# start logger
logger = start_logger(args.debug, args.info, new_dir, name)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# we want to download read relevant cif files
protein_list = set()
with open(args.protein_list) as fd:
    line = fd.readline()
    while line:
        protein_list.add(line[0:4])
        line = fd.readline()

cif_dir = f"{new_dir}/cif"
os.makedirs(cif_dir, exist_ok=True)
utils_pdbe_api.download_protein_list_cif(cif_dir, protein_list)
detect_dir = f"{new_dir}/detect"
os.mkdir(detect_dir)

# we want create detection files
for p in protein_list:
    cif_file = f"{cif_dir}/{p}.cif"
    utils_annotation.detect_list(detect_dir, cif_dir, [p], new_dir, whole_protein="pf")
    fn = f"{detect_dir}/{p}-detected.sses.json"
    p_data = {}

# find out chain -> domain mapping
chains = {}
for p in protein_list:
    detect_file_name = f"{detect_dir}/{p}-detected.sses.json"
    with open(detect_file_name, "r") as fd:
        protein_detect_json = json.load(fd)
        # chains contains dict of all chains - their value will be the start and end residuum
        chains[p] = {}
        for sse in protein_detect_json[p]["secondary_structure_elements"]:
            if utils_sses_types.is_nontrivial_sse(sse["label"], sse["end"] - sse["start"]):
                if sse["chain_id"] not in chains[p]:
                    chains[p][sse["chain_id"]] = [[sse["start"], sse["end"]]]
                else:
                    chains[p][sse["chain_id"]].append([sse["start"], sse["end"]])
logging.info(f"INFO: Number of chains {len(chains)} {chains}")

protein_list_omit = set()
fam = {}
trans = {}
for p in protein_list:
    fam[p], trans[p] = utils_pdbe_api.get_protein_chain_domains_and_trans(p, chains[p])
    for ch in fam[p]:
        if "Rest" in fam[p][ch]:
            protein_list_omit.add(p)

pc_numbers = {}
for p in protein_list:
    if len(fam[p]) not in pc_numbers:
        pc_numbers[len(fam[p])] = 1
    else:
        pc_numbers[len(fam[p])] += 1

max_number = None
max_number_value = 0
for n in pc_numbers:
    if (not max_number) or (pc_numbers[n] > max_number_value):
        max_number = n
        max_number_value = pc_numbers[n]

for p in protein_list:
    if len(fam[p]) != max_number:
        protein_list_omit.add(p)
protein_list = protein_list - protein_list_omit

chain_sets = {}
first = True
for p in protein_list:
    number = 0
    if first:
        for c in fam[p]:
            chain_sets[str(number)] = {"etalon": fam[p][c], p: c}
            number += 1
    else:
        for n in chain_sets:
            for c in fam[p]:
                if chain_sets[str(n)]["etalon"] == fam[p][c]:
                    chain_sets[str(n)][p] = c
                    fam[p][c] = "None"
                    break
    first = False

logging.info(f"INFO: Number of chains {len(chains)} {chains}")

domains = {}

for p in protein_list:
    domains[p] = utils_pdbe_api.get_protein_domains(p)

layout_dir = f"{new_dir}/layouts"
os.makedirs(layout_dir, exist_ok=True)
for p in protein_list:
    for d in domains[p]:
        # for all domains find out whether we have computed its layout yet
        chain_new = d[4:-2]
        dom = d[-2:]
        fn = f"layout-{p}_{chain_new}{dom}.json"
        if args.family_dbs == "dev":
            wfn = f"http://147.251.21.23/files/domain/{p}{chain_new}{dom}/layout.json"
        elif args.family_dbs == "2dprots":
            wfn = f"http://2dprots.ncbr.muni.cz/files/domain/{p}{chain_new}{dom}/layout.json"
        else:
            print(f'wrong 2dprot family dbs specification ({args.family_dbs}) should be "dev" or "2dprots"')
            sys.exit()
        ret = utils_common.curl_file(wfn, layout_dir, f"{layout_dir}/{fn}")

for i in chain_sets:
    # create overprot
    overprot_input = f"{new_dir}/overprot_input_{i}.json"
    with open(overprot_input, "w") as fd:
        data = []
        json.dump(data, fd, indent=4, separators=(",", ": "))
    for protein in protein_list:
        append_to_overprot_whole_input_file(overprot_input, protein, chain_sets[i], f"{i}")
    set_ubertemplate(f"{new_dir}/overprot_{i}", overprot_input, new_dir)

    os.makedirs(f"{new_dir}/annotated_{i}", exist_ok=True)
    shutil.copyfile(f"{new_dir}/overprot_{i}/results/consensus.cif", f"{new_dir}/overprot_{i}/cif_cealign/consensus.cif")
    shutil.copyfile(f"{new_dir}/overprot_{i}/results/consensus.sses.json", f"{new_dir}/overprot_{i}/cif_cealign/consensus-template.sses.json")
    # use template to annotate all domains in this bunch
    for fn in os.listdir(f"{new_dir}/overprot_{i}/cif_cealign"):
        if fn.endswith(".sses.json") and len(fn) == 17:
            os.symlink(f"../../cif/{fn[:4]}.cif", f"{new_dir}/overprot_{i}/cif_cealign/{fn[:4]}.cif")
            utils_annotation.annotate_with_given_template_list(f"{new_dir}/overprot_{i}/cif_cealign/", [(fn[:4], chain_sets[fn[4]][fn[:4]], "", "")], f"{new_dir}", domain_boundaries=False, prefix_secstr_dir="../overprot/OverProt/SecStrAnnot2/")
            shutil.copyfile(f"{new_dir}/overprot_{i}/cif_cealign/{fn[:4]}{chain_sets[fn[4]][fn[:4]]}00-annotated.sses.json", f"{new_dir}/annotated_{i}/{fn[:4]}{chain_sets[fn[4]][fn[:4]]}00-annotated.sses.json")
            shutil.copyfile(f"{new_dir}/overprot_{i}/results/diagram.json", f"{new_dir}/annotated_{i}/diagram.json")

annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)

for p in protein_list:
    a_list = []
    for c in chain_sets:
        a_list.append(f"{new_dir}/annotated_{c}/{p}{chain_sets[c][p]}00-annotated.sses.json")
    ret = utils_aggregate_annotations.aggregate_annotations(a_list, f"{new_dir}/annotation/{p}-annotated.sses.json", detect_dir, chain_sets)

# compute layout for all chains
layout_chain_list = []
for p in protein_list:
    for chain in chains[p].keys():
        struct_asym_id = trans[p][chain][0]
        prefix = f"layout-{p}_{struct_asym_id}"
        chain_layout_file = f"{new_dir}/layout-{p}_{chain}.json"
        chain_domain_layout_files = ""
        lc = 0
        for filename in os.listdir(layout_dir):
            if filename.startswith(prefix) and \
               ((str(filename[len(prefix)]).isdigit()) or (filename[len(prefix)] == "r")):
                chain_domain_layout_files += f" {layout_dir}/{filename}"
                lc = lc + 1
        # generate layout even for 1 chain - we have to choose colors for only one and only rest layout as well
        cif_file = f"{cif_dir}/{p}.cif"
        for i in chain_sets:
            if chain_sets[i][p] == chain:
                chain_set_name = i
        generate_chain_layout(cif_file, chain_layout_file, chain_domain_layout_files, new_dir, f"{new_dir}/annotated_{chain_set_name}/{p}{chain}00-annotated.sses.json", p, True, annotation_dir)
        layout_chain_list.append(f"{p}_{chain}")
    last_time = utils_common.time_log("4 per chains layouts computed", last_time)


# copy statistics form diagram-{family}. json file to sses-stats.json statistics,
logging.info("INFO: Starting save statistics and counting average protein")
stats_file = new_dir + "/sses-stats.json"
d_list = []
for c in chain_sets:
    d_list.append(f"{new_dir}/annotated_{c}/diagram.json")

ret = utils_pf_stats.copy_and_rename_statistics(d_list, stats_file)

# find universal start layout
logging.info("INFO: Find universal start layout")
add_list = []
for p in protein_list:
    universal_start = p  # utils_find_univ_start.find_universal_start(add_list, annotation_dir, 7)
    add_list.append([p, "ALL", "ALL"])
if len(protein_list) == 0:
    print("no layout")
    sys.exit()

# TODO: rewrite
logging.info(f"INFO: Universal_start {universal_start}")

last_time = utils_common.time_log("3 (stats, average protein, universal start layout counted)", last_time)
# generate layout for universal start domain
logging.debug("generate universal start {universal_start} layout")
###############
# compute layout for whole protein together
full_layout_file = f"{new_dir}/layout-{universal_start}.json"
start_layout_filename = f"{new_dir}/start-layout.json"
full_layout_file = start_layout_filename
#chain_domain_layout_files = f"{layout_dir}/layout-{universal_start}_*.json"
# if (len(full_chains) == 1) and : # draw ligands as well - they are only in this phase of generation - TODO - is this faster??
#    chain_layout_file = f"{new_dir}/layout-{args.protein}_{list(full_chains)[0]}.json"
#    logging.debug(f"DEBUG: generate full layout - short proces, copy the only chain layout {chain_layout_file} to full one")
#    shutil.copyfile(chain_layout_file, full_layout_file)
# else:
chain_layout_files = f"{new_dir}/layout-{universal_start}_*.json"
generate_chain_layout(cif_file, full_layout_file, chain_layout_files, new_dir, f"{annotation_dir}/{universal_start}-annotated.sses.json", universal_start, True, annotation_dir)

#########
#res = generate_layout_list(annotation_dir, [[universal_start, "ALL", "ALL"]], None, None, None, new_dir, args.not_ligands)
#if res == []:
#   logging.critical("can't create start layout")
#    sys.exit()
#start_layout_filename = f"{new_dir}/start-layout.json"
#shutil.move(f"{new_dir}/layout-{universal_start[0:4]}ALLALL.json", start_layout_filename)
# generate layouts for all annotated_list
logging.debug("generate all layouts")
color_layout_filename = f"{annotation_dir}/consensus-template.sses.json"
layout_list = []
for p in protein_list:
    chain_layout_files = f"{new_dir}/layout-{p}_*.json"
    generate_chain_layout(cif_file, f"{new_dir}/layout-{p}.json", chain_layout_files, new_dir, f"{annotation_dir}/{p}-annotated.sses.json", p, True, annotation_dir)
    layout_list.append(p)

last_time = utils_common.time_log("4 (all layouts computed)", last_time)

# generate images ------------------------------------------
# if it is not update wget rotation.json and files from previous version of generating this family (to preserve the rotation)
# generate multiple images
all_multiple_diagrams(args.protein_list, new_dir, layout_list, args.not_ligands)

# generate single images
svg_args = argparse.Namespace(
    add_descriptions="chain",
    layouts=[f"{new_dir}"],
    ligands=True,
    opac='0',
    debug=False,
    irotation=f"{new_dir}/rotation.json",
    itemplate=None,
    iligands=None,#f"{new_dir}/ligands.json",
    orotation=None,
    otemplate=None,
    oligands=None,
    output=f"{new_dir}/{args.protein_list}_2.svg",
    data=False,
)
if args.not_ligands:
    svg_args.ligands = False

c = 0
svg_list = []
for p in layout_chain_list:
    svg_args.layouts = [f"{new_dir}/layout-{p}.json"]
    svg_args.data = None
    svg_args.ligands = False
    svg_args.iligands = None
    svg_args.output = f"{new_dir}/{p}.svg"
#    r = utils_generate_svg.compute(svg_args)
#    if r is not None:
#        c += 1
#        svg_list.append(i)
#    svg_args.layouts = [f"{new_dir}/layout-{pcd}.json"]
#    svg_args.data = f"{new_dir}/{prot}_{chain}{domain}.json"
#    svg_args.ligands = True
#    svg_args.iligands = f"{new_dir}/ligands.json"
#    svg_args.output = f"{new_dir}/{prot}_{chain}{domain}l.svg"
#    r = utils_generate_svg.compute(svg_args)

for p in layout_list:
    svg_args.layouts = [f"{new_dir}/layout-{p}.json"]
    svg_args.data = None
    svg_args.ligands = False
    svg_args.iligands = None
    svg_args.output = f"{new_dir}/{p}.svg"
    r = utils_generate_svg.compute(svg_args)
    if r is not None:
        c += 1
        svg_list.append(i)


# final cleanup
now_time = time.time()
if (not args.debug) and (len(add_list) == c):
    shutil.rmtree(f"{annotation_dir}")
    if os.path.isfile(f"{new_dir}/rotation2.json"):
        os.remove(f"{new_dir}/rotation2.json")
    if os.path.isfile(f"{new_dir}/rotation_new.json"):
        os.remove(f"{new_dir}/rotation_new.json")
    if os.path.isfile(f"{new_dir}/rotation_aux.json"):
        os.remove(f"{new_dir}/rotation_aux.json")
    if os.path.isfile(f"{new_dir}/sses-stats.json"):
        os.remove(f"{new_dir}/sses-stats.json")

    last_time = utils_common.time_log(f"5 OK all diagrams created {last_time}", last_time)
    last_time = utils_common.time_log(f"6 OK sum {start_time}", start_time)
else:
    last_time = utils_common.time_log(f"5 BUG all diagrams created {last_time}", last_time)
    last_time = utils_common.time_log(f"6 BUG sum {start_time}", start_time)

if (len(add_list) == c):
    with open(f"{new_dir}/domain_list.txt", "w") as fdl:
        for i in svg_list:
            fdl.write(f"{i}\n")

print("generating succesfully complete\n")

