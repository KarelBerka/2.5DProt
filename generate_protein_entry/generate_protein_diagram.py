#!/usr/bin/env python3
# tool which generate for given protein its directory with
# * protein images
# and image per chain


import logging
import sys
import argparse
import time
import os
import shutil
import subprocess
import json
import copy

import utils_common
import utils_annotation
import utils_pdbe_api
import utils_read_pdb
import utils_read_cif
import utils_sses_buffer
import utils_ligands_input

sys.path.append('..')
import utils_generate_svg
import utils_links


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


def copy_ligand_json(ligand_filename, args, annotation_dir, new_dir):
    # we want to download read relevant ligand cif file
    if args.without_ligands:
        return None

    if (args.ligands_file):
        logging.info(f"INFO: Copy ligand file from {args.ligands_file}")
        ligand_url = args.ligands_file
    else:
        ligand_url = utils_links.link_ligands(args.protein)

    logging.info("INFO: Download ligand json data")
    if not utils_common.curl_file(ligand_url, new_dir, ligand_filename, option="--insecure"):
        # we can't download pre computed ligands_data file
        # we have to create it
        cif_file = f"{annotation_dir}/{args.protein}.cif"
        utils_ligands_input.create_ligands_input(cif_file, annotation_dir, ligand_filename)
#        utils_common.error_report(new_dir, "copy ligand json problem - see curl fail", None, "", f"{new_dir}/missing_ligand_file.txt")
#        logging.error("generate_domain_layout_from_scratch.py failed - report save to missing_ligand_file.txt")
#        return 1
    return 0


def copy_pdb_cif(args, annotation_dir, new_dir):
    # we want to download read relevant pdb file

    pdb_filename = f"{annotation_dir}/{args.protein}.pdb"
    logging.info(f"INFO: Copy pdb file from {args.pdb_file}")
    pdb_url = args.pdb_file

    if not utils_common.download_file(pdb_url, new_dir, pdb_filename):
        utils_common.error_report(new_dir, "copy pdb file problem - see curl fail", None, "", f"{new_dir}/missing_pdb_file.txt")
        logging.error("generate failed - report save to missing_pdb_file.txt")
        print("input protein structure file can't be downloaded", file=sys.stderr)
        sys.exit(1)
    shutil.copyfile(f"{annotation_dir}/{args.protein}.pdb", f"{annotation_dir}/{args.protein}.cif")
    return pdb_filename


#def copy_cif(args, annotation_dir, new_dir):
#    # we want to download read relevant cif file
#
#    cif_filename = f"{annotation_dir}/{args.protein}.cif"
#    if args.cif_file:
#        logging.info(f"INFO: Copy cif file from {args.cif_file}")
#        cif_url = args.cif_file
#    else:
#        cif_url = utils_links.link_cif(args.protein)
#
#    logging.info("INFO: Download ligand json data")
#    if not utils_common.curl_file(cif_url, new_dir, cif_filename, option="--insecure"):
#        utils_common.error_report(new_dir, "copy cif file problem - see curl fail", None, "", f"{new_dir}/missing_cif_file.txt")
#        logging.error("generate failed - report save to missing_cif_file.txt")
#        print("input cif file can't be downloaded")
#        sys.exit(2)
#    return cif_filename
#

def agregate_ligands_and_mo_points(ligands_file, mo_points_file, annotation_dir, new_dir):
    # we want to agregate mo_point_file to ligands_file
    mpf_local = mo_points_file
    logging.info("INFO: Download mop json data")
    if not utils_common.curl_file(mo_points_file, new_dir, mpf_local, option="--insecure"):
        utils_common.error_report(new_dir, "copy mop file problem - see curl fail", None, "", f"{new_dir}/missing_mop_file.txt")
        logging.error("generate failed - report save to missing_mop_file.txt", file=sys.stderr)
        return ligands_file
    # now we have mol_online points local
    # read it and convert it to ligands
    ligands_old = f"{ligands_file}.orig"
    shutil.move(ligands_file, ligands_old)

    with open(mpf_local, "r") as fm:
        data_mpf = json.load(fm)
        with open(ligands_old, "r") as fl:
           data_lig = json.load(fl)

    cc = 0
    for rec in data_mpf:
        d = {"chain": "?",
             "auth_chain": "?",
             "n_atoms": rec['n_atoms'],
             "center": rec['center'],
             "entity": f"s{cc}",
             "entity_type": rec['entity_type'],
             "entity_comp": rec['entity_comp'],
             "entity_description": rec['entity_description']
             }
        data_lig[f"s{cc}"]=d
        cc += 1
    with open(ligands_file, "w") as fo:
        json.dump(data_lig, fo, separators=(',', ':'))
    return ligands_file


def copy_channel_json_to_dir(args, annotation_dir, new_dir):
    # we want to download read relevant channel file
    channels_filename = f"{annotation_dir}/channel-{args.protein}.json"

    logging.info("INFO: Download channel json data")
    # we have to curl file
    # find the url location
    if args.channels_file:
        shutil.copyfile(f"{args.channels_file}", channels_filename)
        return channels_filename, True
    else:
     channels_url = utils_links.link_channels(args.protein)

    # curl ligand file
    if not utils_common.curl_file(channels_url, new_dir, channels_filename, option="--insecure"):
        utils_common.error_report(new_dir, "copy channel json problem - see curl fail", None, "", f"{new_dir}/missing_channel_file.txt")
        logging.error("generate_domain_layout_from_scratch.py failed - report save to missing_channel_file.txt")
        if args.mol_online:
            print("input channel file can't be downloaded", file=sys.stderr)
            sys.exit(3)
        else:
            shutil.copy("channel.json", channels_filename)
        channel = False
    else:
        channel = True
    return channels_filename, channel


def copy_cif_to_dir(args, cif_dir):
    # we want to download read relevant ligand cif files
    if (args.cif_dir):
        logging.info("INFO: Copy protein cif")
        shutil.copyfile(f"{args.cif_dir}/{args.protein}.cif", ligand_filename)
    else:
        logging.info("INFO: Download protein cif")
        ret = utils_pdbe_api.download_protein_cif(cif_dir, args.protein)
        if ret is None:
            sys.exit()


def download_all_precomputed_domains(args, new_dir):
    data = {}
    try:
        chains = utils_pdbe_api.get_protein_chains(args.protein)
    except:
        # if there is no chain, there still could be a set of ligands
        log_fn = f"{new_dir}/missing_protein_fail.txt"
        # the diagram without ligands is empty - so copy the empty diagram
        logging.error(f"protein not in https://www.ebi.ac.uk/pdbe/mappings - report save to {log_fn}")
        utils_common.error_report(new_dir, "protein not in https://www.ebi.ac.uk/pdbe/mappings or has no domain", None, "", log_fn)
        # and continue with the ligands
        chains = set()

    logging.info(f"INFO: Number of chains {len(chains)} {chains}")
    for chain in chains:
        data[chain] = ""
        chain_new = utils_pdbe_api.transform_cath_struct_asym_id_to_chain_id(args.protein, chain)

        # for all chains find its domians

        domains = utils_pdbe_api.get_protein_chain_domains(args.protein, chain)
        for dom in domains:
            # for all domains find out whether we have computed its layout yet
            fn = f"layout-{args.protein}_{chain}{dom}.json"
            if args.family_dbs:
                if args.family_dbs == "dev":
                    wfn = f"http://147.251.21.23/files/domain/{args.protein}{chain_new}{dom}/layout.json"
                else: # args.family_dbs == "2dprots":
                    wfn = f"http://2dprots.ncbr.muni.cz/files/domain/{args.protein}{chain_new}{dom}/layout.json"
#                else:
#                    wfn = f"{args.family_dbs}/domain/{args.protein}{chain_new}{dom}/layout.json"
            else:
                wfn = f"http://2dprots.ncbr.muni.cz/files/domain/{args.protein}{chain_new}{dom}/layout.json"
            ret = utils_common.curl_file(wfn, layout_dir, f"{layout_dir}/{fn}")
            if ret:
                # we have computed them thus we have to omit this interval in the notcomputed section
                db = utils_pdbe_api.get_domain_boundaries(args.protein, chain, dom)
                # 3d coordinates can't be processed and changed  - use the original per protein one
                fix_sses_3d_coordinates(f"{layout_dir}/{fn}", detect_prot_file)
                if data[chain] != "":
                    data[chain] = data[chain]+","+db
                else:
                    data[chain] = db
    return data, chains


def generate_domain_layout(protein, chain, working_dir, annotation_dir, layout_dir):
    # compute single per domain layout
    layout_file = f"{layout_dir}/layout-{protein}_{chain}re.json"
    annotation_file = f"{annotation_dir}/{protein}_{chain}re-annotated.sses.json"
    command = f"./generate_domain_layout_from_scratch.py {annotation_file} {chain} {layout_file} {protein}"
    logging.debug(f"DEBUG: Generating layout for {protein} {chain} re - the rest of chain which is not in any domain computed yet")
    logging.info(f"INFO: Command {command}")
    p = subprocess.Popen(command, shell=True,  stdout=subprocess.PIPE,
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


def generate_chain_layout(chain_layout_file, chain_domain_layout_files, working_dir, annot_file, prot, annot_dir, channel_filename, ligands_filename, mol_online):
    # aggregate several layouts to one multi-layout - either domains to chain or chain to protein ones
    command = f"./generate_chain_layout.py --annotation_file {annot_file} --output {chain_layout_file} --protein {prot} --chains {chain_domain_layout_files}"
    if ligands_filename != None:
         command += f" --ligands {ligand_filename}"
    if channel_filename != None:
         command += f" --channel {channel_filename}"
    if mol_online:
         command += " --mol_online"
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


def generate_sses_ligands_channel_with_sl_layout_file(full_layout_file, sl_layout, new_dir, annotation_dir, prot, channel_filename):
    # aggregate several layouts to one multi-layout - either domains to chain or chain to protein ones
    command = f"./generate_slc_layout_from_sl_layout.py --output {full_layout_file} --protein {prot} --sl_layout {sl_layout} --annotation_dir {annotation_dir}"
    if channel_filename != None:
         command += f" --channel {channel_filename}"
    logging.debug(f"DEBUG: generating chain layout {full_layout_file}")
    logging.debug(f"DEBUG: command {command}")

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0) or (not os.path.isfile(full_layout_file)):
        log_fn = f"{new_dir}/generate_chain_layout_fs-fail.txt"
        utils_common.error_report(new_dir, command, None, err, log_fn)
        logging.error(f"generate_domain_layout.py failed - report save to {log_fn}")
        sys.exit()
        return False
    else:
        return True


def fix_sses_3d_coordinates(annot_fn, detect_fn):
    with open(annot_fn, "r") as fa, open(detect_fn, "r") as fd:
        adata = json.load(fa)
        ddata = json.load(fd)

    adata_hr = adata['helices_residues']
    ddata_sses = ddata[list(ddata.keys())[0]]["secondary_structure_elements"]
    for sse in adata_hr:
        hr = adata_hr[sse]
        for dsse in ddata_sses:
            if ((hr[0] == dsse["start"] and hr[1] == dsse["end"]) or \
                (hr[1] == dsse["start"] and hr[0] == dsse["end"])) and \
                adata["chain_id"][sse] == dsse["chain_id"]:
                adata['3dlayout'][sse] = [dsse["start_vector"], dsse["end_vector"]]

    with open(annot_fn, "w") as fa:
        json.dump(adata, fa, indent=4, separators=(",", ": "))


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a protein structure')
parser.add_argument('--protein', required=True, action="store",
                    help="generate diagram database for given protein")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result protein directory")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')
parser.add_argument('--family_dbs', action='store',
                    help='"dev" - build against development version, or "2dprots" - build against actual 2dprots dbs')
parser.add_argument('--cif_dir', action='store',
                    help='set the directory containing protein cif files used for input')
parser.add_argument('--ligands_json_dir', action='store',
                    help='set the directory containing protein cif files used for input')

# --- mol_online ---
# mol_online have to be set with either --cif_file or --pdb_file
parser.add_argument('--mol_online', action='store_true',
                    help='2DProt calculation for mol_online tool - omit processes and outputs, which are not necessary')
parser.add_argument('--pdb_file', action='store',
                    help='set the pdb input file, ')
parser.add_argument('--cif_file', action='store',
                    help='set the cif input file')
parser.add_argument('--sses_layout_buffer', action='store',
                    help='use SSEs layout buffer to store SSEs positions for given input pdb/cif file')

# mo_poins_file is optional and it is a set of points in 3d in protein model space which will be
# displayed like ligands (usually some special possible tunel start points
parser.add_argument('--mo_points_file', action='store',
                    help='set the ligands input file')
parser.add_argument('--channels_file', action='store',
                    help='set the channels input file')

# --- ligands options ---
# by defaut ligands are read from utils_links.link_ligands link
# but they can set directly
parser.add_argument('--ligands_file', action='store',
                    help='set the ligands file')

# if there is no need to calculate ligands at all use this option
parser.add_argument('--without_ligands', action='store_true',
                    help='do not create ligands output')

args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
new_dir = f"{args.protein[0]}/{args.protein[1]}/generated-{args.protein}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"

if (args.mol_online) and (not args.pdb_file) and (not args.cif_file):
    print("--mol_online needs --pdb_file of cif_file option setting")
    print("wrong argument format", file=sys.stderr)
    sys.exit(20)

#if (args.pdb_file) and not (args.mol_online):
#    print("--pdb_file option needs --mol_online tag set")
#    sys.exit()

if (args.mo_points_file) and not (args.mol_online):
    print("--mo_points_file option needs --mol_online tag set")
    print("wrong argument format", file=sys.stderr)
    sys.exit(20)

if (args.sses_layout_buffer) and not (args.mol_online):
    print("--sses_layout_buffer option needs --mol_online tag set")
    print("wrong argument format", file=sys.stderr)
    sys.exit(20)

#if (args.cif_file) and  (args.pdb_file):
#    print("--cif_file and --pdb_file can't be combined")
#    sys.exit()

os.makedirs(new_dir, exist_ok=True)
annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)
cif_dir = annotation_dir
layout_dir = f"{new_dir}/layouts"
os.mkdir(layout_dir)

# start logger
logger = start_logger(args.debug, args.info, new_dir, args.protein)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

# Detect SSEs for the input protein
logging.info("INFO: Detect list starts")
detect_prot_file = f"{annotation_dir}/{args.protein}-detected.sses.json"
ligand_filename = f"{annotation_dir}/ligands_data-{args.protein}.json"
sses_ligands_layout_file = False

if (args.mol_online):
    if (args.pdb_file):
        # pdb_file was set
        pdb_filename = copy_pdb_cif(args, annotation_dir, new_dir)
        # create -detected.sses.json file from the given pdb one
        data = {}

        # if we use sses_layout_buffer we have to test whether we already do not conpute this layout
        if (args.sses_layout_buffer):
            sses_ligands_layout_file = utils_sses_buffer.is_file_in_buffer(args.protein, args.sses_layout_buffer, pdb_filename)
        if not sses_ligands_layout_file:
            # we do not have precomputed SSEs and ligands data
            # we have to generate ligand json 
            if utils_read_pdb.is_pdb(pdb_filename):
                try:
                    chains = utils_read_pdb.detect_sses(pdb_filename, detect_prot_file, args.protein)
                    utils_read_pdb.generate_ligand_json(pdb_filename, ligand_filename)
                except:
                    logging.exception("wrong protein structure file format")
                    print("wrong pdb protein structure file format", file=sys.stderr)
                    sys.exit(6)
            else:
                try:
                    utils_read_cif.generate_ligand_json(pdb_filename, ligand_filename)
                except:
                    logging.exception("wrong protein structure file format")
                    print("wrong cif protein structure file format", file=sys.stderr)
                    sys.exit(6)
                # we have to generate annotation file as well
                # create -detected.sses.json file from the given pdb one
                ret = utils_annotation.detect_list(annotation_dir, cif_dir, [args.protein], new_dir, whole_protein=True)
                if (ret == 0):
                    print("wrong cif protein structure file format", file=sys.stderr)
                    sys.exit(6)
                chains = set()
                with open(detect_prot_file, "r") as fd:
                    data = json.load(fd)
                    for rec in data[list(data.keys())[0]]["secondary_structure_elements"]:
                         chains.add(rec['chain_id'])

    else:
        print("--pdb_file have to be set", file=sys.stderr)
        sys.exit(20)

    if (args.mo_points_file):
        # mol_online special points was set
        # agregate mo_points file and ligand_file together
        agregate_ligands_and_mo_points(ligand_filename, args.mo_points_file, annotation_dir, new_dir)

else:
    # cif file can be parsed by SecStrAnnotator
    copy_cif_to_dir(args, cif_dir)
    ret = utils_annotation.detect_list(annotation_dir, cif_dir, [args.protein], new_dir, whole_protein=True)
    if (ret == 0):
        sys.exit()

    # go through all chains of the input protein and aggregate all domains
    # which have the diagram already comuted in args.family_dbs database
    data, chains = download_all_precomputed_domains(args, new_dir)

    # get relevant ligand json and cif file - we have to check all relevant ligand options as well
    copy_ligand_json(ligand_filename, args, annotation_dir, new_dir)

channel_filename, channel_presence = copy_channel_json_to_dir(args, annotation_dir, new_dir)
last_time = utils_common.time_log("1 all ssed detected, ligands find, channels downloaded", last_time)

if sses_ligands_layout_file == False:
    # we have to compue SSEs layout
    missing_sses = {}
    all_sses = {}
    # find all noncomputed sections (remove computed ones from the whole)
    # and put them to the rest annotation
    full_chains = set(chains)
    for chain in full_chains:
        missing_sses[chain] = []
        all_sses[chain] = []

    with open(detect_prot_file, "r") as fp:
        prot_data = json.load(fp)
        prot_all = copy.deepcopy(prot_data)
        for sse in prot_data[args.protein]["secondary_structure_elements"]:
            if sse["chain_id"] not in missing_sses:
                missing_sses[sse['chain_id']] = []
                all_sses[sse['chain_id']] = []
            if (sse["chain_id"] in data) and (len(data[sse["chain_id"]]) != 0):
                sint = data[sse["chain_id"]].split(",")
                if sse['chain_id'] not in missing_sses:
                    all_sses[sse['chain_id']] = [sse]
                else:
                    all_sses[sse['chain_id']].append(sse)
                if sse['start'] != sse['end']:
                    is_in = False
                    for interval in sint:
                        inse = interval.split(":")
                        if interval != '' and int(inse[0]) <= sse["start"] and sse["end"] <= int(inse[1]):
                            is_in = True
                    if is_in is False:
                        if sse['chain_id'] not in missing_sses:
                            missing_sses[sse['chain_id']] = [sse]
                        else:
                            missing_sses[sse['chain_id']].append(sse)
            else:
                missing_sses[sse['chain_id']].append(sse)
    prot_data[args.protein]["secondary_structure_elements"] = missing_sses
    prot_all[args.protein]["secondary_structure_elements"] = all_sses

    full_chains = missing_sses.keys()
    if len(full_chains) == 0:
        # if the protein has no chain, no need to compute anything
        logging.info("INFO: The pdbe api protein file is empty - it has no chain with sse")
        with open(f"{new_dir}/missing_data_fail.txt", "w") as f:
            f.close()
            last_lime = utils_common.time_log("8 OK", start_time)
        # draw an empty diagram
        shutil.copyfile("missing_data.svg", f"{new_dir}/{args.protein}.svg")
        # continue to draw ligand diagram which could be notempty

    data = []
    for chain in full_chains:
        rec = [args.protein, f"{args.protein}{chain}", chain, ":"]
        data.append(rec)
    with open(f"{new_dir}/chains.json", "w") as f:
        json.dump(data, f, indent=4, separators=(',', ': '))

    try:
        # for all chains annotate the non-layout sses and draw layout for them
        for i in full_chains:
            prot_data[args.protein]["secondary_structure_elements"] = missing_sses[i]
            prot_data[args.protein]["hydrogen_bonds"] = []
            fn = f"{annotation_dir}/{args.protein}_{i}re-annotated.sses.json"
            with open(fn, "w") as fre:
                json.dump(prot_data, fre, indent=4, separators=(',', ': '))
            prot_data[args.protein]["secondary_structure_elements"] = all_sses[i]
            fn = f"{annotation_dir}/{args.protein}_{i}-annotated.sses.json"
            with open(fn, "w") as fre:
                json.dump(prot_data, fre, indent=4, separators=(',', ': '))
            logging.info(f"INFO: generating rest annotation file {fn}")
            generate_domain_layout(args.protein, i, new_dir, annotation_dir, layout_dir)
            last_time = utils_common.time_log("3 rest layouts computed", last_time)
    except:
        print("chain rest parts generating problem", file=sys.stderr)
        sys.exit(7)
# compute layout for all chains

try:
  if not sses_ligands_layout_file:
    # mol_online with sl layout file was not set or layout was not computed yet
    # thus standart generation proces have to start
    for chain in full_chains:
        prefix = f"layout-{args.protein}_{chain}"
        chain_layout_file = f"{new_dir}/layout-{args.protein}_{chain}.json"
        chain_domain_layout_files = ""
        lc = 0
        for filename in os.listdir(layout_dir):
            if filename.startswith(prefix) and \
               ((str(filename[len(prefix)]).isdigit()) or (filename[len(prefix)] == "r")):
                chain_domain_layout_files += f" {layout_dir}/{filename}"
                lc = lc + 1
        # generate layout even for 1 chain - we have to choose colors for only one and only rest layout as well
        generate_chain_layout(chain_layout_file, chain_domain_layout_files, new_dir, f"{annotation_dir}/{args.protein}_{chain}-annotated.sses.json", args.protein, annotation_dir, None, None, args.mol_online)
except:
    print("chains generating problem", file=sys.stderr)
    sys.exit(8)

last_time = utils_common.time_log("4 per chains layouts computed", last_time)
# compute layout for whole protein together
full_layout_file = f"{new_dir}/layout-{args.protein}.json"
if sses_ligands_layout_file:
    generate_sses_ligands_channel_with_sl_layout_file(full_layout_file, sses_ligands_layout_file, new_dir, annotation_dir, args.protein, channel_filename)
else:
    chain_layout_files = f"{new_dir}/layout-{args.protein}_*.json"
    try:
        generate_chain_layout(full_layout_file, chain_layout_files, new_dir, detect_prot_file, args.protein, annotation_dir, channel_filename, ligand_filename, args.mol_online)
    except:
        print("full layout generating problem", file=sys.stderr)
        sys.exit(9)

last_time = utils_common.time_log("5 all domains and chains layouts computed", last_time)

# draw single image for whole protein
svg_args = argparse.Namespace(
    add_descriptions="protein",
    layouts=[f"{new_dir}/layout-{args.protein}.json"],
    ligands=True,
    transparency=False,
    debug=False,
    irotation=None,
    itemplate=None,
    iligands=None,
    otemplate=None,
    orotation=None,
    oligands=None,
    output=f"{new_dir}/{args.protein}.svg",
    data=None,
    channels=False,
    only="sses, channels, ligands",
    fill=True,
    groups=True,
)

utils_generate_svg.compute(svg_args)

try:
    utils_generate_svg.compute(svg_args)
except:
    print("diagram generating problem", file=sys.stderr)
    sys.exit(10)

last_time = utils_common.time_log("8 OK", start_time)
