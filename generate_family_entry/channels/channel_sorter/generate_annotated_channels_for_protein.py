#!/usr/bin/env python3
# tool which generate for given protein list and mole config
# the set of tunnels and try to parse them (for 1.10.630.10, 1.20.1070.10 are configuration files tuned)

import logging
import sys
import argparse
import time
import os
import shutil
import subprocess
import json

import utils_common
import utils_pdbe_api
import utils_af_api
import utils_read_cif
import utils_tunnel_pdb

import utils_1_20_1070_10
import utils_1_10_630_10
import utils_PIN
import utils_SLC2A3
import utils_P02920

import utils_molespecypic_output


def split_pcd(line):
    # split line record - (3zevA00, 6za8AAA, 1akd, P08684)  to protein, chain, domain list
    # remove obsolete newline
    if line[-1] == "\n":
        pcd_string = line[:-1]
    else:
        pcd_string = line

    if pcd_string[0].isdigit():
        # the record is Pdb domain
        counter1 = 0
        while counter1 < len(pcd_string):
            if pcd_string[counter1].isupper():
                break
            else:
                counter1 += 1
        if counter1 == len(pcd_string):
            # there is no chain tag add A as a default
            pcd = [pcd_string, "A", "00", pcd_string]
        else:
            counter2 = 0
            while counter2 < len(pcd_string[counter1:]):
                if pcd_string[counter1+counter2].isdigit():
                    break
                else:
                    counter2 += 1
            if counter2 == len(pcd_string[counter1:]):
                # there is no chain domain tag add 00 as a default
                pcd = [pcd_string[:counter1], pcd_string[counter1:], "00", pcd_string[:counter1]]
            else:
                # all tags are present
                pcd = [pcd_string[:counter1], pcd_string[counter1:counter1+counter2], pcd_string[counter1+counter2:], pcd_string[:counter1]]
    elif pcd_string[0].isupper():
        # the record is AF domain - chain "A", domain "00" should be added
        pcd = [pcd_string, "A", "00", pcd_string]
    else:
        # the record is unknown domain
        print(f"{pcd_string} is not PDB and AF domain - do not start to digit or big letter")
        sys.exit()
    return pcd


def read_input_pdb_cifs(args, annotation_dir):
    # if we have mixed or pdb input, we use pdb dbs

    pdb_pcds = []                 # found pcds
    not_in_pdb_dbs = []           # not found proteins

    logging.info("INFO: reading input cifs from pdb database")
    # read protein list
    with (open(args.protein_list,"r")) as fpl:
        line = fpl.readline()
        while line:
            # split line to protein, chain, domain tags
            pcd = split_pcd(line)
            # try to download protein from pdb
            ret = copy_or_download_cif_from_pdb_to_dir(args, annotation_dir, protein_name=pcd[0])
            if ret == False:
                # protein was not found - put it on notfound list
                not_in_pdb_dbs.append(pcd)
            else:
                # protein was found - write domain to pcd list
                pdb_pcds.append(pcd)
            line = fpl.readline()

    return pdb_pcds, not_in_pdb_dbs


def read_input_af_cifs(args, annotation_dir, not_in_pdb_dbs):
    # if we have mixed or af input, we use af dbs
    af_pcds = []
    not_in_af_dbs = []           # not found proteins

    logging.info("INFO: reading input cifs from af database")

    if not_in_pdb_dbs == []:
        # if we have mixed set we have not_pdb_protein list already downloaded
        # otherwise we have to download it now
        with (open(args.protein_list,"r")) as fpl:
            # download whole protein list
            line = fpl.readline()
            while line:
                pcd = split_pcd(line)
                not_in_pdb_dbs.append(pcd)
                line = fpl.readline()

    # try to download all proteins in af dbs
    for pcd in not_in_pdb_dbs:
        ret = copy_or_download_cif_from_af_to_dir(args, annotation_dir, protein_name=pcd[0])
        if ret == False:
            not_in_af_dbs.append(pcd)
        else:
            # alphafold has no chain and domain tag
            af_pcds.append(pcd)
    return af_pcds, not_in_af_dbs


def read_input_cif_from_file(protein_file, annotation_dir):
    # copy cif from input file and create pcds structure
    # copy file
    file_name = os.path.basename(protein_file)
    split_file_name= os.path.splitext(file_name)
    file_path = f"{annotation_dir}/{file_name}"
    shutil.copy(protein_file, file_path)

    # create pcds structure
    pcds = []
    chain_list = utils_read_cif.find_all_auth_chains(file_path)
    for chain in chain_list:
        pcds.append([split_file_name[0], chain, "", split_file_name[0]])
    return pcds

def avg_3d_coord(three_d_list):
    x,y,z = 0,0,0
    counter = 0
    for i in three_d_list:
        x += i[0]
        y += i[1]
        z += i[2]
        counter += 1
    return [x/counter, y/counter, z/counter]


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


def copy_or_download_cif_from_pdb_to_dir(args, cif_dir, protein_name=""):
    # we want to download read relevant ligand cif files
    if (protein_name == ""):
        protein = args.protein
    else:
        protein = protein_name
    if (args.cif_dir):
        # copy cif from dir
        logging.info(f"INFO: Copy protein {protein} cif")
        if os.path.isfile(f"{args.cif_dir}/{protein}.cif"):
            shutil.copyfile(f"{args.cif_dir}/{protein}.cif", f"{cif_dir}/{protein}.cif")
        else:
            return False
    else:
        # download cif from net
        logging.info(f"INFO: Download protein {protein} cif")
        ret = utils_pdbe_api.download_protein_cif(cif_dir, protein)
        if ret is None:
            if args.type == "mix":
                # we don't find cif, but we have to continue - it could be in af dbs
                # just inform we are not succesfull
                return False
            else:
                # we only go through pdb dbs thus there is no chance to find this protein
                # stop here
                print(f"ERROR: cif file {protein} can't be found in pdb dbs" )
                sys.exit(1)
    return True


def copy_or_download_cif_from_af_to_dir(args, cif_dir, protein_name=""):
    # we want to download read relevant ligand cif files
    if (protein_name == ""):
        protein = args.protein
    else:
        protein = protein_name
    if (args.cif_dir):
        logging.info(f"INFO: Copy protein {protein} cif")
        shutil.copyfile(f"{args.cif_dir}/{protein}.cif", f"{cif_dir}/{protein}.cif")
    else:
        logging.info(f"INFO: Download protein {protein} cif")
        ret = utils_af_api.download_protein_cif(cif_dir, protein)
        if ret is None:
            print(f"ERROR: cif file {protein} can't be found in af dbs")
            sys.exit(1)


def create_overprote_input(overprot_file, pcds, pcd_type, annotation_dir, family=""):
    data = []
    counter = 0
    for pcd in pcds:
        # create overprot input json
        counter += 1
        rec = {}
        if pcd_type == "pdb" or pcd_type == "mix":
            # pdbmain record
            rec["domain"] = pcd[0] + pcd[1] + pcd[2]
            rec["pdb"] = pcd[0]
            # have to be transformed from label_asym (cath output) to auth_asym (overprot input)
            cif_file=f"{annotation_dir}/{pcd[0]}.cif"
            rec["chain_id"] = utils_read_cif.transform_auth_asym_id_to_label_asym_id(cif_file, pcd[1])
            rec["domain"] = pcd[0] + rec["chain_id"] + pcd[2]
            try:
                res_range = utils_pdbe_api.get_domain_boundaries(pcd[0], pcd[1], pcd[2], family)
                if res_range == "":
                    res_range = "0:20000"
            except:
                res_range = "0:20000"
            rec["ranges"] = res_range
            # transformed pcd from label_asym (cath output) to auth_asym (overprot input) as well
            pcd[1] = rec["chain_id"]
        else:
            # alphafold domain record
            rec["domain"] = str(pcd[0]) + "A00"
            rec["pdb"] = str(pcd[0])
            rec["chain_id"] = "A"
            rec["ranges"] = "0:20000"
        rec["auth_chain_id"] = rec["chain_id"]
        rec["auth_ranges"] = rec["ranges"]
        data.append(rec)
        if counter > 2000:
            break

    with (open(overprot_file, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ": "))


def set_ubertemplate(out_dir, d_file, log_dir, cif_dir, new_dir, fast=False):
    cwd = os.getcwd()
    # prepend cwd only for relative paths - absolute args (e.g. --output from
    # job runner) must stay untouched, os.path.join keeps absolute second arg
    out_dir = os.path.join(cwd, out_dir)
    d_file = os.path.join(cwd, d_file)
    cif_dir = os.path.join(cwd, cif_dir)
    os.chdir("../../../overprot/OverProtCore")

    activ_co = "source venv/bin/activate"
    if fast:
        overprot_co = f"python3 overprot.py 1.1.1.1 {out_dir} --sample_size 200 --domains {d_file} --structure_source file://{cif_dir}/"+"{pdb}.cif"
    else:
        overprot_co = f"python3 overprot.py 1.1.1.1 {out_dir} --sample_size 2000 --domains {d_file} --structure_source file://{cif_dir}/"+"{pdb}.cif"
    command = f"bash -c '{activ_co}; {overprot_co}'"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(new_dir), f"overprot-fail-{time.time()}.txt")
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"overprot on 1.1.1.1 failed, report save to {log_fn} r je {r} log_dir {log_dir}")
        logging.error(f"command {command}")
        logging.error(f"out {out}")
        logging.error(f"err {err}")

    os.chdir(cwd)
    return r


def is_domain_complete_1_10_630_10(model_annotation, prot_annotation, percentage=80):
    mapping = utils_1_10_630_10.sses_model_mapping_1og2(model_annotation)
    test = utils_1_10_630_10.test_domain_completeness(prot_annotation, mapping, percentage)
    return test


def is_domain_complete_1_20_1070_10(model_annotation, prot_annotation, percentage=80):
    mapping = utils_1_20_1070_10.sses_model_mapping_2vt4(model_annotation)
    test = utils_1_20_1070_10.test_domain_completeness(prot_annotation, mapping, percentage)
    return test


def compute_protein_tunnels(tunnels_dir, orig_xml_file_name, pcd, cif_file, start=None, end=None, middle=None):
    # use molonline to create protein tunels dbs

    xml_file_name=f"{tunnels_dir}/customized.xml"
    # at first customize mole input json file
    with open(orig_xml_file_name, "r") as fr, open(xml_file_name, "w") as fw:
        line = fr.readline()
        while line:
            if "PredictedPoint" in line:
                if middle != None and middle != 0:
                    fw.write(f'<Point X="{middle[0]}" Y="{middle[1]}" Z="{middle[2]}" />\n')
            if "Query" in line and start != None and start != 0:
                sline = line.split('"')
                if "HEM" in line:
                    fw.write(f'<Point X="{start[0]}" Y="{start[1]}" Z="{start[2]}" />\n')
                else:
                    auth_chain = utils_read_cif.transform_label_asym_id_to_auth_asym_id(cif_file, pcd[1])
                    ret = utils_read_cif.find_all_ligands_coords(cif_file, auth_chain, "HOH", [sline[1]])
                    if len(ret) > 0:
                        k = list(ret.keys())[0]
                        coord = avg_3d_coord(ret[k]["coords"])
                        fw.write(f'<Point X="{coord[0]}" Y="{coord[1]}" Z="{coord[2]}" />\n')
            elif "SPoint" in line:
                 if start != None and start != 0:
                    fw.write(f'<Point X="{start[0]}" Y="{start[1]}" Z="{start[2]}" />\n')
            elif "EPoint" in line:
                 if end != None and end != 0:
                    fw.write(f'<Point X="{end[0]}" Y="{end[1]}" Z="{end[2]}" />\n')
            elif "XXX" not in line:
                fw.write(line)
            else:
                if "XXXCIF" in line:
                    c = line.find("XXXCIF")
                    fw.write(line[:c]+f"{cif_file}"+line[c+6:])
                if "XXXID" in line:
                    c = line.find("XXXID")
                    fw.write(line[:c]+f"{pcd[0]}"+line[c+5:])
                if "XXXWD" in line:
                    c = line.find("XXXWD")
                    fw.write(line[:c]+f"{tunnels_dir}"+line[c+5:])
            line = fr.readline()

    # run moleonline on it
    command = f"mono ../../../Mole2_cmd/mole2.exe {xml_file_name}"
    print("command", command)
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    ret = p.wait()

    if (ret != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(new_dir), f"mole-fail-{time.time()}.txt")
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"mole on {pcd} failed , report save to {log_fn}")

    counter = 0
    for name in os.listdir(f"{tunnels_dir}/pdb/profile/"):
        if name.startswith("tunnel"):
            counter += 1
    return counter


def compute_1_10_630_10_protein_tunnels_in_several_rounds(new_dir, pcd, cif_file, start, tunnels_dir, tunnels_dir_w, args):
    # create directory for aggregation of all tunnels
    multi_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
    os.mkdir(multi_tunnels_dir)
    multi_tunnels_dir_w= f"{tunnels_dir_w}/{pcd[0]}{pcd[1]}{pcd[2]}"
    os.mkdir(multi_tunnels_dir_w)
    counter = 0
    if args.rounds == "6":
        round_set = ["1", "2", "3", "4", "5", "6"]
    elif args.rounds == "2":
        round_set = ["1", "6"]
    else:
        round_set = ["1", "2", "3", "4", "5", "6"]
    for i in round_set:
#    for i in ["1", "2", "3", "4", "5"]:
#    for i in ["1"]: #]#["3", "6"]:
        # compute tunnels for each setting (maximize the probability to find the tunnel)
        single_tunnels_dir= f"{new_dir}/tunnels{i}/{pcd[0]}{pcd[1]}{pcd[2]}"
        os.makedirs(single_tunnels_dir)
        compute_protein_tunnels(single_tunnels_dir, f"mole_config/1.10.630.10/orig{i}.xml", pcd, cif_file, start)
        # agregare solvent part tunnels (1-5) to one directory
        if i != "6":
            for fn in os.listdir(f"{single_tunnels_dir}/pdb/profile"):
                if fn.startswith("tunnel"):
                    counter +=1
                    shutil.copyfile(f"{single_tunnels_dir}/pdb/profile/{fn}", f"{multi_tunnels_dir}/{fn[:7]}{i}{fn[7:]}")
        else:
            for fn in os.listdir(f"{single_tunnels_dir}/pdb/profile"):
                if fn.startswith("tunnel"):
                    counter +=1
                    shutil.copyfile(f"{single_tunnels_dir}/pdb/profile/{fn}", f"{multi_tunnels_dir_w}/{fn[:7]}{i}{fn[7:]}")
    print(f"for {pcd[0]}{pcd[1]}{pcd[2]} found {counter} tunelu" )


def compute_PIN_protein_tunnels_in_several_rounds(limit, new_dir, annotation_dir, args):
    # PIN proteins have 4 start tunnel points and 5 end tunnel points
    # we have to create all possible tunnels
    tunnels_dir= f"{new_dir}/tunnels/"
    for i in [1, 2, 3, 4, 5]:
      for j in [1, 2, 3, 4, 5]:
        # start and end chosen
        print(f"pin st{i} <> end{j}")
        for pcd in limit:
            # for each domain - compute its tunnels
            # create file and directories names
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            single_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
            os.mkdir(single_tunnels_dir)
            # compute tunnel start and end
            start = utils_PIN.predict_PIN_tunnel_start(annotation_dir, pcd, i)
            end = utils_PIN.predict_PIN_tunnel_end(annotation_dir, pcd, j)
            # compute tunnels
            compute_protein_tunnels(single_tunnels_dir, "mole_config/PIN/whole_path.xml", pcd, cif_file, start=start, end=end)
            if not args.no_sort:
                sort_tunnels_unknown_family(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir)
        # clean up directory sstructure for the next round 
        # and save it in te proper place
        os.rename(f"{new_dir}/all_tunnels",f"{new_dir}/all_tunnels{i}{j}" )
        os.rename(f"{new_dir}/tunnels",f"{new_dir}/tunnels{i}{j}" )
        os.mkdir(tunnels_dir)
#        print("start, end", pcd,  start, end)
#        cif_file=f"/home/ivarekova/radka/bilkoviny/generate_family_entry/channels/channel_sorter/generated-tunnels-list_1.20.1070.10.small-2026-01-05T16_56_48/cifs/{pcd[0]}{pcd[1]}_rot.cif"
#        last_lime = utils_common.time_log("4 tunnels (in family PIN) computed", last_time)
#        print("4 tunnels (in family PIN) computed")
#        last_lime = utils_common.time_log("5 tunnels are sorted (in family PIN) computed", last_time)

#    print("5 tunnels are sorted (in family PIN) computed")


def compute_PIN_protein_tunnels_from_ligand(limit, new_dir, annotation_dir, args):
    # PIN proteins have 4 start tunnel points and 5 end tunnel points
    # we have to create all possible tunnels
    tunnels_dir= f"{new_dir}/tunnels/"
    for pcd in limit:
        # for each domain - compute its tunnels
        # create file and directories names
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        single_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
        os.mkdir(single_tunnels_dir)
        # compute tunnel start and end
        start = utils_PIN.predict_PIN_ligand_position(annotation_dir, pcd, 1)
        # compute tunnels
        compute_protein_tunnels(single_tunnels_dir, "mole_config/PIN/from_ligand.xml", pcd, cif_file, start=start, end=None)
        if not args.no_sort:
            sort_tunnels_PIN(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir, part="tunnels")
        # clean up directory sstructure for the next round 
        # and save it in te proper place


def compute_SLC2A3_protein_tunnels(limit, new_dir, annotation_dir, args):
    # PIN proteins have 4 start tunnel points and 5 end tunnel points
    # we have to create all possible tunnels
    tunnels_dir= f"{new_dir}/tunnels/"
    for pcd in limit:
        print("pcd", pcd)
        # for each domain - compute its tunnels
        # create file and directories names
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        single_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
        os.mkdir(single_tunnels_dir)
        # compute tunnel start and end
        start, end = utils_SLC2A3.predict_SLC2A3_tunnel_start_end(annotation_dir, pcd)
        print("start", start)
        print("end", end)
        # compute tunnels
        compute_protein_tunnels(single_tunnels_dir, "mole_config/SLC2A3/start_end.xml", pcd, cif_file, start=start, end=end)
#        compute_protein_tunnels(single_tunnels_dir, f"mole_config/SLC2A3/start.xml", pcd, cif_file, start=start, end=end)
#        if not args.no_sort:
        sort_tunnels_unknown_family(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir, part="all", all_ligands=True)
        # clean up directory sstructure for the next round 
        # and save it in te proper place


def compute_P02920_protein_tunnels(limit, new_dir, annotation_dir, args):
    # PIN proteins have 4 start tunnel points and 5 end tunnel points
    # we have to create all possible tunnels
    tunnels_dir= f"{new_dir}/tunnels/"
    for pcd in limit:
        print("pcd", pcd)
        # for each domain - compute its tunnels
        # create file and directories names
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        single_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
        os.mkdir(single_tunnels_dir)
        # compute tunnel start and end
        start, end = utils_P02920.predict_P02920_tunnel_start_end(annotation_dir, pcd)
        print("start", start)
        print("end", end)
        # compute tunnels
        compute_protein_tunnels(single_tunnels_dir, f"mole_config/SLC2A3/start_end.xml", pcd, cif_file, start=start, end=end)
#        compute_protein_tunnels(single_tunnels_dir, f"mole_config/SLC2A3/start.xml", pcd, cif_file, start=start, end=end)
#        if not args.no_sort:
        sort_tunnels_unknown_family(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir, part="all", all_ligands=True)
        # clean up directory sstructure for the next round 
        # and save it in te proper place


def filter_duplicate_tunnels(aux_tunnels_dir, tunnels_dir):
    if not os.path.exists(f"{tunnels_dir}/"):
        os.makedirs(f"{tunnels_dir}/")
    for fn_new in os.listdir(f"{aux_tunnels_dir}/"):
        if fn_new.startswith("tunnel"):
            for fn_old in os.listdir(f"{tunnels_dir}/"):
                if fn_new.startswith("tunnel"):
                    if utils_tunnel_pdb.test_are_tunnels_almost_same(f"{aux_tunnels_dir}/{fn_new}", f"{tunnels_dir}/{fn_old}", toleration = 2):
                        break
            shutil.copyfile(f"{aux_tunnels_dir}/{fn_new}", f"{tunnels_dir}/{fn_new}")
    return 0


def sort_tunnels(tunnels_dir, annotation_dir, pcd, log_dir, args, new_dir, output_file,  part="all"):
    if pcd[2] != "":
       domain = pcd[2]
    else:
       domain = '""'

    if part != "all":
        command=f"./sort_channels_based_on_position_between_SSEs_1.10.630.10.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {domain} --type {args.type} --part_only {part} "#>>{output_file}"
    else:
        command=f"./sort_channels_based_on_position_between_SSEs_1.10.630.10.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {domain} --type {args.type} "#>>{output_file}"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(new_dir), f"sort-fail-{time.time()}.txt")
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"sort tunnels failed, report save to {log_fn}")
        raise RuntimeError(f"sort tunnels failed for {pcd[0]}{pcd[1]}{pcd[2]}: see {log_fn}")
    return r


def sort_tunnels_1_20_1070_10(tunnels_dir, annotation_dir, pcd, log_dir, args, new_dir, output_file,  part="all"):
    if pcd[2] != "":
       domain = pcd[2]
    else:
       domain = '""'

    if part != "all":
        command=f"./sort_channels_based_on_position_between_SSEs_1.20.1070.10.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {domain} --type {args.type} --part_only {part}"
    else:
        command=f"./sort_channels_based_on_position_between_SSEs_1.20.1070.10.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {domain} --type {args.type}"
#    print("command", command)
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(new_dir), f"sort-fail-{time.time()}.txt")
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"sort tunnels failed, report save to {log_fn}")
        raise RuntimeError(f"sort tunnels failed for {pcd[0]}{pcd[1]}{pcd[2]}: see {log_fn}")
    return r


def sort_tunnels_unknown_family(tunnels_dir, annotation_dir, pcd, log_dir, args, new_dir, part="all", all_ligands=False):
    if pcd[2] != "":
       domain = pcd[2]
    else:
       domain = '""'

    if "part_only" in args and (args.part_only != None):
        command=f"./sort_channels_based_on_position_between_SSEs_unknown_family.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {domain} --type {args.type} --pt {args.part_only}"
    else:
        command=f"./sort_channels_based_on_position_between_SSEs_unknown_family.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {domain} --type {args.type} --pt {part}"
    if all_ligands:
        command = command + " --all_ligands"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(new_dir), f"sort-fail-{time.time()}.txt")
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"sort tunnels failed, report save to {log_fn}")
        raise RuntimeError(f"sort tunnels failed for {pcd[0]}{pcd[1]}{pcd[2]}: see {log_fn}")
    return r


def sort_tunnels_PIN(tunnels_dir, annotation_dir, pcd, log_dir, args, new_dir, part="all"):

    if part != "all":
        command=f"./sort_channels_based_on_position_between_SSEs_PIN.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {pcd[2]} --type {args.type} --part_only {part}"
    else:
        command=f"./sort_channels_based_on_position_between_SSEs_PIN.py --tunnels_dir {tunnels_dir} --annotation_dir {annotation_dir} --cif_dir {annotation_dir}  --protein {pcd[0]}  --output_fn pom  --chain {pcd[1]} --domain {pcd[2]} --type {args.type}"
#    print("command", command)
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(new_dir), f"sort-fail-{time.time()}.txt")
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"sort tunnels failed, report save to {log_fn}")
    return r


def do_molspecypic_output_1_10_630_10(new_dir):

    command=f"./molspecypic_input.py --dir {new_dir} --output {new_dir}"
    print("command", command)
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = f"{new_dir}/molspecypic_input-{time.time()}.txt"
        utils_common.error_report(new_dir, command, out, err, log_fn)
        logging.error(f"molspecypic output failed, report save to {log_fn}")
    return r


def output(tunnels_dir, annotation_dir, pcd, log_dir, args, output_file):

    command=f"./generate_pymol_command.py --tunnel_name oo --annotation_file {annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json  --prefix  /home/ivarekova/radka/bilkoviny/generate_family_entry/channels/channel_sorter/generated-tunnels-list_P08684-2025-09-09T12_07_41/cifs/ --tunnel 2b --width 1.6"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = os.path.join(os.path.abspath(log_dir), f"overprot-fail-{time.time()}.txt")
        utils_common.error_report(log_dir, command, out, err, log_fn)
        logging.error(f"overprot on 1.1.1.1 failed, report save to {log_fn}")
    return r


# read arguments
parser = argparse.ArgumentParser(
    description=
    '''for a HEM protein list generate set of tunnels - using combination of MOLE calls and sort the tunnels to standart tunnel types
       e.g.: ./generate_annotated_channels_for_protein.py --protein_list uniprot_mappig/pdb_uniprot_mapping-2026-01-07T19_27_15/list_W8VZW3  --type pdb --family 1.20.1070.10 
    ''')
parser.add_argument('--protein_list', required=False, action="store",
                    help="generate sorted protein list for given protein list")
parser.add_argument('--protein', required=False, action="store",
                    help="generate sorted protein list for given protein list")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result protein directory")
parser.add_argument('--type', required=True, action="store",
                    help="type of input database af(default) pdb(optional) mix(model have to be from pdb)")
parser.add_argument('--family', action='store', default="1.10.630.10",
                    help='print info output and leave debug data')
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('--rounds', action='store', default=6,
                    help='print info output')
parser.add_argument('--no_sort', action='store_true', default=False,
                    help='print info output')
parser.add_argument('--cif_dir', action='store',
                    help='set the directory containing protein cif files used for input')
parser.add_argument('--mole_config_file', action='store',
                    help='MOLE configuration file - used only for not known families (known 1.10.630.10, 1.20.1070.10)')
parser.add_argument('--part_only', action='store',
                    help='path/pore/tunnel/all')

args = parser.parse_args()


if (not args.protein_list and not args.protein):
    print("Either argument --protein_list or --protein have to be set")
    sys.exit(1)
if (args.protein_list and args.protein):
    print("Only one argument --protein_list or --protein have to be set")
    sys.exit(1)

if args.protein_list:
    name = os.path.basename(args.protein_list)
else:
    # just remove .cif suffix from name
    name = os.path.basename(args.protein)[:-4]

# ---------------------------------------------
# create new UniProt protein database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time

if args.output:
    new_dir = f"{args.output}/"
else:
    new_dir = f"generated-tunnels-{name}-{timestamp}"
os.mkdir(new_dir)

# ---------------------------------------------
# choose model domain - defined for processed families,
# the first item on list for new families
if args.family == "1.10.630.10":
    MODEL = utils_1_10_630_10.MODEL #["1akd","A", "00"]
elif args.family == "1.20.1070.10":
    MODEL = utils_1_20_1070_10.MODEL
elif args.family == "PIN":
    MODEL = utils_PIN.MODEL
else:
    # family was not pre-processed, no info about it
    # model is the first domain from domain list
    with (open(args.protein_list,"r")) as fpl:
        line = fpl.readline()
        if args.type == "pdb" or args.type == "mix":
            MODEL = [line[0:4], line[4:-3], line[-3:-1], line[0:4].lower()]

        else:
            MODEL = [line[:-1], "A", "00"]

# start logger
logger = start_logger(args.debug, False, new_dir, name)
last_lime = utils_common.time_log("0 Befor process starts", last_time)

annotation_dir = f"{new_dir}/annotation"
os.mkdir(annotation_dir)

# ---------------------------------------------
# now we have to download cifs for all domains

if args.protein_list:
    # args.protein_list was set
    # download domains list cifs from pdb database
    if args.type == "pdb" or args.type == "mix":
        pdb_pcds, not_in_pdb_dbs = read_input_pdb_cifs(args, annotation_dir)
    else:
        pdb_pcds, not_in_pdb_dbs = [], []

    # download domains list cifs from af dbs
    if args.type == "af" or args.type == "mix":
        af_pcds, not_in_af_dbs = read_input_af_cifs(args, annotation_dir, not_in_pdb_dbs)
    else:
        af_pcds, not_in_af_dbs = [], []

    if (args.type == "pdb" and not_in_pdb_dbs != []) or (args.type == "mixed" and not_in_af_dbs != []) or (args.type == "af" and not_in_af_dbs != []):
        # some protein from input protein list was not found - that is problem for only pdb list
        if args.type == "pdb" and not_in_pdb_dbs != []:
            print(f"Some protain was not found in pdb database ({not_in_pdb_dbs[0]})")
        if args.type == "mixed" and not_in_af_dbs != []:
            print(f"Some protain was not found in pdb and af database ({not_in_af_dbs[0]})")
        if args.type == "af" and not_in_af_dbs != []:
            print(f"Some protain was not found in af database ({not_in_af_dbs[0]})")
        sys.exit()
    pcds = pdb_pcds + af_pcds
else:
    # args.protein was set
    pcds = read_input_cif_from_file(args.protein, annotation_dir)

print 

# download model - if it is not already downloaded
if [MODEL[0], MODEL[1], MODEL[2], MODEL[0]] not in pcds:
    if args.type == "af":
        # model should be in af dbs
        copy_or_download_cif_from_af_to_dir(args, annotation_dir, protein_name=MODEL[0])
    elif args.type == "mix" or args.type == "pdb":
        # model should be in pcd dbs
        copy_or_download_cif_from_pdb_to_dir(args, annotation_dir, protein_name=MODEL[0])
    set_contains_model = False
else:
    set_contains_model = True
pcds.append([MODEL[0], MODEL[1], MODEL[2], MODEL[0]])

last_lime = utils_common.time_log(f"1 protein cif file downloaded {len(pcds)}proteins", last_time)
print("===================================")
print(f"1/ protein cif file downloaded {len(pcds)} domains ({pcds})")

# ---------------------------------------------
# remove ligands - leave only HEMs - not to bother with them
for fn in os.listdir(annotation_dir):
    if args.family == "1.10.630.10":
        # leave HEM ligand - target of tunnels
        cif_file = f"{annotation_dir}/{fn}"
        aux_cif_file = f"{annotation_dir}/{fn[:-4]}_without_ligands.cif"
        utils_read_cif.remove_ligands(cif_file, aux_cif_file, ["HEM", "HEC", "HEB"])
        pass
    elif args.family == "1.20.1070.10":
        # no need to remove ligands
        pass
    elif args.family == "SLC2A3":
        pass
    else:
        # remove all ligands - no need to leave it
        utils_read_cif.remove_ligands(f"{annotation_dir}/{fn}", f"{annotation_dir}/{fn}.aux", [])
        os.remove(f"{annotation_dir}/{fn}")
        shutil.move(f"{annotation_dir}/{fn}.aux", f"{annotation_dir}/{fn}")

# ---------------------------------------------
# annotate all domains
logging.info("INFO: annotation of input cifs ")
overprot_dir = f"{new_dir}/overprot"
os.mkdir(overprot_dir)

# create overprot input json
overprot_file = f"{overprot_dir}/overprot_input.json"
if args.family == "1.10.630.10":
    # find only domains from "1.10.630.10"
    # remove all domains which are shorter then 300 res - they could not be full
    pcds2 = []
    for pcd in pcds:
        res_list = utils_read_cif.find_auth_asym_id_len(f"{annotation_dir}/{pcd[0]}.cif", pcd[1])
        if res_list > 350 and res_list<600:
            # the size of domain have to be appropriate
            coords = utils_read_cif.find_fe_hem_coords_auth(f"{annotation_dir}/{pcd[0]}.cif", pcd[1])
            if coords == 0:
                print("no HEM - ommit")
            if args.type != "pdb" or coords != 0:
                # either pcd can be from af (do not have ligands) or ocd have HEM
                pcds2.append(pcd)
    pcds = pcds2
    if set_contains_model:
        create_overprote_input(overprot_file, pcds[:-1], args.type, annotation_dir, "1.10.630.10")
    else:
        create_overprote_input(overprot_file, pcds, args.type, annotation_dir, "1.10.630.10")
elif args.family == "1.20.1070.10":
    # find only domains from "1.20.070.10"
    if set_contains_model:
        create_overprote_input(overprot_file, pcds[:-1], args.type, annotation_dir, "1.20.1070.10")
    else:
        create_overprote_input(overprot_file, pcds, args.type, annotation_dir, "1.20.1070.10")
else:
    # we do not know the family id thus we find all domains
    print("set_contains_model", set_contains_model)
    if set_contains_model:
        create_overprote_input(overprot_file, pcds[:-1], args.type, annotation_dir)
    else:
        create_overprote_input(overprot_file, pcds, args.type, annotation_dir)

# do annotation
ret = set_ubertemplate(overprot_dir, overprot_file, new_dir, annotation_dir, new_dir, fast=False)
if ret != 0:
    # if there was a problem, just finish
    print(f"ERROR: ubertemplate for {args.protein_list} can't be computed")
    sys.exit()
last_lime = utils_common.time_log("3 ubertemplate computed", last_time)
print("===================================")
print("3 ubertemplate computed")
# TODO: annotation is in chain: auth_asym_id, res auth_seq_id

# move wated files to right directories
for pcd in pcds:
    shutil.copyfile(f"{overprot_dir}/cif_cealign/{pcd[0]}{pcd[1]}{pcd[2]}-clust.sses.json", f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json")
    shutil.copyfile(f"{overprot_dir}/cif_cealign/{pcd[0]}{pcd[1]}{pcd[2]}-rotation.json", f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-rotation.json")
shutil.copyfile(f"{overprot_dir}/results/consensus.cif", f"{annotation_dir}/consensus.cif")
shutil.copyfile(f"{overprot_dir}/results/consensus.sses.json", f"{annotation_dir}/consensus-template.sses.json")

print("===================================")
print("4 annotation and rotation copied")
# ---------------------------------------------
# compute tunnels
# create directory structure - directory for tunnel generation
tunnels_dir= f"{new_dir}/tunnels"
os.mkdir(tunnels_dir)
if args.family == "1.10.630.10":
    # 1.10.630.10 family needs several rounds of
    #  MOLE callculations - thus directory structure have to be bigger
    # aux_tunnels dir is necessary only with several rounds of mole run
    aux_tunnels_dir= f"{new_dir}/tunnels_aux"
    os.mkdir(aux_tunnels_dir)
    # water channel part is 1.10.630.10 speciality
    tunnels_dir_w= f"{new_dir}/tunnels_w"
    os.mkdir(tunnels_dir_w)

# if model is not in our domain list just do not compute its tunnels
if set_contains_model:
    # all pcds should be showed - there is no artificial model
    limit = pcds[:-1]
else:
    # the last instance in pcd is artificial model
    limit = pcds[:-1]

# compute tunnels for all domains on list (model is already removed if necessary)
print("5 =================================================")
if args.family == "1.10.630.10":
    for pcd in limit:
        print(f"=== pcd {pcd[0]}{pcd[1]}{pcd[2]}  (auth_asym_id), family 1.10.630.10) ----- ")
        prot_annotation = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
        rotation_json = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-rotation.json"
        model_annotation = f"{annotation_dir}/{pcds[-1][0]}{pcds[-1][1]}{pcds[-1][2]}-annotated.sses.json"
        model_rotation = f"{annotation_dir}/{pcds[-1][0]}{pcds[-1][1]}{pcds[-1][2]}-rotation.json"
        # at first we have to test weather domain is full, sometimes domains are uncomplete
        # in this cases tunnel generation works bogus and result tunnels are unusefull
        domain_is_ok = is_domain_complete_1_10_630_10(model_annotation, prot_annotation, percentage=75)
        if not domain_is_ok:
            # we have to remove rotation file of this domain,
            # domain is not complete, so we don't want tu use it
            os.remove(f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-rotation.json")
            continue

        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        model_cif_file = f"{annotation_dir}/{pcds[-1][0]}.cif"
        aux_cif_file = f"{annotation_dir}/{pcd[0]}_without_ligands.cif"
        if args.type and args.type == "pdb":
            # pdb proteins contain ligands - no need to predict HEM positions
            auth_chain = utils_read_cif.transform_label_asym_id_to_auth_asym_id(cif_file, pcd[1])
            print(auth_chain, "<->", pcd[1])
#            if auth_chain != pcd[1]:
#               print("!!!", auth_chain, "<->", pcd[1])
            start = utils_read_cif.find_fe_hem_coords_auth(cif_file, auth_chain)
        elif args.type and args.type == "mix":
            auth_chain = utils_read_cif.transform_label_asym_id_to_auth_asym_id(cif_file, pcd[1])
            print(auth_chain, "<->", pcd[1])
#            if auth_chain != pcd[1]:
#               print("!!!", auth_chain, "<->", pcd[1])
            start = utils_read_cif.find_fe_hem_coords_auth(cif_file, auth_chain)
            if start == 0:
                # in case there is not a SSEs which is necessary to compute hem,
                # do not fall, only ommit this pcd
                start = utils_1_10_630_10.predict_hem_position(model_annotation, prot_annotation, f"{annotation_dir}/{pcd[0]}.cif", rotation_json, model_rotation, model_cif_file)
        else:
            # af proteins do not contain ligands - we have to predict HEM positions
            try:
                # in case there is not a SSEs which is necessary to compute hem,
                # do not fall, only ommit this pcd
                start = utils_1_10_630_10.predict_hem_position(model_annotation, prot_annotation)
                if start == None:
                    continue
            except:
                continue

        compute_1_10_630_10_protein_tunnels_in_several_rounds(new_dir, pcd, aux_cif_file, start, aux_tunnels_dir, tunnels_dir_w, args)
        last_lime = utils_common.time_log("4 tunnels (in family 1.10.630.10) computed", last_time)

        # TODO - remove duplicates
        filter_duplicate_tunnels(f"{aux_tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}", f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}")
        solvent_tunnels_dir= f"{new_dir}/tunnels/{pcd[0]}{pcd[1]}{pcd[2]}"
        water_tunnels_dir= f"{new_dir}/tunnels_w/{pcd[0]}{pcd[1]}{pcd[2]}"
        if not args.no_sort:
            sort_tunnels(solvent_tunnels_dir, annotation_dir, pcd, new_dir, args, new_dir, "/dev/null", part="solvent")
            sort_tunnels(water_tunnels_dir, annotation_dir, pcd, new_dir, args, new_dir, "/dev/null", part="water")
        last_lime = utils_common.time_log("5 tunnels are sorted (in family 1.10.630.10) computed", last_time)
        print("    tunnels sorted")

elif args.family == "1.20.1070.10":
    for pcd in limit:
        print(f"pcd {pcd[0]}{pcd[1]}{pcd[2]}  (auth_asym_id), family 1.20.1070.10)")
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"

        prot_annotation = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
        model_annotation = f"{annotation_dir}/{pcds[-1][0]}{pcds[-1][1]}{pcds[-1][2]}-annotated.sses.json"
        # at first we have to test weather domain is full, sometimes domains are not complete or does not belong to the relevant family
        # in this cases tunnel generation works bogus and result tunnels are unusefull
        domain_is_ok = is_domain_complete_1_20_1070_10(model_annotation, prot_annotation, percentage=80)
        if not domain_is_ok:
            print(f"    domain is not complete - tunnels are not computed")
            # we have to remove rotation file of this domain,
            # domain is not complete, so we don't want tu use it
            os.remove(f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-rotation.json")
            continue

        # GPCR family - membrane protein - 7 helices parallel channels should go inside
        tunnels_dir= f"{new_dir}/tunnels/"
        single_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
        os.mkdir(single_tunnels_dir)

        # have to be computed start and end
        l_id = utils_read_cif.transform_auth_asym_id_to_label_asym_id(cif_file, pcd[1])
        start, end = utils_1_20_1070_10.predict_1_20_1070_10_tunnel_start_end(annotation_dir, [pcd[0], pcd[1], pcd[2]])
        middle = utils_1_20_1070_10.predict_1_20_1070_10_ligand_position(annotation_dir, [pcd[0], pcd[1], pcd[2]])
        # TODO: input should be protein, chain, domain in ......
        if args.mole_config_file:
            mole_file = args.mole_config_file
            if mole_file.endswith("extracellular.xml"):
                ec_point = utils_1_20_1070_10.predict_1_20_1070_10_ec_binding_point_from_ligand(annotation_dir, [pcd[0], pcd[1], pcd[2]])
                if not ec_point:
                    ec_point = utils_1_20_1070_10.predict_1_20_1070_10_ec_binding_point_from_ligand_p07700(annotation_dir, [pcd[0], pcd[1], pcd[2]])
                compute_protein_tunnels(single_tunnels_dir, args.mole_config_file, pcd, cif_file, start=ec_point)
            else:
                compute_protein_tunnels(single_tunnels_dir, args.mole_config_file, pcd, cif_file, start=start, end=end, middle=middle)
        else:
            compute_protein_tunnels(single_tunnels_dir, f"mole_config/1.20.1070.10/whole_path.xml", pcd, cif_file, start=start, end=end, middle=middle)
        # TODO: annotation is in chain: auth_asym_id, res auth_seq_id

        last_lime = utils_common.time_log("4 tunnels (in family 1.20.1070.10) computed", last_time)
        print("    tunnels computed")

        if not args.no_sort:
            if args.mole_config_file:
                if args.part_only:
                    sort_tunnels_1_20_1070_10(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir, "solvent", part=args.part_only)
                else:
                    sort_tunnels_1_20_1070_10(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir, "solvent", part="all")
            else:
                sort_tunnels_1_20_1070_10(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir, "solvent")
        last_lime = utils_common.time_log("5 tunnels are sorted (in family 1.20.1070.10) computed", last_time)
        print("    tunnels sorted computed")

elif args.family == "PIN":
    compute_PIN_protein_tunnels_in_several_rounds(limit, new_dir, annotation_dir, args)
    print("5 tunnels are sorted (in family PIN) computed")
elif args.family == "PIN_from_ligand":
    compute_PIN_protein_tunnels_from_ligand(limit, new_dir, annotation_dir, args)
    print("5 tunnels are sorted (in family PIN) computed")
elif args.family == "SLC2A3":
    print("5 tunnels are sorted (in family SLC2A3) computed", args)
    compute_SLC2A3_protein_tunnels(limit, new_dir, annotation_dir, args)
    print("5 tunnels are sorted (in family SLC2A3) computed")
elif args.family == "P02920":
    print("5 tunnels are sorted (in family P02920) computed", args)
    compute_P02920_protein_tunnels(limit, new_dir, annotation_dir, args)
    print("5 tunnels are sorted (in family P02920) computed")


else:
    for pcd in limit:
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        print("UNKNOWN")

        # unknown family - we know nothing about it
        tunnels_dir= f"{new_dir}/tunnels/"
        single_tunnels_dir= f"{tunnels_dir}/{pcd[0]}{pcd[1]}{pcd[2]}"
        os.mkdir(single_tunnels_dir)
        compute_protein_tunnels(single_tunnels_dir, args.mole_config_file, pcd, cif_file)
        last_lime = utils_common.time_log("4 tunnels (in unknown family) computed", last_time)
        print("4 tunnels (in unknown family PIN) computed")

        if not args.no_sort:
            sort_tunnels_unknown_family(f"{single_tunnels_dir}/pdb/profile", annotation_dir, pcd, new_dir, args, new_dir)
        last_lime = utils_common.time_log("5 tunnels are sorted (in unknown family) computed", last_time)
        print("5 tunnels are sorted (in unknown family) computed")

# create aggregation_files per tunnel type
tunel_stats = {}
for tdir in  os.listdir(new_dir):
    if os.path.isdir(f"{new_dir}/{tdir}") and ((tdir != "annotation") and (tdir != "cifs") and (tdir != "overprot") and (not tdir.startswith("tunnel"))):
        tunel_stats[tdir]={}
        counter = set()
        all_tunnel_name=f"tunnel_{tdir}.pdb"
        with open(f"{new_dir}/{tdir}/{all_tunnel_name}", "w") as fw:
            for file_name in os.listdir(f"{new_dir}/{tdir}"):
                if file_name != all_tunnel_name:
                    fnsplit=file_name.split("_")
                    if fnsplit[1] not in tunel_stats[tdir]:
                        tunel_stats[tdir][fnsplit[1]] = set()

                    #print("fff",fnsplit, counter)
                    counter.add(fnsplit[1])
                    with open(f"{new_dir}/{tdir}/{file_name}", "r") as fr:
                        line = fr.readline()
                        while line:
                            fw.write(line)
                            line = fr.readline()
                    rad = utils_tunnel_pdb.get_tunnel_min_radius(f"{new_dir}/{tdir}/{file_name}", ommit_start_part=5, ommit_end_part=0)
                    print(f"--- {new_dir}/{tdir}/{file_name} - radius: {rad}")
                    tunel_stats[tdir][fnsplit[1]].add(rad)

number_of_domains = 0
for cif_file in os.listdir(f"{new_dir}/cifs"):
    number_of_domains += 1

for tt in tunel_stats:
    print(f"{tt} - {len(tunel_stats[tt])}/{number_of_domains}={int(len(tunel_stats[tt])*100/number_of_domains)}")
    for fn in tunel_stats[tt]:
        print(f"TUNNEL {tt} {fn} - {max(tunel_stats[tt][fn])} ({tunel_stats[tt][fn]})")

if (args.family == "1.10.630.10") and (not set_contains_model):
    os.remove(f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json")
    os.remove(f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-rotation.json")

if (args.family == "1.10.630.10"):
    # do molespecific output
    utils_molespecypic_output.do_molespecypic_output(new_dir,f"{new_dir}/{name}.json","-")

    if not args.debug:
        # remove useless directories
        shutil.rmtree(annotation_dir)
        shutil.rmtree(overprot_dir)
        shutil.rmtree(tunnels_dir)
        shutil.rmtree(f"{new_dir}/tunnels1")
        shutil.rmtree(f"{new_dir}/tunnels2")
        shutil.rmtree(f"{new_dir}/tunnels3")
        shutil.rmtree(f"{new_dir}/tunnels4")
        shutil.rmtree(f"{new_dir}/tunnels5")
        shutil.rmtree(f"{new_dir}/tunnels6")
        shutil.rmtree(f"{new_dir}/tunnels_aux")
        shutil.rmtree(f"{new_dir}/tunnels_w")
        shutil.rmtree(f"{new_dir}/cifs")
        shutil.rmtree(f"{new_dir}/garbage")
