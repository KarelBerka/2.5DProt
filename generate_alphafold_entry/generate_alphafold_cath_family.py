#!/usr/bin/env python3
# tool which generate for given alphafold family
# its diagrams

import logging
import sys
import argparse
import time
import os
import subprocess
import shutil
import json
import numpy

import utils_common
import utils_annotation
import utils_alphafold_confidence
import utils_synchronize_rotation_files
import utils_pdbe_api
import color_me

sys.path.append('..')
import utils_links
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


def domain_already_used2(protein, dom, data):
    # find out whether the domain with this name is already in data dbs
    for rec in data:
        if rec[0] == protein and rec[2] == dom:
            return True
    return False


def find_family_domains(family, new_dir, pdb_dir, cath_file, verbose):
    # find list of all domains in the given family
    # actually just read cath/cath_family_domains file
    af_list = []
    if verbose:
        print("find family domains")

    with open(cath_file, "r") as fd:
        line = fd.readline()
        counter = 0
        while line:

            counter += 1
            if (counter % 1000 == 0) and (verbose):
                print(f"* {counter}")

            # example of line:
            # af_A8MT69_40_81	cath_exp	42	homo_sapiens	6	6.10.130.30	CATH_HMM_S95_SSAP	98.25	FALSE	1	23.81	G
            line_split = line.split("\t")
            if line_split[5] == family:
                # if the line belonges to the given family - read the domain and add it to list
                full_name = line_split[0][3:]  # just omit one whitespace and "af_" prefix
                full_name_split = full_name.split("_")
                dname = full_name_split[0]
                drange = full_name_split[1:]
                fam = line_split[3]
                if (fam == "homo_sapiens") and (int(drange[1]) > 1400) and (os.path.isfile(f"{pdb_dir}/{dname[0]}/{dname[1]}/{dname[2]}/{dname[3]}/AF-{dname}-F2-model_v4.pdb")):
                    # if the domain is on "homo_sapiens" organism, there have to be found -F{number}
                    # version of pdb file, if these files are present (-F2 is present test) and pdb is split
                    # set version tag and range tag properly
                    ind = ((int(drange[1]) - 1200) // 200) + 1
                    fversion = f"F{ind}"
                    drange2 = []
                    for i in drange:
                        drange2.append(int(i)-(ind-1)*200)
                    drange = drange2
                else:
                    # default version tag is 1 - it counld not be 1 only in "homo_sapiens"
                    fversion = "F1"

                # domain tag - each af protein can be in several domains - thus if the protein
                # was used in this family, just increase this "domain" number
                dom = 0
                while domain_already_used2(f"{dname}-{fversion}", dom, af_list):
                    dom += 1

                af_list.append([f"{dname}-{fversion}", drange, dom, fam, fversion])
            line = fd.readline()
    return af_list


def simlink_cath_alphafold_list(af_list, pdb_dir, output_dir, pdb_dir2, new_dir):
    # af_location="https://alphafold.ebi.ac.uk/files"
    # cath location https://liveuclac-my.sharepoint.com/personal/ucbtnb4_ucl_ac_uk/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fucbtnb4%5Fucl%5Fac%5Fuk%2FDocuments%2FCATH%2DAlphaFold%2D21organisms&ga=1
    # alphafold list could be only simlinked, not to have several copies
    # https://alphafold.ebi.ac.uk/files/AF-A0A0A1EI90-F1-model_v4.pdb
    # pdb_dir structure:
    #         {pdb_dir}/0/2/Z/AF-A0A7W8CZ20-F1-model_v4_TED01.{pdb, cif}
    # pdb_dir2 expected structure:
    #         {pdb_dir2}/AF-A0A7W8CZ20-F1-model_v4_TED01.{pdb, cif}

    af_addr = "https://alphafold.ebi.ac.uk/files"
    af_list2 = []

    for af in af_list:

        # at first check pdb file then the same for cif file
        res = True
        for suffix in ['pdb', 'cif']:
            # at first check pdb/cif file presence
            if (not os.path.isfile(f"{output_dir}/{af[0]}.{suffix}")) and (not os.path.islink(f"{output_dir}/{af[0]}.{suffix}")):
                # if it is not in output dir we have to get it

                link_pdb = f"{af_addr}/AF-{af[0]}-model_v4.{suffix}"
                pdb1_filename = f"{af[0][-4]}/{af[0][-5]}/{af[0][-6]}/AF-{af[0]}-model_v4.{suffix}"
                pdb2_filename = f"{af[0]}.{suffix}"
                # at first try to browse pdb_dir directory
                if (pdb_dir is not None) and (os.path.isfile(f"{pdb_dir}/{pdb1_filename}") or os.path.islink(f"{pdb_dir}/{pdb1_filename}")):
                    try:
                        os.symlink(f"../../{pdb_dir}/{pdb1_filename}", f"{output_dir}/{af[0]}.{suffix}")
                        logging.info(f"{af[0]}.{suffix} symlinked from {pdb_dir}")
                    except:
                        # if this fails just grab the pdb/cif file from web
                        ret = utils_common.curl_file(link_pdb, new_dir, f"{output_dir}/{af[0]}.{suffix}", err_report=False)
                        if not ret:
                            # files not found - the CATH family <-> domain mapping could not be up-to-date
                            # write an error, but continue to process the rest of family
                            logging.error(f"Alpha fold entry {af} is not found on web")
                            res = False
                            break
                        else:
                            logging.info(f"{af[0]}.{suffix} downloaded")

                # try to browse pdb2_dir directory as well
                elif (pdb_dir2 is not None) and (os.path.isfile(f"{pdb_dir2}/{pdb2_filename}") or os.path.islink(f"{pdb_dir2}/{pdb2_filename}")):
                    try:
                        os.symlink(f"../../{pdb_dir2}/{pdb2_filename}", f"{output_dir}/{af[0]}.{suffix}")
                        logging.info(f"{af[0]}.{suffix} symlinked from {pdb_dir2}")
                    except:
                        # if this fails just grab the pdb/cif file from web
                        ret = utils_common.curl_file(link_pdb, new_dir, f"{output_dir}/{af[0]}.{suffix}", err_report=False)
                        if not ret:
                            # files not found - the CATH family <-> domain mapping could not be up-to-date
                            # write an error, but continue to process the rest of family
                            logging.error(f"Alpha fold entry {af} is not found on web")
                            res = False
                            break
                        else:
                            logging.info(f"{af[0]}.{suffix} downloaded")

                # we have to download pdb/cif from net
                else:
                    ret = utils_common.curl_file(link_pdb, new_dir, f"{output_dir}/{af[0]}.{suffix}", err_report=False)
                    if not ret:
                        # files not found - the CATH family <-> domain mapping could not be up-to-date
                        # write an error, but continue to process the rest of family
                        logging.error(f"Alpha fold entry {af} is not found on web")
                        res = False
                        break
                    else:
                        logging.info(f"{af[0]}.{suffix} downloaded")
        if res:
            # both files pdb and ciff are downloaded correctly
            af_list2.append(af)

    return af_list2


def simlink_cath_alphafold_list_annoation_files(af_list, src_ann_dir, annotation_dir, new_dir):
    # we have to copy files:
    # ../overprot/results/consensus.sses.json
    # {src_ann_dir}/{af[0][:-3]}-F1A0{012}-annotated.sses.json
    # TODO: have to be fixed based on new cath AF dbs format
    os.mkdir(f"{new_dir}/overprot")
    os.mkdir(f"{new_dir}/overprot/results")
    os.symlink(f"../../../{src_ann_dir}/../overprot/results/consensus.sses.json", f"{new_dir}/overprot/results/consensus.sses.json")

    for af in af_list:
        ann_filename = f"{af[0][:-3]}-F1A00-annotated.sses.json"
        if os.path.isfile(f"{src_ann_dir}/{ann_filename}"):
            os.symlink(f"../../{src_ann_dir}/{ann_filename}", f"{annotation_dir}/{ann_filename}")
        ann_filename = f"{af[0][:-3]}-F1A01-annotated.sses.json"
        if os.path.isfile(f"{src_ann_dir}/{ann_filename}"):
            os.symlink(f"../../{src_ann_dir}/{ann_filename}", f"{annotation_dir}/{ann_filename}")
        ann_filename = f"{af[0][:-3]}-F1A02-annotated.sses.json"
        if os.path.isfile(f"{src_ann_dir}/{ann_filename}"):
            os.symlink(f"../../{src_ann_dir}/{ann_filename}", f"{annotation_dir}/{ann_filename}")
    return True


def domain_already_used(domain, chain, dom, data):
    # find out whether the domain with this name is already in data dbs
    domain_name = str(domain[0])+str(chain)+str(dom).zfill(2)
    for rec in data:
        if rec["domain"] == domain_name:
            return True
    return False


def add_channel_file_to_ch_dir(annotation_dir, domain, new_dir, ch_source_dir=False):
    # copy channel file from input directory (ch_source_dir) if it is set
    # or from utils_links.link_channels_af link otherwise
    if ch_source_dir:
        sd0 = domain[0].split("-")
        ch_fn_out = f'{ch_source_dir}/{sd0[0]}_result.json'
        ch_fn_in = f"{annotation_dir}/channel-{domain[0][:-3]}.json"
        if os.path.isfile(f"{ch_fn_out}"):
            shutil.copy(ch_fn_out, ch_fn_in)
    else:
        link_chann = utils_links.link_channels_af(domain[0][:-3])
        ret = utils_common.curl_file(link_chann, new_dir, f"{annotation_dir}/channel-{domain[0][:-3]}.json", err_report=False)
        return ret


def add_2dprot_layout(family, directory, cif_dir):
    # copy one layout from https://2dprots.ncbr.muni.cz/
    # to sync diagrams with 2dprots ones

    # we have to find which files are in
    # 2dprots dbs
    # use /files/family/<family>/domain_list.txt
    rec = None

    link = f"https://2dprots.ncbr.muni.cz/files/family/{family}/domain_list.txt"
    twodprots_dir = f"{directory}/2dprots/"
    os.mkdir(twodprots_dir)
    domain_list_file = f"{twodprots_dir}/domain_list.txt"
    ret = utils_common.curl_file(link, directory, domain_list_file)
    logging.debug(f"2dprots version of domain_list.txt file exists: {ret}")
    if ret is False:
        logging.debug("no domain_list.txt - no diagrams to synchronize with, create diagrams from scratch")
    else:
        # we found domain_list.txt try to download layout based on it
        with open(domain_list_file, "r") as fd:
            line = fd.readline()
            while line:
                # line contains :('2cw7','A', '02')
                sline = line.split("'")
                domain = f"{sline[1]}{sline[3]}{sline[5]}"
                link = f"https://2dprots.ncbr.muni.cz/files/domain/{domain}/layout.json"
                domain_file = f"layout-{domain}.json"
                ret = utils_common.curl_file(link, directory, f"{twodprots_dir}/{domain_file}")
                logging.debug(f"2dprots version of {domain} file exists: {ret}")
                if not ret:
                    # this layout was not found
                    line = fd.readline()
                    continue
                # now we have to download its pdb
                utils_pdbe_api.download_protein_cif(cif_dir, sline[1])
                # we have to found the domain ranges as well
                ranges = utils_pdbe_api.get_domain_ranges(sline[1], sline[3], sline[5])
                # format [(1,2),(3,4)] - 3dpuA01
                rec = [sline[1], sline[3], sline[5], "2dprots", ranges]
                break
    return rec


def transform_list_to_db(dblist):
    # transform list [1,2,3,4] of ranges to the proper 1:15,45:90,... format
    ret = ""
    start = 0
    while len(dblist) >= start + 2:
        if ret == "":
            ret = f"{dblist[start]}:{dblist[start+1]}"
        else:
            ret += f",{dblist[start]}:{dblist[start+1]}"
        start += 2
    return ret


def transform_range_list_to_db(dblist):
    # transform list [(1,2),(3,4)] of ranges to the proper 1:15,45:90,... format
    ret = ""
    start = 0
    for rec in dblist:
        if ret == "":
            ret = f"{rec[0]}:{rec[1]}"
        else:
            ret += f",{rec[0]}:{rec[1]}"
        start += 1
    return ret


def append_to_overprot_input_file(overprot_file, cif_dir, domain, chain="A"):
    # create overprot input file
    with open(overprot_file, "r") as fd:
        data = json.load(fd)
    rec = {}
    if domain[3] != "2dprots":
        # alphafold domain record
        rec["domain"] = str(domain[0]) + str(chain) + str(domain[2]).zfill(2)
        rec["pdb"] = str(domain[0])
        rec["chain_id"] = chain
        rec["ranges"] = transform_list_to_db(domain[1])
    else:
       # 2dprots domain record added to sync rotation with 2dprots dbs
        rec["domain"] = domain[0] + domain[1] + domain[2]
        rec["pdb"] = domain[0]
        sa_id = utils_pdbe_api.transform_cath_chain_id_to_struct_asym_id(domain[0], domain[1])
        rec["chain_id"] = sa_id
        rec["ranges"] = transform_range_list_to_db(domain[4])
    rec["auth_chain_id"] = rec["chain_id"]
    rec["auth_ranges"] = rec["ranges"]
    data.append(rec)
    with (open(overprot_file, "w")) as fd:
        json.dump(data, fd, indent=4, separators=(",", ": "))


def set_ubertemplate(out_dir, d_file, log_dir, cif_dir, fast=False):
    # use overprot tool to download set ubertemplate for whole family
    # we have to call overprot for its directory,
    # and read activate file before we call overprot itself
    # fast set if we want to speed up whole process
    cwd = os.getcwd()
    os.chdir("../overprot/OverProtCore")

    activ_co = "source venv/bin/activate"
    if fast:
        overprot_co = f"python3 overprot.py 1.1.1.1 {cwd}/{out_dir} --sample_size 200 --domains {cwd}/{d_file} --structure_source file://{cwd}/{cif_dir}/"+"{pdb}.cif"
    else:
        overprot_co = f"python3 overprot.py 1.1.1.1 {cwd}/{out_dir} --sample_size 1000 --domains {cwd}/{d_file} --structure_source file://{cwd}/{cif_dir}/"+"{pdb}.cif"
    command = f"bash -c '{activ_co}; {overprot_co}'"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = f"{log_dir}/overprot-fail-{time.time()}.txt"
        log_fn = f"overprot-fail-{time.time()}.txt"
        utils_common.error_report(log_dir, command, out, err, log_fn)
        logging.error(f"overprot on 1.1.1.1 failed, report save to {log_fn}")

    os.chdir(cwd)
    return r


def rotate_channel(channel_filename, channel_rotation, aux_channel_filename):
    # rotate channel files to the same rotation
    if not (os.path.isfile(channel_filename)):
        return
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


def generate_cl(cif_dir, af_list, chain, working_dir, annotation_dir, layout_dir, start_layout_filename, args):
    counter = 0
    for pdi in af_list:
        counter += 1
        if (counter % 10000) and (args.verbose):
            print(counter, "generate_cl pcs completed", pdi)
        if pdi[3] == "2dprots":
            protein, domain, chain, aux, i = pdi
            interval = chain
            layout_file = f"{layout_dir}/layout-{pdi[0]}{pdi[1]}{pdi[2]}.json"
            dom = domain
        else:
            protein, domain, interval, org, f = pdi
            layout_file = f"{layout_dir}/layout-{protein}-{domain[0]}-{domain[1]}.json"
            dom = "A"
        annotation_file = f"{annotation_dir}/{protein}{dom}{str(interval).zfill(2)}-annotated.sses.json"

        if not args.channels_tmp or not os.path.isfile(args.channel_filename):
            command = f"./generate_domain_layout_from_scratch.py {annotation_file} {chain} {layout_file} {protein} alpha_fold {start_layout_filename}"
        else:
            if len(domain) > 1:
                channel_filename = f"{annotation_dir}/channel-{protein[:-3]}.json"
                aux_channel_filename = f"{annotation_dir}/{protein}{chain}{domain[0]}_{domain[1]}-channel.json"
            else:
                channel_filename = f"{annotation_dir}/channel-{protein}.json"
                aux_channel_filename = f"{annotation_dir}/{protein}{chain}{domain}-channel.json"
            rotation = f"{working_dir}/overprot/cif_cealign/{protein}{dom}{str(interval).zfill(2)}-rotation.json"
            rotate_channel(channel_filename, rotation, aux_channel_filename)
            command = f"./generate_domain_layout_from_scratch.py {annotation_file} {chain} {layout_file} {protein} alpha_fold {start_layout_filename} {aux_channel_filename}"

        logging.debug(f"DEBUG: Generating layout for {protein}")
        logging.info(f"INFO: Command {command}")
        if args.layout_dir:
            old_layout_file = f"{args.layout_dir}/layout-{protein}-{domain[0]}-{domain[1]}.json"
            if os.path.isfile(old_layout_file):
                os.symlink(f"../{old_layout_file}", f"layout-{protein}-{domain[0]}-{domain[1]}.json")
                continue
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()
        if (r != 0) or (not os.path.isfile(layout_file)):
            log_fn = f"{working_dir}/generate_protein_layout_fs-fail.txt"
            utils_common.error_report(working_dir, command, None, err, log_fn)
            logging.error(f"generate_domain_layout_from_scratch.py failed - report save to {log_fn}")
    return True


def multiple_diagrams(detection_dir, layout_dir, first_svg, not_ligands):
    svg_args = argparse.Namespace(
        only="template",
        add_descriptions="multi",
        layouts=[f"{layout_dir}"],
        ligands=False,
        opac='1',
        debug=False,
        irotation=f"{layout_dir}/rotation.json",
        itemplate=f"{layout_dir}/{args.family}.json",
        iligands=None,
        orotation=None,
        otemplate=None,
        output=f"{layout_dir}/{args.family}_template.svg",
        oligands=None,
        data=None,
        groups=None,
        fill=False
    )
    if not_ligands:
        svg_args.ligands = False
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")
    svg_args.only = "sses"
    svg_args.output = f"{layout_dir}/{args.family}_sses2.svg"
    svg_args.opac = "0"
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{detection_dir}")


def averager(rlist):
    s = 0
    for i in rlist:
        s += i
    return (float(s))


def create_layout_am_color(layout, am_file, new_layout):
    # TODO: comment
    with open(am_file, "rb") as fd:
        am_res = {}
        line = fd.readline()
        while line:
            strline = str(line)[2:-2]
            ssline = strline.split(" ")
            sline = []
            for ss in ssline:
                if ss != "":
                    sline.append(ss) 
            # just skip the intro rows
            if len(sline) > 2 and sline[0]=="ATOM":
                res = int(sline[5])
                if res not in am_res:
                    am_res[res] = set()
                am_res[res].add(float(sline[10]))
            line = fd.readline()

    with open(layout, "r") as fd:
        data = json.load(fd)
    am = {}

    for sses in data["helices_residues"]:
        am[sses] = {}
        res_int = data["helices_residues"][sses]
        r = res_int[1]
        while r < res_int[0]:
            av = sum(am_res[r]) / len(am_res[r])
            am[sses][r] = int(av*100)
            r += 1

    data["am"] = am
    with open(new_layout, "w") as fd:
        json.dump(data, fd, indent=4, separators=(',', ':'))


def create_layout_sift_color(layout, sift_file, new_layout):
    # TODO: comment
    with open(sift_file, "r") as fd:
        sift_res = {}
        line = fd.readline()
        while line:
            strline = str(line)
            sline = strline.split('\t')
            # just skip the intro rows
            if len(sline) > 1 and (sline[0]!= "pos"):
                res = int(sline[0][:-1])
                sift_res[res] =float(sline[1][:-1])
            line = fd.readline()

    with open(layout, "r") as fd:
        data = json.load(fd)
    sift = {}

    for sses in data["helices_residues"]:
        sift[sses] = {}
        res_int = data["helices_residues"][sses]
        r = res_int[1]
        while r < res_int[0]:
            av = sift_res[r]
            if av < 0.085:
                sift[sses][r] = 70
            elif av < 0.16:
                sift[sses][r] = 50
            else:
                sift[sses][r] = 20
            r += 1

    data["am"] = sift
    with open(new_layout, "w") as fd:
        json.dump(data, fd, indent=4, separators=(',', ':'))



def create_layout_primateai_color(layout, sift_file, new_layout):
    # TODO: comment
    if not os.path.isfile(sift_file):
        return False
    with open(sift_file, "r") as fd:
        sift_res = {}
        line = fd.readline()
        while line:
            strline = str(line)
            sline = strline.split('\t')
            # just skip the intro rows
            if len(sline) > 1 and (sline[0]!= "gene_name"):
                res = int(sline[1])
                sift_res[res] = float(sline[3][:-1])
            line = fd.readline()

    with open(layout, "r") as fd:
        data = json.load(fd)
    sift = {}
    for sses in data["helices_residues"]:
        sift[sses] = {}
        res_int = data["helices_residues"][sses]
        r = res_int[1]
        while r < res_int[0]:
#            print("layout", layout,r, sift_res.keys())
            if r not in sift_res:
                return False
            av = sift_res[r]
            if av <= 0.50:
                sift[sses][r] = 20
            elif av < 0.60:
                sift[sses][r] = 50
            else:
                sift[sses][r] = 70

            r += 1

    data["am"] = sift
    with open(new_layout, "w") as fd:
        json.dump(data, fd, indent=4, separators=(',', ':'))
    return True


def create_layout_charges_color(layout, ac_file, new_layout):
    with open(ac_file, "rb") as fd:
        ac_res = {}
        counter = 0
        line = fd.readline()
        start = False
#        if "Q9LU77" in ac_file:
#            start = True
        while line:
            counter += 1
            strline = str(line)[2:-2]
            if strline.startswith("TER"):
                if not start:
                    start = True
                else:
                    break
            if start and strline.startswith("ATOM") and len(line) < 80:
                if int(line[21:26]) not in ac_res:
                    ac_res[int(line[21:26])] = []
                ac_res[int(line[21:26])].append(float(line[55:62]))
            line = fd.readline()
    with open(layout, "r") as fd:
        data = json.load(fd)
    charges = {}
    for sses in data["helices_residues"]:
        charges[sses] = {}
        res_int = data["helices_residues"][sses]
        r = res_int[1]

        while r < res_int[0]:
            av = averager(ac_res[r])
            av = averager(ac_res[r])
            charges[sses][r] = int(av*100)
            r += 1
    data["charges"] = charges
    with open(new_layout, "w") as fd:
        json.dump(data, fd, indent=4, separators=(',', ':'))


def draw_am_diagrams(lf, prot, chain, domain, family, directory, suffix):
    svg_args = argparse.Namespace(
        add_descriptions="chain",
        layouts=[lf],
        ligands=False,
        opac='1',
        debug=False,
        irotation=None,
        itemplate=None,
        iligands=None,
        orotation=f"{directory}/rotation.json",
        otemplate=f"{directory}/{family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{new_dir}/image-{prot}_{chain}:{domain}_{suffix}.svg",
        channels=False,
        am=True,
        data=None,
        groups=None,
    )
    utils_generate_svg.compute(svg_args)


def draw_multi_am_diagrams(ld, family, directory, suffix):
    svg_args = argparse.Namespace(
        add_descriptions="multi, template",
        layouts=[ld],
        ligands=False,
        opac='0',
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        itemplate=None,
        iligands=None,
        orotation=None,
        otemplate=f"{directory}/{family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{new_dir}/{family}_{suffix}.svg",
        channels=False,
        am=True,
        data=None,
        groups=None,
        fill=None,
    )
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{directory}")


def draw_charges_diagrams(lf, prot, chain, domain, family, directory):
    svg_args = argparse.Namespace(
        add_descriptions="chain",
        layouts=[lf],
        ligands=False,
        opac='1',
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        itemplate=None,
        iligands=None,
        orotation=None,
        otemplate=f"{directory}/{family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{new_dir}/image-{prot}_{chain}:{domain}_charges.svg",
        channels=False,
        charges=True,
        data=None,
        groups=None,
        fill=None,
    )
    utils_generate_svg.compute(svg_args)


def draw_multi_charges_diagrams(ld, family, directory):
    svg_args = argparse.Namespace(
        add_descriptions="multi, template",
        layouts=[ld],
        ligands=False,
        opac='0',
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        itemplate=None,
        iligands=None,
        orotation=None,
        otemplate=f"{directory}/{family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{new_dir}/{family}_charges.svg",
        channels=False,
        charges=True,
        data=None,
        groups=None,
        fill=None,
    )
    utils_generate_svg.compute(svg_args)
    utils_common.convert_svg_to_png(svg_args.output, f"{directory}")


def sync_diagrams(family, directory, layout_list, not_ligands, verbose):
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
        otemplate=f"{directory}/{args.family}.json",
        oligands=f"{directory}/ligands.json",
        output=f"{directory}/{args.family}_sses.svg",
        data=None,
        groups=None,
        fill=None,
    )
    if not_ligands:
        svg_args.ligands = False
    link = f"https://2dprots.ncbr.muni.cz/files/family/{family}/rotation.json"
    ret = utils_common.curl_file(link, directory, f"{directory}/2dprots/rotation.json")
    logging.debug(f"Old version of rotation file exists: {ret}")
    if ret is False:
        logging.debug("no diagrams to synchronize with, create diagrams from scratch")
        if verbose:
            print("no diagrams to synchronize with, create diagrams from scratch")
        # we have to sync the new versions of diagrams with the previous one
        utils_generate_svg.compute(svg_args)
        utils_common.convert_svg_to_png(svg_args.output, f"{directory}")
    else:
        logging.debug("Try to synchronize with old version")
        if verbose:
            print("Try to synchronize with old version")
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
            if verbose:
                print("INFO: there is no common layout")
            command = f"curl from layout list {layout_list} (trying to sync layout rotation with old version)"
            utils_common.error_report(directory, command, None, ret, f"{directory}/curl-warning-{time.time()}.txt")
            utils_generate_svg.compute(svg_args)
            utils_common.convert_svg_to_png(svg_args.output, f"{directory}")
        else:
            logging.info(f"common layout was found {found}")
            if verbose:
                print(f"common layout was found {found}")
            # sync can be done
            aux_layout_fn = "layout-aaaaA00.json"
            aux_layout_full_fn = f"{directory}/{aux_layout_fn}"
            shutil.move(f"{directory}/{fn}", aux_layout_full_fn)
            # {directory}/{fn} will be the {directory}/2dprots/{fn} one with sses named the same way as layouts
            # in {directory}/layout-* files
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
    description='generate diagram for a alphafold structure')
parser.add_argument('--family', required=True, action="store",
                    help="family name")
parser.add_argument('--pdb_dir', action="store",
                    help="pdb and cif files source directory (tree structure e.g. ./0/2/Z/AF-A0A7W8CZ20-...{pdb,cif})")
parser.add_argument('--pdb_dir2', action="store",
                    help="pdb and cif files source directory (no structure e.g. ./AF-A0A7W8CZ20-...{pdb,cif})")
parser.add_argument('--output', action="store",
                    help="path to which will be put the result alphafold directory")
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debugging output')
parser.add_argument('-i', '--info', action='store_true',
                    help='print info output')
parser.add_argument('--am', required=False, action="store_true",
                    help="do alphamissense diagram as well")
parser.add_argument('--am_tmp', required=False, action="store",
                    help="cache directory for input and output alphmissense files")
parser.add_argument('--sift_tmp', required=False, action="store",
                    help="cache directory for input and output sift files")
parser.add_argument('--primateai_tmp', required=False, action="store",
                    help="cache directory for input and output sift files")
parser.add_argument('--only_am', required=False, action="store_true",
                    help="generte domains in with alphmissense files only")
parser.add_argument('--charges', required=False, action="store_true",
                    help="do alphamissense dyeing as well")
parser.add_argument('--charges_tmp', required=False, action="store",
                    help="cache directory for input and output charges files")
parser.add_argument('--channels_tmp', required=False, action="store",
                    help="do channels diagrams as well, grab the data from channels_tmp dir NOT WORKING NOW")
parser.add_argument('--detect', required=False, action="store_true",
                    help="do only sses detection")
parser.add_argument('--fast', required=False, action="store_true",
                    help="do optimisation to be as fast as possible (data quolity may be decreased)")
parser.add_argument('--cath_file', required=False, action="store",
                    help="optional cath - cath_assignments.tsv file location")
parser.add_argument('--annotation_dir', required=False, action="store",
                    help="annotation directory - like generated-.../annotated")
parser.add_argument('--layout_dir', required=False, action="store",
                    help="annotation directory - like generated-.../annotated")
parser.add_argument('--verbose', required=False, action="store_true",
                    help="verbose")
args = parser.parse_args()

# create new protein diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
name = args.family
if args.detect:
    new_dir = f"detected-{name}-{timestamp}"
else:
    new_dir = f"generated-{name}-{timestamp}"
if args.output:
    new_dir = f"{args.output}/{new_dir}"
os.makedirs(new_dir, exist_ok=True)
family = args.family

# start logger
logger = start_logger(args.debug, args.info, new_dir, name)
last_time = utils_common.time_log("0 Befor process starts", last_time)

# we want to download/read relevant cif files
cif_dir = f"{new_dir}/cif"
os.mkdir(cif_dir)

if not args.cath_file:
    cath_file = "cath/cath_assignments.tsv"
else:
    cath_file = args.cath_file
af_list = find_family_domains(args.family, new_dir, args.pdb_dir, cath_file, args.verbose)

if (len(af_list) == 0):
    last_time = utils_common.time_log("EMPTY", last_time)
    with open(f"{new_dir}/empty_{args.family}", "w") as fd:
        fd.close()
    sys.exit()

last_time = utils_common.time_log("01 family domains found", last_time)
if args.am and args.only_am:
    # go through all af_list and find which are in am database
    new_af_list = []
    for cd in af_list:
        protein = cd[0][:-3]
        am_fn = f'{args.am_tmp}/{protein}.am'
        if os.path.isfile(f"{am_fn}"):
            new_af_list.append(cd)
    if args.only_am:
        af_list = new_af_list
last_time = utils_common.time_log("02 am downloaded", last_time)

# download pdb and cif files if it is needed
if args.annotation_dir:
    logging.info("Skip pdb and cif downloading, no need we have annotation files")
    if args.verbose:
        print("Skip pdb and cif downloading, no need we have annotation files")
else:
    if args.pdb_dir:
        pdb_dir = args.pdb_dir
    else:
        pdb_dir = None

    if args.pdb_dir2:
        pdb_dir2 = args.pdb_dir2
    else:
        pdb_dir2 = None
    # download pdb and cif protein files from local dir cifdir, error messages go to new_dir
    af_list = simlink_cath_alphafold_list(af_list, args.pdb_dir, cif_dir, pdb_dir2, new_dir)
last_time = utils_common.time_log("03 cif and pdb downloaded", last_time)

# we have to download one layout from 2dprots to synchronize rotation
twodprots_layout = add_2dprot_layout(args.family, new_dir, cif_dir)

if twodprots_layout:
    af_list.append(twodprots_layout)

last_time = utils_common.time_log("031 cif and pdb downloaded", last_time)
if args.annotation_dir:
    annotation_dir = f"{new_dir}/annotation"
    os.mkdir(annotation_dir)

    res = simlink_cath_alphafold_list_annoation_files(af_list, args.annotation_dir, annotation_dir, new_dir)
    if res is False:
        logging.error("AF annotation files not found completely")
else:
    # create overprot input json
    overprot_input = f"{new_dir}/overprot_input.json"
    with open(overprot_input, "w") as fd:
        data = []
        json.dump(data, fd, indent=4, separators=(",", ": "))
    counter = 0
    for domain in af_list:
        # create overprot input json
        # + create input pdb and cifs
        counter += 1
        append_to_overprot_input_file(overprot_input, cif_dir, domain)
        if counter > 2000:
            break
    last_time = utils_common.time_log("032 cif and pdb downloaded", last_time)
    annotation_dir = f"{new_dir}/annotation"
    os.mkdir(annotation_dir)
    if args.channels_tmp:
        # channels data should be generated as well
        for domain in af_list:
            print("add_channel_file_to_ch_dir", domain)
            add_channel_file_to_ch_dir(annotation_dir, domain, new_dir)#, args.channels_tmp)
    last_time = utils_common.time_log("033 cif and pdb downloaded", last_time)
    # Annotate download protein structure
    if args.detect:
        counter = 0
        for domain in af_list:
            counter += 1
            if (counter % 1000) == 0:
                print(f"detected {counter} {domain}")
            utils_annotation.detect_chain_list(f"{new_dir}/overprot/cif_cealign/",
                                               cif_dir,
                                               [(f"{domain[0]}", "A", "", f"{domain[0]}A{domain[2]}")],
                                               f"{new_dir}")
    last_time = utils_common.time_log("034 cif and pdb downloaded", last_time)

    set_ubertemplate(f"{new_dir}/overprot", overprot_input, new_dir, cif_dir, args.fast)

    last_time = utils_common.time_log("04 ubertemplate set", last_time)

    shutil.copyfile(f"{new_dir}/overprot/results/consensus.cif", f"{new_dir}/overprot/cif_cealign/consensus.cif")
    shutil.copyfile(f"{new_dir}/overprot/results/consensus.sses.json", f"{new_dir}/overprot/cif_cealign/consensus-template.sses.json")
    counter = 0
    for domain in af_list:
        counter += 1
        if (counter % 1000) == 0:
            print(f"annotated {counter}")

        if domain[3] != "2dprots":
            # af domain

            shutil.copyfile(f"{cif_dir}/{domain[0]}.cif", f"{new_dir}/overprot/cif_cealign/{domain[0]}.cif")
            domain_boundaries_tag = transform_list_to_db(domain[1])
            utils_annotation.annotate_with_given_template_list(f"{new_dir}/overprot/cif_cealign/",
                                                               [(f"{domain[0]}", "A", domain_boundaries_tag, f"{domain[0]}A{domain[2]}")],
                                                               f"{new_dir}",
                                                               domain_boundaries="Manual",
                                                               prefix_secstr_dir="../overprot/OverProt/SecStrAnnot2/")
            shutil.copyfile(f"{new_dir}/overprot/cif_cealign/{domain[0]}A00-annotated.sses.json", f"{annotation_dir}/{domain[0]}A{str(domain[2]).zfill(2)}-annotated.sses.json")
        else:
           # 2dprots domain record added to sync rotation with 2dprots dbs
            shutil.copyfile(f"{cif_dir}/{domain[0]}.cif", f"{new_dir}/overprot/cif_cealign/{domain[0]}.cif")
            domain_boundaries_tag = transform_range_list_to_db(domain[4])
            sa_id = utils_pdbe_api.transform_cath_chain_id_to_struct_asym_id(domain[0], domain[1])

            ret = utils_annotation.annotate_with_given_template_list(f"{new_dir}/overprot/cif_cealign/",
                                                                     [(f"{domain[0]}", sa_id, domain_boundaries_tag, f"{domain[0]}{domain[1]}{domain[2]}")],
                                                                     f"{new_dir}",
                                                                     domain_boundaries="Manual",
                                                                     prefix_secstr_dir="../overprot/OverProt/SecStrAnnot2/")

            shutil.copyfile(f"{new_dir}/overprot/cif_cealign/{domain[0]}{sa_id}00-annotated.sses.json", f"{annotation_dir}/{domain[0]}{domain[1]}{str(domain[2]).zfill(2)}-annotated.sses.json")
    shutil.copyfile(f"{new_dir}/overprot/results/diagram.json", f"{annotation_dir}/diagram.json")

last_time = utils_common.time_log("05 annotation complete", last_time)

# now we have to go through all chains
# found out all domain ranges which are already their diagrams
# and aggregate this intervals
# count layoutn of chain  "A" of imput alpha fold enntity
generate_cl(cif_dir, [af_list[0]], "A", new_dir, annotation_dir, new_dir, None, args)
start_layout_filename = f"{new_dir}/start-layout.json"

shutil.move(f"{new_dir}/layout-{af_list[0][0]}-{af_list[0][1][0]}-{af_list[0][1][1]}.json", start_layout_filename)

generate_cl(cif_dir, af_list, "A", new_dir, annotation_dir, new_dir, start_layout_filename, args)
last_time = utils_common.time_log("06 layout computed", last_time)
# for af in af_list:
#    if af[3] != "2dprots":
#        r = utils_alphafold_confidence.switch_color_to_confidence(f"{new_dir}/layout-{af[0]}-{af[1][0]}-{af[1][1]}.json", f"{cif_dir}/{af[0]}.pdb", f"{new_dir}/overprot/cif_cealign/consensus.sses.json")

last_time = utils_common.time_log("07 generation complete", last_time)
####
# color sses based on sses variance in family
#
# for file_name in os.listdir(new_dir):
#     if file_name.startswith("layout-"):
#         color_me.switch_color_to_stats(f"{new_dir}/{file_name}", f"{new_dir}/sses-stats.json")
# overprot_link = "https://overprot.ncbr.muni.cz/data/db/"
# link_variance_data = f"http://overprot.ncbr.muni.cz/data/db/family/consensus_sses//consensus-{args.family}.sses.json"
# consensus_file = f"{new_dir}/consensus-{args.family}.sses.json"
# ret = utils_common.curl_file(link_variance_data, new_dir, consensus_file)
# if not ret:
#
#     sys.exit()

# consensus_file = f"{new_dir}/overprot/results/consensus.sses.json"
# for file_name in os.listdir(new_dir):
#     if file_name.startswith("layout-"):
#         color_me.switch_color_to_variance(f"{new_dir}/{file_name}", consensus_file)
####
# for file_name in os.listdir(new_dir):
#     if file_name.startswith("layout-"):
#         utils_alphafold_confidence.switch_htype_to_one(f"{new_dir}/{file_name}")
#
consensus_file = f"{new_dir}/overprot/results/consensus.sses.json"
# ret = utils_common.curl_file(link_variance_data, new_dir, consensus_file)
for file_name in os.listdir(new_dir):
    if file_name.startswith("layout-"):
        utils_alphafold_confidence.switch_htype(f"{new_dir}/{file_name}", consensus_file)

# TODO: we have to use the same rotation as 2dprot cath dbs
# synchronize names in old and new layout
sync_diagrams(args.family, new_dir, af_list, False, args.verbose)

if twodprots_layout:
    os.remove(f"{new_dir}/layout-{twodprots_layout[0]}{twodprots_layout[1]}{twodprots_layout[2]}.json")
    af_list.remove(twodprots_layout)

new_layout_dir = f"{new_dir}/new_layout"
os.mkdir(new_layout_dir)

# use this rotation file
multiple_diagrams(annotation_dir, new_dir, False, False)

# draw single image per chain and whole together
for af in af_list:
    if af[3] == "2dprots":
        layout = f"{new_dir}/layout-{af[0]}{af[1]}{af[2]}.json"
    else:
        layout = f"{new_dir}/layout-{af[0]}-{af[1][0]}-{af[1][1]}.json"
    svg_args = argparse.Namespace(
        add_descriptions="1",
        data=False,
        debug=False,
        irotation=f"{new_dir}/rotation.json",
        otemplate=None,
        layouts=[layout],
        opac='0',
        orotation=None,
        itemplate=None,
        iligands=None,
        ligands=None,
        oligands=None,
        groups=None,
        output=f"{new_dir}/image-{af[0]}-{transform_list_to_db(af[1])}.svg"
    )
    utils_generate_svg.compute(svg_args)
last_time = utils_common.time_log("3 image created", last_time)

counter = 0
counter_all = 0
if args.am:
    logging.info("Generate AlphaMissense diagrams")

    # if alphamissense output is set, generate am diagrams as well
    for cd in af_list:
        counter_all += 1
        protein = cd[0][:-3]
        am_fn = f'{args.am_tmp}/{protein}.am'
        if os.path.isfile(f"{am_fn}"):
            counter += 1
            # we have charges file - create layout with per res coloring (charges tab)
            layout_fn = f"layout-{protein}-F1-{cd[1][0]}-{cd[1][1]}.json"

            if os.path.isfile(f"{new_dir}/{layout_fn}"):
                create_layout_am_color(f"{new_dir}/{layout_fn}", am_fn, f"{new_layout_dir}/{layout_fn}")
                draw_am_diagrams(f"{new_layout_dir}/{layout_fn}", protein, cd[1][0], cd[1][1], family, new_layout_dir, "am")
            else:
                logging.error(f"no {new_dir}/{layout_fn} entry")

    if counter > 0:
        draw_multi_am_diagrams(f"{new_layout_dir}", args.family, new_layout_dir, "am")

if args.sift_tmp:
    logging.info("Generate sift diagrams")
    stif_layout_dir = f"{new_dir}/stif_layout"
    os.mkdir(stif_layout_dir)

    # if sift output is set, generate am diagrams as well
    for cd in af_list:
        counter_all += 1
        protein = cd[0][:-3]
        sift_fn = f'{args.sift_tmp}/{protein}.txt'
        if os.path.isfile(f"{sift_fn}"):
            counter += 1
            # we have charges file - create layout with per res coloring (charges tab)
            layout_fn = f"layout-{protein}-F1-{cd[1][0]}-{cd[1][1]}.json"

            if os.path.isfile(f"{new_dir}/{layout_fn}"):
                create_layout_sift_color(f"{new_dir}/{layout_fn}", sift_fn, f"{stif_layout_dir}/{layout_fn}")
                draw_am_diagrams(f"{stif_layout_dir}/{layout_fn}", protein, cd[1][0], cd[1][1], family, stif_layout_dir, "sift")
            else:
                logging.error(f"no {new_dir}/{layout_fn} entry")

    if counter > 0:
        draw_multi_am_diagrams(f"{stif_layout_dir}", args.family, stif_layout_dir, "sift")

if args.primateai_tmp:
    logging.info("Generate primateai diagrams")
    prim_layout_dir = f"{new_dir}/prim_layout"
    os.mkdir(prim_layout_dir)

    # if primateai output is set, generate am diagrams as well
    for cd in af_list:
        counter_all += 1
        protein = cd[0][:-3]
        primateai_fn = f'{args.primateai_tmp}/{protein}.txt'
        if os.path.isfile(f"{sift_fn}"):
            counter += 1
            # we have charges file - create layout with per res coloring (charges tab)
            layout_fn = f"layout-{protein}-F1-{cd[1][0]}-{cd[1][1]}.json"

            if os.path.isfile(f"{new_dir}/{layout_fn}"):
                ret = create_layout_primateai_color(f"{new_dir}/{layout_fn}", primateai_fn, f"{prim_layout_dir}/{layout_fn}")
                if ret: 
                    draw_am_diagrams(f"{prim_layout_dir}/{layout_fn}", protein, cd[1][0], cd[1][1], family, prim_layout_dir, "primateai")
            else:
                logging.error(f"no {new_dir}/{layout_fn} entry")

    if counter > 0:
        draw_multi_am_diagrams(f"{prim_layout_dir}", args.family, new_layout_dir, "primateai")

last_time = utils_common.time_log("08 am diagrams complete", last_time)

if args.charges:
    logging.info("Generate charges diagrams")
    charges_dir = f"{new_dir}/charges"
    os.mkdir(charges_dir)
    charges_layout_dir = f"{new_dir}/charges-layout"
    os.mkdir(charges_layout_dir)
    ccc = 0
    for cd in af_list:
        protein = cd[0][:-3]
        if args.charges_tmp:
            charges_fn = f'{args.charges_tmp}/{protein}.alphacharges'
        else:
            charges_fn = f'{new_dir}/{protein}.alphacharges'
        if not args.charges_tmp or not os.path.isfile(f"{charges_fn}"):
            # charges buffer was not set - we have to download charges file ourself
            # try to curl them
            link_fn = f"https://alphacharges.ncbr.muni.cz/download_files?ID={protein}_7.2_4"
            utils_common.curl(new_dir, err_report=True, options=f'-L --output /dev/null -X POST "https://alphacharges.ncbr.muni.cz/" -H "Content-Type: application/x-www-form-urlencoded" -d "ph=7.2" -d "prediction_version=4" -d "code={protein}" -d "action=calculate charges"')
            utils_common.curl(new_dir, err_report=True, options=f'-L --max-time 10 -X POST "https://alphacharges.ncbr.muni.cz/calculation?ID={protein}_7.2_4"')
            counter = 0
            while True:
                utils_common.curl(new_dir, err_report=False, options=f'--output {new_dir}/pomXXX  "https://alphacharges.ncbr.muni.cz/progress?ID={protein}_7.2_4"')
                is_in = False
                with open(f"{new_dir}/pomXXX", "r") as fd:
                    line = fd.readline()
                    while line:
                        if line.startswith("<p><strong> Step 6/6:</strong> Partial atomic charges calculated. "):
                            is_in = True
                        line = fd.readline()
                if is_in:
                    ret = True
                    break
                counter += 1
                time.sleep(10)
                if counter > 1:
                    ret = False
                    break
            if ret:
                utils_common.curl_file(link_fn, "new_dir", charges_fn)
                shutil.copyfile(charges_fn, f"{charges_dir}/{protein}.alphacharges")
        else:
            # it is in buffer copy it
            ret = True
            shutil.copyfile(charges_fn, f"{charges_dir}/{protein}.alphacharges")
        ccc = 0
        if ret:
            # we have charge file now - just draw the diagram based on charges
            layout_fn = f"layout-{protein}-F1-{cd[1][0]}-{cd[1][1]}.json"
            create_layout_charges_color(f"{new_dir}/{layout_fn}", f"{charges_dir}/{protein}.alphacharges", f"{charges_layout_dir}/{layout_fn}")
            draw_charges_diagrams(f"{charges_layout_dir}/{layout_fn}", protein, cd[1][0], cd[1][1], family, new_layout_dir)
            ccc += 1
    if ccc > 0:
        draw_multi_charges_diagrams(f"{charges_layout_dir}", args.family, new_layout_dir)

last_time = utils_common.time_log("09 charges diagrams complete", last_time)

if (not args.debug):
    shutil.rmtree(f"{cif_dir}")
    shutil.rmtree(f"{annotation_dir}")
last_time = utils_common.time_log("4 OK", start_time)
