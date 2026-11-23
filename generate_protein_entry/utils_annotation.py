#!/usr/bin/env python3

import os
import subprocess
import time
import logging
import sys
import utils_common
import utils_pdbe_api
import shutil


def detect_list(directory, cif_directory, pcd_list, log_directory, whole_protein=False, prefix=".."):
    # detect protein domain list
    # call SetStrAnnotator to detect sses in all sse from protein_list
    # pcd_list is a list of triples protein, chain, domain
    os.environ['DOTNET_CLI_TELEMETRY_OUTPUT'] = "1"

    detected_pcd = 0
    for pcd in pcd_list:
        if whole_protein is False:
            prot, chain, domain = pcd
            # find domain boundaries
            db = utils_pdbe_api.get_domain_boundaries(prot, chain, domain)
            pcd_string = f"{prot},{chain},{db}"
            cif_file = f"{cif_directory}/{prot}.cif"
            if not os.path.isfile(cif_file):
                os.symlink(f"../cif/{prot}.cif", cif_file)
        elif whole_protein == "cath_alpha_fold":
            prot, chain, domain = pcd
            cif_file = f"{cif_directory}/{prot}.cif"
            if not os.path.isfile(cif_file):
                os.symlink(f"../cif/{prot}.cif", cif_file)
            pcd_string = f"{prot},{chain},{domain}"
        elif whole_protein == "cath_alpha_fold_family":
            prot = pcd[0]
            chain = "A"
            domain = f"{pcd[1][0]}:{pcd[1][1]}"
            cif_file = f"{cif_directory}/{prot}.cif"

            if not os.path.isfile(cif_file):
                os.symlink(f"../cif/{prot}.cif", cif_file)
            pcd_string = f"{prot},{chain},{domain}"
        elif whole_protein == "pdbe_best_structures":
            prot = pcd[0]
            chain = pcd[1]
            domain = f"{pcd[2][0]}:{pcd[2][1]}"
            cif_file = f"{cif_directory}/{prot}.cif"
            if not os.path.isfile(cif_file):
                os.symlink(f"../cif/{prot}.cif", cif_file)

            pcd_string = f"{prot},{chain},{domain}"
        elif whole_protein == "pf":
            cif_file = f"{cif_directory}/{pcd}.cif"
            if not os.path.isfile(cif_file):
                os.symlink(f"../cif/{pcd}.cif", cif_file)
            pcd_string = pcd
        else:
            cif_file = f"{cif_directory}/{pcd}.cif"

            if not os.path.isfile(cif_file):
                os.symlink(f"../cif/{pcd}.cif", cif_file)
            pcd_string = pcd

        # run SecStrAnnotator.ddl
        command = f"dotnet {prefix}/overprot/OverProt/SecStrAnnot2/SecStrAnnotator.dll --onlyssa {directory} {pcd_string} -a cealign --fallback 120 --verbose"
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            outs, errs = p.communicate(timeout=3600)
        except subprocess.TimeoutExpired:
            p.kill()
            outs, errs = proc.communicate()

        if p.returncode != 0:
            # if the first version is not succesfull try the second version
            command = f"dotnet {prefix}/overprot/OverProt/SecStrAnnot2/SecStrAnnotator.dll --onlyssa {directory} {pcd_string} -a none --fallback 120 --verbose"
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            try:
                outs, errs = p.communicate(timeout=3600)
            except subprocess.TimeoutExpired:
                p.kill()
                outs, errs = proc.communicate()
            if p.returncode != 0:
                log_fn = f"{log_directory}/SecStrAnnotator-fail-{time.time()}.txt"
                utils_common.error_report(directory, command, outs, errs, log_fn)
                logging.error(f"SecStrAnnotator on {pcd_string} failed, report save to {log_fn}")
                continue

        detected_pcd += 1

    return detected_pcd


def detect_chain_list(directory, cif_directory, pcd_list, log_directory):
    # detect protein chain list
    # call SetStrAnnotator to detect sses in all sse from protein_list
    os.environ['DOTNET_CLI_TELEMETRY_OUTPUT'] = "1"

    detected_pcd = 0
    for pcd in pcd_list:
        cif_file = f"{cif_directory}/{pcd[0]}.cif"
        if not os.path.isfile(cif_file):
            os.symlink(f"../cif/{pcd[0]}.cif", cif_file)
        pcd_string = f"{pcd[0]},{pcd[1]},:"

        # run SecStrAnnotator.ddl
        command = f"dotnet ../SecStrAnnotator.dll --onlyssa {directory} {pcd_string} -a none --fallback 120 --verbose"
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        r = p.wait()
        if (r != 0):
            out, err = p.communicate()
            log_fn = f"{log_directory}/SecStrAnnotator-fail-{time.time()}.txt"
            utils_common.error_report(directory, command, out, err, log_fn)
            logging.error(f"SecStrAnnotator on {pcd_string} failed, report save to {log_fn}")
            continue
        shutil.move(f"{directory}/{pcd[0]}-detected.sses.json", f"{directory}/{pcd[0]}{pcd[1]}00-annotated.sses.json")
        detected_pcd += 1

    return detected_pcd


def annotate_list(directory, pcd_list, log_dir, family, domain_boundaries=True, prefix_secstr_dir="../"):

    # call SetStrAnnotator to annotate sses in all sse from protein_list
    SSAnnot_metrictype = "3"
    # TODO: for family 3.20.20.80 we have to use dp here, mom is too slow (2021-10-30) - maybe will be fixed, but have to be tested
    if (family == "3.20.20.80") or (family == "3.40.50.300") \
        or (family == "3.40.720.10") or (family == "2.140.10.30") \
        or (family == "1.50.10.10") or (family == "2.130.10.10") \
        or (family == "3.80.10.10") or (family == "3.40.50.620") \
        or (family == "3.40.640.10") or (family == "3.50.50.60") \
        or (family == "3.20.20.470") or (family == "2.160.20.10") \
        or (family == "3.30.470.20") or (family == "3.40.190.10") \
        or (family == "3.40.50.1820") or (family == "1.10.640.10") \
        or (family == "3.40.50.2000") or (family == "short"):
        SSAnnot_matching = "dp"
    else:
        SSAnnot_matching = "mom"

    annotated_pcd = []
    for pcd in pcd_list:
        prot, chain, domain, domain_name = pcd
        if domain_boundaries:
            domain_range = utils_pdbe_api.get_domain_boundaries(prot, chain, domain)
            domain_range_tag = f",{domain_range}"
        else:
            domain_range_tag = ""
        logging.info(f"INFO: domain {domain_name} is annotating ")
        # run SecStrAnnotator.ddl
        # use -a cealign only if input cif files are not from ubertemplate directly otherwise use -a none
        command = f"dotnet {prefix_secstr_dir}SecStrAnnotator.dll --metrictype {SSAnnot_metrictype} --unannotated -m {SSAnnot_matching} {directory} consensus,A {domain_name},{chain}{domain_range_tag} -a cealign --fallback 120"
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        r = p.wait()
        if (r != 0):
            out, err = p.communicate()
            log_fn = f"{log_dir}/SecStrAnnotator-fail-{time.time()}.txt"
            utils_common.error_report(directory, command, out, err, log_fn)
            logging.error(f"SecStrAnnotator on {pcd} failed, report save to {log_fn}")
            sys.exit()
        else:
            annotated_pcd.append(pcd)
    return annotated_pcd


def annotate_with_given_template_list(directory, pcd_list, log_dir, domain_boundaries="Automatic", prefix_secstr_dir="../"):
    # call SetStrAnnotator to annotate sses in all sse from protein_list
    SSAnnot_metrictype = "3"
    SSAnnot_matching = "dp"     # shorter
    # SSAnnot_matching="mom"

    annotated_pcd = []
    counter = 0
    for pcd in pcd_list:
        counter += 1
        prot, chain, domain, domain_name = pcd
        if domain_boundaries == "Automatic":
            domain_range = utils_pdbe_api.get_domain_boundaries(prot, chain, domain)
            domain_range_tag = f",{domain_range}"
        elif domain_boundaries == "Manual":
            domain_range_tag = domain
        else:
            domain_range_tag = ":"
        logging.info(f"INFO: domain {domain_name} is annotating ")
        # run SecStrAnnotator.ddl
        # use -a cealign only if input cif files are not from ubertemplate directly otherwise use -a none
        command = f"dotnet {prefix_secstr_dir}SecStrAnnotator.dll --metrictype {SSAnnot_metrictype} --unannotated -m {SSAnnot_matching} {directory} consensus,A {prot},{chain},{domain_range_tag} -a cealign --fallback 120"
        if counter % 10000:
            print(counter, "command", command)
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()
        if (r != 0):
            log_fn = f"{log_dir}/SecStrAnnotator-fail-{time.time()}.txt"
            utils_common.error_report(directory, command, out, err, log_fn)
            logging.error(f"SecStrAnnotator on {pcd} failed, report save to {log_fn}")
#            sys.exit()
        else:
            shutil.move(f"{directory}/{pcd[0]}-annotated.sses.json", f"{directory}/{pcd[0]}{pcd[1]}00-annotated.sses.json")
            annotated_pcd.append(pcd)
    return annotated_pcd
