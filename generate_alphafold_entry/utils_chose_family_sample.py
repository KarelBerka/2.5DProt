#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import utils_common
import json
import os


def chose_family_sample(file_name, family, err_dir, max_dom, cut=0):
    # choose only up to max_dom doains from each organism (taxId), andomit (cut) shortest and longest domains as well (probably bugs)
    # if cut is set at first read the sizes of all domains to set the max and min values
    res_sizes = {}
    all_dom = 0
    if cut != 0:
        # find out the min_limit and max_limit - the limits of domain length
        # find the res sizes distribution
        with open(file_name, "r") as fd:
            line = fd.readline()
            # loooks like
            # "af_W1Q4I0_327_638	311	cath_ifs	ccatfp product	7	3.20.20.70	CATH_HMM_S95_SSAP	98.31	FALSE	31	39.53	G..."
            while line:
                sline = line.split("\t")  # record like af_W1Q4I0_327_638
                ssline = sline[0].split("_")
                if len(ssline) > 3:
                    all_dom += 1
                    res_s = int(ssline[3]) - int(ssline[2])
                    if res_s not in res_sizes:
                        res_sizes[res_s] = 1
                    else:
                        res_sizes[res_s] += 1
                line = fd.readline()
        res_sizes_values = list(res_sizes.keys())
        res_sizes_values.sort()

        limit = int(cut * all_dom / 100)

        # find the minimal limit ("cut" percent of all)
        i = 0
        min_sum = 0
        while ((i+1) < len(res_sizes_values)):
            min_sum += res_sizes[res_sizes_values[i]]
            if min_sum > limit:
                min_limit = res_sizes_values[i]-1
                break
            i += 1

        # find the maximal limit ("cut" percent of all)
        i = 1
        max_sum = 1
        while (i < len(res_sizes_values)):
            max_sum += res_sizes[res_sizes_values[len(res_sizes_values)-i]]
            if max_sum > limit:
                max_limit = res_sizes_values[len(res_sizes_values)-i]-1
                break
            i += 1

    domain_list = {}
    counter = 0
    aa = 0
    # go through the data set and ommit unnecessary ones
    # the minimal and maximal limits
    # only the max_dom pcs from one organism
    with open(file_name, "r") as fd:
        line = fd.readline()
        # loooks like
        # "af_W1Q4I0_327_638	311	cath_ifs	ccatfp product	7	3.20.20.70	CATH_HMM_S95_SSAP	98.31	FALSE	31	39.53	G..."
        while line:

            aa += 1
            sline = line.split("\t")  # record like af_W1Q4I0_327_638
            ssline = sline[0].split("_")
            if (sline[5] != family) and (family != "ALL") and (family != ""):
                # the record is about doain from another family
                line = fd.readline()
                continue

            if (cut != 0):
                # we have to check residue size if it is in limit
                if (len(ssline) > 3):
                    res = int(ssline[3]) - int(ssline[2])
                    if (res < min_limit) or (res > max_limit):
                        line = fd.readline()
                        continue
                else:
                    line = fd.readline()
                    continue

            # the size of domain is ok - we have to test organism
            link = f"https://alphafold.ebi.ac.uk/api/prediction/{ssline[1]}"
            lfn = f"{ssline[1]}.json"

            ret = utils_common.curl_file(link, err_dir, lfn, err_report=True)
            if ret:
                with open(lfn, "r") as fl:
                    data = json.load(fl)
                    taxid = data[0]["taxId"]
                    if (taxid not in domain_list):
                        domain_list[taxid] = [sline[0]]
                        counter += 1
                    elif (len(domain_list[taxid]) < max_dom):
                        domain_list[taxid].append(sline[0])
                        counter += 1
                os.remove(lfn)
            line = fd.readline()
    return domain_list


print("start")
chose_family_sample("af_3.20.20.70_adolases/copy/cath_assignment-3.20.20.700.tsv", "3.20.20.700", "./", 2, cut=2)
# chose_family_sample("af_3.20.20.70_adolases/copy/cath_assignment-3.20.20.7000s.tsv", "3.20.20.7000s", "./", 2, cut = 2)
# chose_family_sample("af_3.20.20.70_adolases/copy/cath_assignment-3.20.20.7000s.tsv", "3.20.20.7000s", "./", 2, cut = 10)
print("end")
