#!/usr/bin/env python3

# go through start_layout list, compare it with given protein-domain-chain name
# if there is a protein-domain-chain with rmsd < c, then return pxdxc layout
# else return -1
#
# Example:
# ./sort_to_prefamilies.py ./processed/generated-1.10.10.10

import os
import contextlib
import logging
import urllib
import shutil
import json
import utils_common


def sort_to_s35clusters(directory, family):
    # download the last version of s35 dbs
    url = f"https://2dprots.ncbr.muni.cz/pub/s35/2022-07-24/{family}_s35.json"
    path = f'{directory}/{family}_s35.json'
    logging.debug("Downloading %s to %s" % (url, path))
    utils_common.curl_file(url, directory, path, False, "")
    # parse per family cath s35 dbs and create all s35 clusters
    with open(path, 'r') as fd:
        data = json.load(fd)

    for s35 in data:
        path = directory + "/" + s35
        if not os.path.isdir(path):
                os.mkdir(path)
        for p in data[s35]:
            if os.path.isfile(f"{path}/../layout-{p}.json"):
                os.symlink(f"../layout-{p}.json", f"{path}/layout-{p}.json")
            else:
                logging.warning(f"Domain {p} is not in family {family}, but it is in s35 {s35}")


def sort_to_s35clusters_cath(directory, layout_list):
    # for all pcd in layout_list find out, in which s35 cluster it belonges

    # TODO - fix the last release
    # download the last cath version of s35 dbs
    cath_release = 'http://download.cathdb.info/cath/releases/all-releases/v4_3_0'

    url = os.path.join(cath_release, 'cath-classification-data', 'cath-domain-list-v4_3_0.txt')
    path = os.path.join(directory, 'cath-domain-list-v4_3_0.txt')

    logging.debug("Downloading %s to %s" % (url, path))
    try:
        with contextlib.closing(urllib.request.urlopen(url)) as r:
            with open(path, 'wb') as f:
                shutil.copyfileobj(r, f)
    except urllib.error.URLError as e:
        if '550 Failed to change directory.' in e.reason:
            logging.error("Failed to download %s because it is not present on CATH servers" % url)
            raise FileNotFoundError(e.reason)
        else:
            raise

    # parse the cath s35 dbs
    logging.debug("Parsing file %s" % path)
    cath_domain_list = {}
    with open(path, 'r') as f:
            for row in f.read().split('\n'):
                if row == '' or row[0] == '#':
                    continue
                row_split = row.split()
                cath_domain_list[row_split[0]] = ".".join(row_split[1:6])
    logging.debug("Parsed file %s with %s entries" % (path, len(cath_domain_list)))
    os.remove(path)

    for l in layout_list:
        # for all pcds find out its s35 cluster in parsed list
        protein = l[0]
        cath_chain = l[1]
        domain = l[2]
        f = f"layout-{l[0]}{l[1]}{l[2]}.json"

        s35 = cath_domain_list[protein + cath_chain + domain]
        if s35 is not None:
            path = directory + "/" + s35
            if not os.path.isdir(path):
                os.mkdir(path)
            os.symlink("../"+f, path+"/"+f)

         # TODO: test s35 starts by family name, is not None

