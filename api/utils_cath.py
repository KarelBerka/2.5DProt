#!/usr/bin/env python3

# go through start_layout list, compare it with given protein-domain-chain name
# if there is a protein-domain-chain with rmsd < c, then return pxdxc layout
# else return -1
#
# Example:
# ./sort_to_prefamilies.py ./processed/generated-1.10.10.10

import os
import json
import contextlib
import logging
import urllib
import shutil

import utils_common

cath_version = "v4_3_0"
cath_lv_release = f'http://download.cathdb.info/cath/releases/all-releases/{cath_version}'
cath_latest_release = f'http://download.cathdb.info/cath/releases/latest-release'
cath_daily_release = f"http://download.cathdb.info/cath/releases/daily-release/newest/" # cath-b-newest-latest-release.gz


def cath_s35_cluster_dbs_gen(directory):
    # download the last cath version of s35 dbs and generate per family dbs
    # to "directory directory - it is necessary to do this after new cath version
    url = os.path.join(cath_lv_release, 'cath-classification-data', f'cath-domain-list-{cath_version}.txt')
    path = os.path.join(directory, f'cath-domain-list-{cath_version}.txt')
    logging.debug("Downloading %s to %s" % (url, path))
    utils_common.curl_file(url, directory, path, False, "")

    # parse the cath s35 dbs
    logging.debug("Parsing file %s" % path)
    cath_domain_list = {}
    with open(path, 'r') as f:
        for row in f.read().split('\n'):
            if row == '' or row[0] == '#':
                continue
            row_split = row.split()
            family = ".".join(row_split[1:5])
            s35 = ".".join(row_split[1:6])
            if family not in cath_domain_list:
                cath_domain_list[family] = {}
            if s35 not in cath_domain_list[family]:
                cath_domain_list[family][s35] = [row_split[0]]
            else:
                cath_domain_list[family][s35].append(row_split[0])
    logging.debug("Parsed file %s with %s entries" % (path, len(cath_domain_list)))
    os.remove(path)
    os.mkdir(f"{directory}/s35")
    for fam in cath_domain_list:
        with open(f"{directory}/s35/{fam}_s35.json", "w") as fd:
            json.dump(cath_domain_list[fam], fd, indent=4, separators=(",", ":"))

    return cath_domain_list


def cath_s35clusters(directory):
    # return actual dbs of s35 clusters

    # download the last cath version of s35 dbs
    url = os.path.join(cath_lv_release, 'cath-classification-data', f'cath-domain-list-{cath_version}.txt')
    path = os.path.join(directory, f'cath-domain-list-{cath_version}.txt')

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

    return cath_domain_list


def cath_family_list(directory, leave=False, cfile=False):
    # return actual dbs of cath families clusters
    # download the last cath version of dbs

    if cfile == False:
        dbs_name = 'cath-b-newest-latest-release.gz'

        url = os.path.join(cath_daily_release, dbs_name)
        path = os.path.join(directory, dbs_name)

        logging.debug("Downloading %s to %s" % (url, path))
        try:
            with contextlib.closing(urllib.request.urlopen(url)) as r:
                with open(path, 'wb') as f:
                    shutil.copyfileobj(r, f)
        except urllib.error.URLError as e:
            if 'family dbs Failed to change directory.' in e.reason:
                logging.error("Failed to download %s because it is not present on CATH servers" % url)
                raise FileNotFoundError(e.reason)
            else:
                raise

        # we havbe to unzip cath dbs
        utils_common.gunzip_archive(path, directory)
        path = path[:-3]

        # parse the cath s35 dbs
        logging.debug("Parsing file %s" % path)
    else:
         dbs_name = 'cath-b-newest-latest-release'
         path = os.path.join(directory, dbs_name)
         shutil.copyfile(cfile, path)

    cath_family_list = set()

    with open(path, 'r') as f:
        for row in f.read().split('\n'):
            if row == '' or row[0] == '#':
                continue
            row_split = row.split(" ")
            cath_family_list.add(row_split[2])
    logging.debug("Parsed file %s with %s entries" % (path, len(cath_family_list)))
    if not leave:
        os.remove(path)

    return cath_family_list


def cath_domains_of_given_family(directory, family, leave=False, cfile=False):
    # return actual dbs of cath families clusters
    # download the last cath version of dbs

    if cfile == False:
        dbs_name = 'cath-b-newest-latest-release.gz'

        url = os.path.join(cath_daily_release, dbs_name)
        path = os.path.join(directory, dbs_name)

        logging.debug("Downloading %s to %s" % (url, path))
        try:
            with contextlib.closing(urllib.request.urlopen(url)) as r:
                with open(path, 'wb') as f:
                    shutil.copyfileobj(r, f)
        except urllib.error.URLError as e:
            if 'family dbs Failed to change directory.' in e.reason:
                logging.error("Failed to download %s because it is not present on CATH servers" % url)
                raise FileNotFoundError(e.reason)
            else:
                raise

        # we have to unzip cath dbs
        utils_common.gunzip_archive(path, directory)
        path = path[:-3]

        # parse the cath s35 dbs
        logging.debug("Parsing file %s" % path)
    else:
         dbs_name = 'cath-b-newest-latest-release'
         path = os.path.join(directory, dbs_name)
         shutil.copyfile(cfile, path)

    cath_domain_list = set()

    with open(path, 'r') as f:
        for row in f.read().split('\n'):
            if row == '' or row[0] == '#':
                continue
            row_split = row.split(" ")
            if row_split[2] == family:
                cath_domain_list.add(row_split[0])
    logging.debug("Parsed file %s with %s entries" % (path, len(cath_domain_list)))
    if not leave:
        os.remove(path)

    return cath_domain_list
