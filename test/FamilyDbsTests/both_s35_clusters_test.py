#!/usr/bin/env python3

import time
import argparse
import os
import subprocess
import json
import filecmp
import sys

sys.path.append('..')
import loader
import two_d_prot.pdbe_api


def curl_file(link, directory, filename, err_report=True, option=""):
    # curl "link" file to given "filename"
    command = f"curl {option} -f {link} -o {filename} >/dev/null 2>/dev/null"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    r = p.wait()

    if (r != 0):
        if err_report is True:
            out, err = p.communicate()
            log_fn = f"{directory}/curl-fail-{time.time()}.txt"
            error_report(directory, command, out, err, log_fn)
            print(f"curl of {link} failed, report save to {log_fn}")
        return False
    return True


def read_cath_family_list():
    # ./cath-family-updates.py --update-putative-only --update-from-release $(date -d "2020-10-20" "+%Y-%m-%d" )
    # ./cath-family-updates.py --update-from-release $(date -d "2020-10-20" "+%Y-%m-%d" )         - all families
    # command = f'./cath-family-updates.py --update-from-release $(date -d "2020-10-20" "+%Y-%m-%d" )'
    filename = "cath-b-newest-all.gz"
    link = f"http://download.cathdb.info/cath/releases/daily-release/newest/{filename}"

    ret = curl_file(link, "", filename, True, "")
    if ret is True:
        command = f"gunzip {filename}"
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        out, err = p.communicate()
        r = p.wait()
        if (r != 0):
            return False
        else:
            # parsovat
            flist = set()
            with open("cath-b-newest-all") as fd:
                for line in fd.readlines():
                    s = line.split(" ")
                    flist.add(s[2])
            os.remove("cath-b-newest-all")
            return flist
    else:
        print(f"ERROR: problem in getting list - address {link}")
        return set()


def read_fdir_family_list(fdir):
    ret = []
    for famdir in os.listdir(fdir):
        fs = famdir.split("-")
        fam = fs[1]
        ret.append(fam)
    return ret


def error_report(directory, command, out, err, log_fn):
    # unified error report to log_fn file
    with open(log_fn, "w") as fd:
        timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
        fd.write(f"* date: {timestamp}\n")
        fd.write(f"* commandline: {command}\n")
        fd.write(f"* dir: {directory}\n")
        fd.write(f"* stdout:\n")
        if out is not None:
            if type(out) == bytes:
                fd.write(out.decode("utf-8"))
            else:
                fd.write(out)
        else:
            fd.write("None\n")
        fd.write(f"* stderr:\n")
        if err is not None:
            if type(err) == bytes:
                fd.write(err.decode("utf-8"))
            else:
                fd.write(str(err))
        else:
            fd.write("None\n")


def change_comparison(n_fn, o_fn):
    return not filecmp.cmp(n_fn, o_fn)


def color_sses_comparison(err, n_data, o_data):
    # test there is compared amount of prime sses in new and old family diagram file
    n_sses = 0
    if "nodes" in n_data:
        for sses in n_data["nodes"]:
            if sses["occurrence"] > 0.8:
                n_sses += 1

    o_sses = 0
    if "nodes" in o_data:
        for sses in o_data["nodes"]:
            if sses["occurrence"] > 0.8:
                o_sses += 1

    if not n_sses == o_sses:
        err.append(f"family changed prime sses number from {o_sses} to {n_sses}")
        print(f"family changed prime sses number from {o_sses} to {n_sses}")
    return err


def all_sses_comparison(err, n_data, o_data):
    # test there is compared amount of prime sses in new and old family diagram file
    n_sses = 0
    if "nodes" in n_data:
        for sses in n_data["nodes"]:
            n_sses += sses["occurrence"]

    o_sses = 0
    if "nodes" in o_data:
        for sses in o_data["nodes"]:
            o_sses += sses["occurrence"]

    if n_sses < 0.95*o_sses:
        err.append(f"family downgrade all sses number from {o_sses} to {n_sses}")
        print(f"family downgrade all sses number from {o_sses} to {n_sses}")
    if o_sses < 0.95*n_sses:
        err.append(f"family upgrade all sses number from {o_sses} to {n_sses}")
        print(f"family upgrade all sses number from {o_sses} to {n_sses}")
    return err


def read_p_sses(filename, counter):
    with open(filename) as fd:
        data = json.load(fd)
        for i in data["h_type"]:
            if data["h_type"][i] == 1 and i not in counter:
                counter.append(i)
    return counter


def rotation_file_test(flist):
    counter_notfound = 0
    counter_emptyfamily = 0
    for fam in flist:
        link = f"https://2dprots.ncbr.muni.cz/files/family/{fam}/rotation.json"
        ret = curl_file(link, "", "/dev/null", err_report=False, option="-I")
        if ret is not True:
            # we want to find whether the family is not empty
            try:
                loader.Family(fam).domains
            except two_d_prot.pdbe_api.PDBeApiReturnCodeError as e:
                if int(str(e)) == 404:
                    # there is no domains in the family - we have to create empty domain list
                    print("empty family", fam)
                    counter_emptyfamily += 1
                    continue
            print("nonempty family with no rotation file", fam)
            counter_notfound += 1
    print(f"empty family {counter_emptyfamily} times")
    return counter_notfound


# Read input arguments
parser = argparse.ArgumentParser(
    description="Test input family dbs directory, check wheather its layout contains the same or higner",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument("--fdir", required=False, action="store",
                    help="go through all families which has generated directory in this directory")
parser.add_argument("--all", action="store_true",
                    help="go through all cath families")
parser.add_argument("--rotation", required=False, action="store_true",
                    help="test the presence of rotation file in nonempty family estimated time 30 min")
args = parser.parse_args()

flist = []
if args.fdir:
    flist = read_fdir_family_list(args.fdir)

if args.all:
    flist = read_cath_family_list()

if args.rotation:
    rotation_file_test(flist)
