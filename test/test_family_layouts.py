#!/usr/bin/env python3

import time
import argparse
import os
import subprocess
import json
import filecmp
import sys

def error_report(directory, command,  out, err, log_fn):
    # unified error report to log_fn file
    with open(log_fn, "w") as fd:
         timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
         fd.write(f"* date: {timestamp}\n")
         fd.write(f"* commandline: {command}\n")
         fd.write(f"* dir: {directory}\n")
         fd.write(f"* stdout:\n")
         if out != None:
             if type(out) == bytes:
                 fd.write(out.decode("utf-8"))
             else:
                 fd.write(out)
         else:
             fd.write("None\n")
         fd.write(f"* stderr:\n")
         if err != None:
             if type(err) == bytes:
                 fd.write(err.decode("utf-8"))
             else:
                 fd.write(str(err))
         else:
             fd.write("None\n")


def curl_file(link, directory , filename, err_report = True, option = ""):
    # curl "link" file to given "filename"
    command = f"curl {option} -f {link} -o {filename} >/dev/null 2>/dev/null"
    p = subprocess.Popen(command, shell = True, stdout=subprocess.PIPE)
    r = p.wait()
 
    if (r != 0):
        if (err_report == True):
            out, err = p.communicate()
            log_fn = f"{directory}/curl-fail-{time.time()}.txt"
            error_report(directory, command, out, err, log_fn)
            print(f"curl of {link} failed, report save to {log_fn}")
        return False
    return True


def change_comparison(n_fn, o_fn):
    return not filecmp.cmp(n_fn,o_fn)


def color_sses_comparison(err, n_data, o_data):
    #test there is compared amount of prime sses in new and old family diagram file
    n_sses = 0
    if "nodes" in n_data:
        for sses in n_data["nodes"]:
            if sses["occurrence"] >0.8:
                n_sses += 1

    o_sses = 0
    if "nodes" in o_data:
        for sses in o_data["nodes"]:
            if sses["occurrence"] >0.8:
                o_sses += 1

    if not n_sses == o_sses:
        err.append(f"family {line[:-1]} changed prime sses number from {o_sses} to {n_sses}")
        print(f"family {line[:-1]} changed prime sses number from {o_sses} to {n_sses}")
    return err


def all_sses_comparison(err, n_data, o_data):
    #test there is compared amount of prime sses in new and old family diagram file
    n_sses = 0
    if "nodes" in n_data:
        for sses in n_data["nodes"]:
            n_sses += sses["occurrence"]

    o_sses = 0
    if "nodes" in o_data:
        for sses in o_data["nodes"]:
            o_sses += sses["occurrence"]

    if n_sses < 0.95*o_sses:
        err.append(f"family {line[:-1]} downgrade all sses number from {o_sses} to {n_sses}")
        print(f"family {line[:-1]} downgrade all sses number from {o_sses} to {n_sses}")
    if o_sses < 0.95*n_sses:
        err.append(f"family {line[:-1]} upgrade all sses number from {o_sses} to {n_sses}")
        print(f"family {line[:-1]} upgrade all sses number from {o_sses} to {n_sses}")
    return err


def read_p_sses(filename, counter):
    with open(filename) as fd:
        data = json.load(fd)
        for i in data["h_type"]:
            if data["h_type"][i] == 1 and i not in counter:
               counter.append(i)
    return counter

# Read input arguments
parser = argparse.ArgumentParser(
    description="Test input family dbs directory, check wheather its layout contains the same or higner",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--fdir", required=True, action="store",
                    help="location of the families to be tested against currenr 2DProt dbs on web (TODO: Link)")
args = parser.parse_args()

# Read each family
data = {}
for famdir in os.listdir(args.fdir):
    counter = 0
    fs = famdir.split("-")
    fam = fs[1]
    webpsses = []
    dirpsses = []
    for ffile in os.listdir(f"{args.fdir}/{famdir}"):
        if ffile.startswith("layout-") and ffile.endswith(".json"):
            data[famdir] = ffile
            webfile = f"https://2dprots.ncbr.muni.cz/static/web/layouts/generated-{fam}/{ffile}"
            ret = curl_file(webfile, ffile, f"aux.json", err_report = False)
            if ret:
                webpsses = read_p_sses("aux.json", webpsses)
                dirpsses = read_p_sses(f"{args.fdir}/{famdir}/{ffile}", dirpsses)
            else:
                print("--", fam, ffile, "not in web") 
            counter += 1
            if counter > 5:
                break
    if len(webpsses) > len(dirpsses):
        print("--", len(webpsses)- len(dirpsses), fam, ffile, len(webpsses), "<=", len(dirpsses), "error error error", webpsses, dirpsses)
    else:
        print("--", len(webpsses)- len(dirpsses), fam, ffile, len(webpsses), "<=", len(dirpsses))


