#!/usr/bin/env python3
# generate family list, which are in CATH dbs
# (and eventually are not generated in the given directory)
import argparse
import sys
import os
import time
import logging
import requests

sys.path.append("../api")
import utils_cath
sys.path.append("..")
import utils_common
import utils_pdbe_api


def family_in_directory(fam, fdir):
    # test wheather the famili is generated in the given directory
    fstring = f"generated-{fam}-"
    for gdir in os.listdir(fdir):
        if gdir.startswith(fstring):
            return True
    return False


def family_in_blacklist(fam, blacklist):
    # test wheather the famili is generated in the blacklist
    # blacklist format is:
    # <...>-<familyname>-<...>
    # e.g.
    # generated-1.10.10.10-2024-02-24T20_02_23
    fnstring = f"-{fam}-"
    with open(blacklist, "r") as fb:
        line = fb.readline()
        while(line):
            if fnstring in line:
                print(f"{fam} skipped")
                return True
            line = fb.readline()

    return False

def equiv(l1,l2):
    for i in l1:
        if i not in l2:
            return False
    for i in l2:
        if i not in l1:
            return False
    return True


def family_on_web(link, l_dir, fam, new_domain_list):
    # test the presence of family on 2DProts web
    try: 
        response = requests.head(link)
        if (response.status_code == 200) or (response.status_code == 302):
            print(f"T fam {fam} {link}")
            return True
#          try:
#            family_domain_list = f"{l_dir}/{fam}_domain_list.txt"
#            ndl = utils_pdbe_api.get_family_domains(fam)
#            ndl2 = []
#            for i in ndl:
#                ndl2.append(f"{i[0]}{i[1]}{i[2]}")
#
#            # download family domain list from 2dprot web
#            utils_common.get_url_to_file(link, f"{l_dir}/{fam}_domain_list.txt")
#            with open(family_domain_list, "r") as fd:
#              # read family domain list
#              line = fd.readline()
#              domain_list=[]
#              while (line):
#                  sline = line.split("'")
#                  domain_list.append(f"{sline[1]}{sline[3]}{sline[5]}")
#                  line = fd.readline()
#              # test the difference in actual and downloaded list
#              if not equiv(domain_list,ndl2):
#                  print(len(ndl2), fam, len(domain_list))
#              else:
#                  print("OK", fam)
#            return True
#          except:
#            print(len(new_domain_list[fam]), fam, 0)
        else:
            print(f"F fam {fam} {response.status_code} - {link}")
            return False
    except:
        return False


parser = argparse.ArgumentParser(description='Create 2DProt family dbs update list')
parser.add_argument('--all', action='store_true',
                    help='list all families, which are in cathdb')
parser.add_argument('--dir', action='store',
                    help='list all families, which are in cathdb, and have not generated-* file in <dir> directory')
parser.add_argument('--twodprots', action='store_true',
                    help='list all families, which are in cathdb, and have not record in 2DProts')
parser.add_argument('--dev', action='store_true',
                    help='update all families, which are in cathdb, and have not record in development version of 2DProts')
parser.add_argument('--debug', action='store_true',
                    help='create debug output')
parser.add_argument('--cath_dbs', action='store',
                    help='cath-b-newest-latest-release file')
parser.add_argument('--blacklist', action='store',
                    help='blacklist file - families, which could not be on list - have some major problem')
args = parser.parse_args()


# create new family diagram database directory
timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
start_time = time.time()
last_time = start_time
new_dir = "generated-" + timestamp
os.mkdir(new_dir)

if args.debug:
    logging.basicConfig(level=logging.DEBUG)


# ommit calls when more then 1 argument is used
# and set dbs comparison name
name = ""
if args.all:
    if args.dir or args.twodprots or args.dev:
        print("Use only one option")
    else:
        name = "all"
if args.dir:
    if args.twodprots or args.dev:
        print("Use only one option")
    else:
        name = "dir"
if args.twodprots:
    if args.dev:
        print("Use only one option")
    else:
        name = "2dprots"
if args.dev:
    name = "dev"

list_fn = f"{new_dir}/list-{name}-{timestamp}.txt"

if (args.all) or (args.dir) or (args.twodprots) or (args.dev):
    # find out all families in the current CATH dbs
    if args.cath_dbs:
        f_list = utils_cath.cath_family_list(new_dir, True, args.cath_dbs)
    else:
        f_list = utils_cath.cath_family_list(new_dir, True)
    new_domain_list = {}
    # download actual family list and list of domains
    with open(f"{new_dir}/cath-b-newest-latest-release") as fd:
        line = fd.readline()
        while line:
            sline = line.split(" ")
            if sline[2] not in new_domain_list:
                new_domain_list[sline[2]] = []
            new_domain_list[sline[2]].append(sline[0])
            line = fd.readline()

    with open(list_fn, "w") as fd:
        families = 0
        for fam in f_list:
            if ("." in fam):
                if (args.all):
                    fd.write(f"{fam}\n")
                    families += 1
                if (args.dir) and not family_in_directory(fam, args.dir) and (not args.blacklist or not family_in_blacklist(fam, args.blacklist)):
                    fd.write(f"{fam}\n")
                    families += 1
                if (args.twodprots) and not family_on_web(f"https://2dprots.ncbr.muni.cz/files/family/{fam}/domain_list.txt", new_dir, fam, new_domain_list):
                    fd.write(f"{fam}\n")
                    families += 1
#                if (args.dev) and not family_on_web(f"http://147.251.21.23//files/family/{fam}/domain_list.txt"):
                if (args.dev) and not family_on_web(f"http://2dprots-stage.ncbr.muni.cz:5173//files/family/{fam}/domain_list.txt", new_dir, fam, new_domain_list):
                    fd.write(f"{fam}\n")
                    families += 1
else:
    print("no option set")

print(f"List filename is {list_fn}, contains {families} families.")
