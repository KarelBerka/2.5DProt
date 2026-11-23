#!/usr/bin/env python3

import time
import argparse
import os
import subprocess


test_set = [{"name": "1.10.10.140",
             "is": ["5w97h00"],
             "is_not": ["5w97U00"],
             },
            {"name": "1.10.10.1020",
             "is": ["1w36D01", "1w36G01", "3k70D01", "3k70G01", "5ld2D01"],
             "is_not": ["1w36C01", "1w36F01", "3k70C01", "3k70F01", "5ld2C01"],
             },
            {"name": "6.20.360.10",
             "is": ["4dx9q00", "4dx9w00"],
             "is_not": ["4dx9S00", "4dx9EB00"],
             }]


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


def test_dir_family(gdir, family, dis, disnot):
    # for each family instance find whether this domain svg are in the given
    # family dbs and disnot are not there

    errors = {"family": 0, "domains": 0}
    ok = {"family": 1, "domains": 0}

    family_found = False
    for fdir in os.listdir(gdir):
        if fdir.startswith(f"generated-{family}-"):
            family_found = True
            # go through the list of domains which should be in the family
            # test the presence of relevant *.svg file. e.g. 4dx9_q00.svg
            # in 6.20.360.10
            for issvg in dis:
                present = False
                for ffile in os.listdir(f"{gdir}/{fdir}"):
                    if ffile == f"{issvg[0:4]}_{issvg[4:]}.svg":
                        present = True
                        ok["domains"] += 1
                        break
                if present is False:
                    print(f"BUG: family {family} in dir {gdir}/{fdir} should contain {issvg} ({issvg[0-4]}_{issvg[4:]}.svg) but it is not there")
                    errors["domains"] += 1

            # go through the list of domains which should not be in the family
            # test there is not the relevant *.svg file. e.g. 4dx9_S00.svg
            # in 6.20.360.10
            for isnotsvg in disnot:
                present = False
                for ffile in os.listdir(f"{gdir}/{fdir}"):
                    if ffile == f"{isnotsvg[0:4]}_{isnotsvg[4:]}.svg":
                        present = True
                        print(f"BUG: family {family} in dir {gdir}/{fdir} should not contain {isnotsvg} ({isnotsvg[0-4]}_{isnotsvg[4:]}.svg) but it is there")
                        errors["domains"] += 1
                        break
                if present is False:
                    ok["domains"] += 1

    if family_found is False:
        errors["family"] = 1
        ok["family"] = 0

    return (ok, errors)


def test_web_family(gdir, family, dis, disnot):
    # test the presence of domains in domain_list.txt file of the given family
    errors = {"family": 0, "domains": 0}
    ok = {"family": 1, "domains": 0}
    link = f"https://2dprots.ncbr.muni.cz/files/family/{family}/domain_list.txt"
    filename = f"{family}_domain_list.txt"
    ret = curl_file(link, "", filename, True, "-L")
    if ret is not True:
        # file was not found family error
        ok["family"] = 0
        errors["family"] = 1
        return (ok, errors)

    with open(filename, "r") as fd:
        listis = {}
        listisnot = {}
        for i in dis:
            listis[i] = 0
        for i in disnot:
            listisnot[i] = 0
        pcd_line = fd.readline()
        while pcd_line:
            sp = pcd_line.split("'")
            if len(sp) > 8:
                pcd = sp[7]
            else:
                pcd = f"{sp[1]}{sp[3]}{sp[5]}"
            for i in listis:
                if i == pcd:
                    listis[i] = 1
            for i in listisnot:
                if i == pcd:
                    listisnot[i] = 1
            pcd_line = fd.readline()
    for i in listis:
        if listis[i] == 1:
            ok["domains"] += 1
        else:
            errors["domains"] += 1
    for i in listisnot:
        if listisnot[i] == 1:
            errors["domains"] += 1
        else:
            ok["domains"] += 1

    os.remove(filename)
    return (ok, errors)


def test_chain_names(gdir, web):
    # go through whole tes_set and run test for each family instance
    errors = {"family": 0, "domains": 0}
    ok = {"family": 0, "domains": 0}

    for test in test_set:
        if web is not False:
            res = test_web_family(web, test["name"], test["is"], test["is_not"])
        else:
            res = test_dir_family(gdir, test["name"], test["is"], test["is_not"])
        ok["family"] += res[0]["family"]
        ok["domains"] += res[0]["domains"]
        errors["family"] += res[1]["family"]
        errors["domains"] += res[1]["domains"]
    if web is not False:
        print("2DProt web chain names test: ok", ok, "errors", errors)
    else:
        print(gdir, "directory dbs chain names test: ok", ok, "errors", errors)


# Read input arguments
parser = argparse.ArgumentParser(
    description="Test whether the chains in family proteins are named in the right way",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument("--dir", required=False, action="store",
                    help="test the dbs in input directory")
parser.add_argument("--web", action="store_true", help="test the dbs in 2DProt web")
args = parser.parse_args()

flist = []
if args.dir:
    test_chain_names(args.dir, False)

if args.web:
    test_chain_names(None, args.web)
