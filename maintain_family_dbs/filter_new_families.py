#!/usr/bin/env python3
import sys
import os
import shutil
import datetime
import argparse


def move_family(family, dir_name):
    prefix_dir = family[:family.rfind("/")]
    path_directory = prefix_dir + "/" + dir_name
    if not os.path.isdir(path_directory):
        os.mkdir(path_directory)
    shutil.move(family, path_directory)


def copy_new_fam_to_dbs(family_dir, old_dir):
    # copy new family generation to dbs directory
    # TODO if there already is a generated dir of this family, leave there only the newer one
    # ommit unnecessary files
    dir_structure = family_dir.split("/")
    s = dir_structure[-1].split("-")

    dir_1 = family_dir + "/" + s[1]
    if os.path.isdir(dir_1):
        # metacentrum remove this directory itself
        for ffile in os.listdir(dir_1):
            if (ffile.endswith("-detected.sses.json")) or \
               (ffile.endswith(".cif") and not ffile.startswith("consensus")):
                os.remove(dir_1 + "/" + ffile)

    if os.path.isdir(dir_1):
        # metacentrum remove this directory itself
        dir_2 = family_dir + "/" + s[1]+"/consensus_A00"
        for ffile in os.listdir(dir_2):
            if (ffile.endswith("-detected.sses.json")) or \
               (ffile.endswith("-annotated.sses.json")) or \
               (ffile.endswith("-alignment.json")) or \
               (ffile.endswith("cif") and not ffile.startswith("consensus")):
                os.remove(dir_2 + "/" + ffile)

    # move the rest
    shutil.move(family_dir, old_dir)


parser = argparse.ArgumentParser(
    description="filter properly finished family diagrams which are already without solvable problems")
parser.add_argument("--output", required=True, action="store",
                    help="directory to which will be the input family diagram directories parsed")
parser.add_argument("--input", required=True, action="store",
                    help="directory to which contains families to be parsed")
args = parser.parse_args()

new_dir = args.input
if not os.path.isdir(new_dir):
    print(f"ERROR: Input directory {new_dir} does not exist, nothing parsed")
    sys.exit()

dbs_dir = args.output
if not os.path.isdir(dbs_dir):
    print("INFO: No output directory {dbs_dir}, the directory was created")
    os.mkdir(dbs_dir)

not_finished = 0
fail = 0
updated = 0
added = 0

today = datetime.date.today()
ts_dir = today.strftime("old-%Y_%m_%d")
path_ts_dir = new_dir+"/"+ts_dir
if not os.path.isdir(path_ts_dir):
    os.mkdir(path_ts_dir)

# go throught the directory structure to which group
# (not finished | failed | updated | added) they belonged
for inner_dir in os.listdir(new_dir):
    if inner_dir.startswith("generated"):
        family_dir = new_dir+"/"+inner_dir
        fam_fail = 0

        # test whether the generating was finished
        # (domain_list.txt file was generated)
        f_dir_split = family_dir.split("-")
        family = f_dir_split[1]
        if not os.path.isfile(family_dir+"/domain_list.txt"):
            fam_finished = False
            not_finished += 1
        elif os.path.isfile(f"{family_dir}/main-{family}.log"):
            end_is_there = False
            with open(f"{family_dir}/main-{family}.log","r") as fd:
                line = fd.readline()
                while line:
                    if line.startswith("--- Timestamp: 6 OK"):
                        end_is_there = True
                    line = fd.readline()
            if end_is_there:
                fam_finished = True
            else:
                fam_finished = False
                not_finished += 1
        else:
            fam_finished = False
            not_finished += 1

        for filename in os.listdir(family_dir):
            # test whether there is no fail
            if filename.find("-fail-") != -1:
                fam_fail += 1

        if fam_fail > 0:
            # there was a fail
            move_family(family_dir, "failed")
            fail += 1
        elif not fam_finished:
            # there is no fail, but the run was not finished
            move_family(family_dir, "not_finished")
        else:
            # the family was generated properly
            # move the old one(s) to new_dbs directory
            # if there is not the same generated directory,
            # otherwise just remove it (it is already there)
            if not os.path.isdir(f"{dbs_dir}/{os.path.basename(family_dir)}"):
                copy_new_fam_to_dbs(family_dir, dbs_dir)
            else:
                shutil.rmtree(family_dir)
                print(f"WARNING: {dbs_dir}/{os.path.basename(family_dir)} already exists in old directory")
            # find out wheather there are not duplicares of the moved family
            # if "yes" just remove the old one
            split_inner_dir = inner_dir.split("-")
            prefix_old_fam = split_inner_dir[0] + "-" + split_inner_dir[1] + "-"
            dbs_fam = []
            for filename in os.listdir(dbs_dir):
                if filename.startswith(prefix_old_fam):
                    dbs_fam.append(dbs_dir + "/" + filename)

            if len(dbs_fam) > 1:
                # there is more then one family generated dir
                # leave only the newest one
                leave = dbs_fam[0]
                for d in dbs_fam:
                    if leave < d:
                        leave = d
                for d in dbs_fam:
                    if d != leave:
                        if not os.path.isdir(f"{path_ts_dir}/{os.path.basename(d)}"):
                            copy_new_fam_to_dbs(d, path_ts_dir)
                        else:
                            shutil.rmtree(d)
                            print(f"WARNING: {dbs_dir}/{os.path.basename(d)} already exists in old directory")
                updated += 1
            else:
                print("No older version of", prefix_old_fam)
                added += 1
            # move the new one to database directory
            # remove unneccessary files as well

print("Not finished: ", not_finished)
print("Failed families: ", fail)
print("Updated: ", updated)
print("Added: ", added)
print("All families:", not_finished + fail + updated + added)
