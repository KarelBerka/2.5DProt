#!/usr/bin/env python3
import os
import shutil
import datetime
import argparse
import sys


def move_prot(family, dir_name):
    # move protein directory to given folder
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    i = family.rfind("/")

    if i != -1:
        gener_dir_name = f"{dir_name}/{family[i:]}"
    else:
        gener_dir_name = f"{dir_name}/{family}"

    if not os.path.isdir(gener_dir_name):
        shutil.move(family, gener_dir_name)
    else:
        print(f"directory {gener_dir_name} already existes")
        shutil.rmtree(family)


def copy_new_prot_to_dbs(input_dir, prot_dir, output_dir):
    # copy new family generation to dbs directory

    # ommit unnecessary files
    ds = prot_dir.split("/")
    s = ds[-1].split("-")
    full_prot_dir = f"{input_dir}/{prot_dir}"
    log_file = f"{full_prot_dir}/main-{s[1]}"
    if os.path.isfile(log_file):
        os.remove(log_file)

    annotation_dir = full_prot_dir + "/" + "annotation"
    if os.path.isdir(annotation_dir):
        shutil.rmtree(annotation_dir)

    cif_dir = full_prot_dir + "/" + "cif"
    if os.path.isdir(cif_dir):
        shutil.rmtree(cif_dir)

    # move the rest
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(f"{output_dir}/{ds[0]}"):
        os.mkdir(f"{output_dir}/{ds[0]}")
    if not os.path.isdir(f"{output_dir}/{ds[0]}/{ds[1]}"):
        os.mkdir(f"{output_dir}/{ds[0]}/{ds[1]}")
    # test whether there already are not any other version of this protein, prefer the newer one

    shutil.move(full_prot_dir, f"{output_dir}/{ds[0]}/{ds[1]}")


def sort_protein_directory(missing_l_f, prot_fail, prot_finished, protein_dir, input_dir, move_prot_f, problematic_fail, updated, added):
    if missing_l_f == 1:
        # there is not relevant ligand file
        move_prot(protein_dir, f"{output_dir}/missing_ligand_file")
        move_prot_f += 1
    elif prot_fail > 0:
        # there was a fail
        # TODO test whether the fail was not present in the old version as well
        move_prot(protein_dir, f"{input_dir}/failed")
        problematic_fail += 1
    elif not prot_finished:
        # there is no fail, but the run was not finished
        print(f"move_prot(protein_dir, f'{input_dir}/not_finished')")
        move_prot(protein_dir, f"{input_dir}/not_finished")
    else:
        # the protein was generated properly
        # move the old one to timestamp directory
        split_inner_dir = inner_dir3.split("-")
        prefix_old_prot = split_inner_dir[0] + "-" + split_inner_dir[1] + "-"
        dbs_prot = []
        output_dir_inner = f"{output_dir}/{inner_dir1}/{inner_dir2}"
        if os.path.isdir(output_dir_inner):
            for filename in os.listdir(output_dir_inner):
                print(prefix_old_prot, filename)
                if filename.startswith(prefix_old_prot):
                    dbs_prot.append(output_dir_inner + "/" + filename)

        # TODO: test whether there is not several occurances of protein diagram dir
        if len(dbs_prot) != 0:
            for i in dbs_prot:
                shutil.move(i, path_ts_dir)
            updated += 1
        else:
            print("No older version of", prefix_old_prot)
            added += 1

        # move the new one to database directory
        # remove unneccessary files as well
        copy_new_prot_to_dbs(input_dir, f"{inner_dir1}/{inner_dir2}/{inner_dir3}", output_dir)


parser = argparse.ArgumentParser(
    description="filter properly finished potein diagrams which are already without solvable problems")
parser.add_argument("--output", required=True, action="store",
                    help="directory to which will be the input proteins parsed")
parser.add_argument("--input", required=True, action="store",
                    help="directory to which contains proteins to be parsed")
args = parser.parse_args()

input_dir = args.input
if not os.path.isdir(input_dir):
    print(f"ERROR: Input directory {input_dir} does not exist, nothing parsed")
    sys.exit()

output_dir = args.output
if not os.path.isdir(output_dir):
    print(f"INFO: No output directory {input_dir}, the directory was created")
    os.mkdir(output_dir)

not_finished = 0
problematic_fail = 0
move_prot_f = 0

updated = 0
added = 0
all_proteins = 0

today = datetime.date.today()
ts_dir = today.strftime("old-%Y_%m_%d")
path_ts_dir = input_dir+"/"+ts_dir
if not os.path.isdir(path_ts_dir):
    os.mkdir(path_ts_dir)

# protein directory structure has three levels of directories:
# 1/t/generated-1tqn-timestamp/
for inner_dir1 in os.listdir(input_dir):
    if not os.path.isdir(f"{input_dir}/{inner_dir1}") or len(inner_dir1) != 1:
        continue
    for inner_dir2 in os.listdir(f"{input_dir}/{inner_dir1}"):
        if not os.path.isdir(f"{input_dir}/{inner_dir1}/{inner_dir2}") or len(inner_dir1) != 1:
            continue
        for inner_dir3 in os.listdir(f"{input_dir}/{inner_dir1}/{inner_dir2}"):
            if len(inner_dir3) < 15:
                # bogus file name format
                continue

            # we are in per protein generated directory
            # e.g. 1/t/generated-1tqn-timestamp/...
            protein_dir = f"{input_dir}/{inner_dir1}/{inner_dir2}/{inner_dir3}"
            protein = inner_dir3[10:14]
            prot_finished = True

            # test whether the generating was finished
            # thus contains {protein}.svg and {protein}_l.svg images
            if not os.path.isfile(f"{protein_dir}/{protein}.svg") or \
               not os.path.isfile(f"{protein_dir}/{protein}_l.svg"):
                prot_finished = False
                not_finished += 1

            prot_fail = 0
            missing_l_f = 0
            if os.path.isdir(protein_dir):
                for filename in os.listdir(protein_dir):
                    # test whether there is no fail
#                   if filename == "missing_ligand_file.txt":
#                       # if there are missing ligand files, they are filtered to the separate directory
#                       missing_l_f = 1
#                       break
                    if filename.find("-fail-") != -1 and \
                        filename != "not_in_mapping_fail.txt" and \
                        filename != "empty protein_fail.txt" and \
                        not filename.startswith("secstrcutter-fail") and \
                        not filename.startswith("curl-fail"):
                        prot_fail += 1

            # sort prodein dir based on set flags
            sort_protein_directory(missing_l_f, prot_fail, prot_finished,
                                   protein_dir, input_dir,
                                   move_prot_f, problematic_fail, updated, added)
            all_proteins += 1

print("Missing ligand file", move_prot_f)
print("Not finished: ", not_finished - move_prot_f - problematic_fail)
print("Failed families: ", problematic_fail)
print("Updated: ", updated)
print("Added:", added)
print("All proteins: ", all_proteins)
