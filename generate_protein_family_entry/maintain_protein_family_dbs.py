#!/usr/bin/env python3

import argparse

def omit_characters(sline, characters):
    res = []
    for i in sline:
        while len(i) > 0 and i[0] in characters:
            i = i[1:]
        while len(i) > 0 and i[-1] in characters:
            i = i[:-1]
        if len(i)> 0:
            res.append(i)
    return(res)

parser = argparse.ArgumentParser(
    description='todo')
parser.add_argument('--dbs_file', required=True, action="store", help="todo")
parser.add_argument('--min_fam_size', action="store", help="minimal number of domains in family")
parser.add_argument('--homog', action="store_true", help="only protein families with the same set of famili 'domais'")
parser.add_argument('--output_dir', action="store", help="only protein families with the same set of famili 'domais'")

args = parser.parse_args()


with open(args.dbs_file, "r") as fd:
    line = fd.readline()
    fam_clusters = {}
    pf_sum = {}
    while line:
        sline = line[:].split(" ")
        sline = omit_characters(sline, "[',]\n")
#        print(sline)
        number = sline[0]
        line = fd.readline()
        fam_set = []
        prot_set = []
        for i in sline[1:]:
            if "." in i:
                fam_set.append(i)
            else:
                prot_set.append(i)
        if len(fam_set) == 0:
            # proteins does not have domains in cath
            continue
        homog = True
        if len(fam_set) != 0:
            for i in fam_set:
                if i != fam_set[0]:
                    homog = False
                    break

        if fam_set[0] not in pf_sum:
            pf_sum[fam_set[0]] = int(number)
        else:
            pf_sum[fam_set[0]] += int(number)


        if args.min_fam_size and int(number) < int(args.min_fam_size):
            continue
        if not homog and args.homog:
            continue


        if homog:
            if fam_set[0] not in fam_clusters:
                fam_clusters[fam_set[0]] = {len(fam_set): (number, prot_set)}
            else:
                fam_clusters[fam_set[0]][len(fam_set)] = (number, prot_set)

for i in fam_clusters:
#    print(i , pf_sum[i]/4, end=" ")
    for j in fam_clusters[i]:
       if pf_sum[i] / 4 < int(fam_clusters[i][j][0]):
          with open(f"{args.output_dir}/{j}_{i}.txt","w") as fd:
              for p in fam_clusters[i][j][1]:
                  fd.write(f"{p}\n")
#          print(f"{i} ({j}, {fam_clusters[i][j]})")
#       else:
#          print(f"[{j}, {fam_clusters[i][j]}]", end = " ")

#          print("pass - limit", pf_sum[i]/5, fam_clusters[i][j])
    print()

#for i in fam_clusters:
#    print("i",i, pf_sum[i], fam_clusters[i])