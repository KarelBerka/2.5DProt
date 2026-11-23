#!/usr/bin/env python3
import json
import sys
import math


def file_stats(file_name):
    # output the variance statistics for the input file
    # TODO write here the possible output
    with (open(file_name, "r")) as fd:
        data = json.load(fd)
    res = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for sses in data["consensus"]["secondary_structure_elements"]:
        # we calculate the sses confidence as the minimum of confidence of all sses res
        minconf = (15 - math.sqrt(sses["variance"]))*7

        # based of average sses confidence we use the same color as alfa fold
        if minconf > 90:
            res[1] += 1
            res[6] += sses["occurrence"]
        elif minconf > 70:
            res[2] += 1
            res[7] += sses["occurrence"]
        elif minconf > 50:
            res[3] += 1
            res[8] += sses["occurrence"]
        else:
            res[4] += 1
            res[9] += sses["occurrence"]
        res[0] += 1
        res[5] += sses["occurrence"]
        res[10] += sses["occurrence"]*minconf

    res[1] = int(100*res[1]/res[0])
    res[2] = int(100*res[2]/res[0])
    res[3] = int(100*res[3]/res[0])
    res[4] = int(100*res[4]/res[0])
    res[9] = int(100*res[9]/res[5])
    res[6] = int(100*res[6]/res[5])
    res[7] = int(100*res[7]/res[5])
    res[8] = int(100*res[8]/res[5])
    res[10] = int(res[10]/res[5])
    return res


number = 1
stats = []
while number < len(sys.argv):
    file_name = sys.argv[number]
    stats.append(file_stats(file_name))
    number += 1


avg = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

rec = 0
for i in stats:
    num = 0
    while num < 11:
        avg[num] += i[num]
        num += 1
    rec += 1

num = 0

while num < 10:
    avg[num] = avg[num]/rec
    num += 1

print("sses_number    in 90      in 70    in 50    less then 50 ")

print("-----------")
print(avg)
