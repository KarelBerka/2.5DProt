#!/usr/bin/env python3

# This script go through database and do statistics regarding familie members change in time
# input file list containing statistics (EXAMPLE: update_2021-08-28T10:10:10,0000000000+02:00)

import sys
import logging

def read_stats(stat_file, stat_data):
   # read whole stat file and add domain x chains to each family
   with open(stat_file, "r") as fd:
      stat_data[stat_file]={}
      c=fd.readline()
      while (c != ""):
          # go through the file until you reach tehe end of file
          if (c.find(" all: ") == -1):
              # read only non-summary rows
              s = c.split()
              if s[0] not in stat_data[stat_file]:
                  stat_data[stat_file][s[0]]=set([s[2]])
              else:
                  stat_data[stat_file][s[0]].add(s[2])
          c = fd.readline()

def week_stats(stat_data):
    dates = list(stat_data.keys())
    print("=== per update file statistics ==================================")
    for i in range(1,len(dates)):
        # for each pair of statistic count the diffs
        st1, st2 = dates[i-1], dates[i]
        diffc = 0

        addc = 0
        addcd = 0

        remc = 0
        remcd = 0

        for j in stat_data[st1].keys():
            # for each family in the time snapshot:
            if (j not in stat_data[st2]):
               # se should have the family in the next snapshot too
               stat_data[st2][j]=set()
#            print(len(stat_data[st2][j]), stat_data[st2][j])
            if (len(stat_data[st2][j]) != 0) and (len(stat_data[st1][j]) != 0):
                if stat_data[st1][j] != stat_data[st2][j]:
                    diffc += 1
                if stat_data[st1][j].difference(stat_data[st2][j]) != set():
                    remc += 1
                    remcd += len(stat_data[st1][j].difference(stat_data[st2][j]))
                if stat_data[st2][j].difference(stat_data[st1][j]) != set():
                    addc += 1
                    addcd += len(stat_data[st2][j].difference(stat_data[st1][j]))
        print(f"{st1} {st2} {i} data, celkem zmeneno familii {remcd+addcd}, celkem odebrano cd  ({remcd}), celkem prideno dc ({addcd})" )
    print("=============================================================================")

def family_stats(stat_data):
    counter_all = {}
    dates = list(stat_data.keys())
    print(dates)

    all_families = set()
    for snapshot in list(stat_data.keys()):
       for fam in stat_data[snapshot].keys():
           if fam not in all_families:
              all_families.add(fam)

    # all changed domain x chains either added or removed
    for fam in all_families:
        counter = {}
        for i in range(1, len(dates)):
            snapshot1, snapshot2 = dates[i-1], dates[i]
            if fam not in stat_data[snapshot1]:
                stat_data[snapshot1][fam] = set()
            if fam not in stat_data[snapshot2]:
                stat_data[snapshot2][fam] = set()
            dom_chain_dif_list = stat_data[snapshot1][fam].symmetric_difference(stat_data[snapshot2][fam])
            if (len(stat_data[snapshot2][fam]) != 0) and (len(stat_data[snapshot1][fam]) != 0):
                for dom_chain in dom_chain_dif_list:
                   if (dom_chain,fam) not in counter:
                      counter[(dom_chain,fam)] = 1
                   else:
                      counter[(dom_chain,fam)] += 1
                   counter_all[(dom_chain, fam)]=counter[(dom_chain, fam)]

    val = set(sorted(list(counter_all.values())))
    number_fails={}
    for i in counter_all.keys():
       if counter_all[i] not in number_fails:
           number_fails[counter_all[i]] = 1
       else:
           number_fails[counter_all[i]] += 1
    print("============================================================")
    print("changed domain x chains (both added, removed", number_fails)
    for i in range(1,len(number_fails)):
        print(f"{i} changes was detected in {number_fails[i]} times ")

    # domain x chains which was missed in 2 snapshot and added then
    counter_all = {}
    for fam in all_families:
        counter = {}
        for i in range(2, len(dates)):
            snapshot1, snapshot2, snapshot3 = dates[i-2], dates[i-1], dates[i]
            dom_chain_dif_list2 = stat_data[snapshot1][fam].difference(stat_data[snapshot2][fam])
            dom_chain_dif_list3 = stat_data[snapshot1][fam].difference(stat_data[snapshot3][fam])
            dom_chain_dif_list = dom_chain_dif_list2.intersection(dom_chain_dif_list3)
            for dom_chain in dom_chain_dif_list:
               if (dom_chain,fam) not in counter:
                  counter[(dom_chain,fam)] = 1
               else:
                  counter[(dom_chain,fam)] += 1
               counter_all[(dom_chain, fam)]=counter[(dom_chain, fam)]
    val = set(sorted(list(counter_all.values())))
    number_fails={}
    for i in counter_all.keys():
       if counter_all[i] not in number_fails:
           number_fails[counter_all[i]] = 1
       else:
           number_fails[counter_all[i]] += 1
    print("domain chains which was missed in 2 snapshot and then added", number_fails)

    # domain x chains which was missed in 2 snapshot and added then
    counter_all = {}
    for fam in all_families:
        counter = {}
        for i in range(3, len(dates)):
            snapshot1, snapshot2, snapshot3, snapshot4 = dates[i-3], dates[i-2], dates[i-1], dates[i]
            dom_chain_dif_list2 = stat_data[snapshot1][fam].difference(stat_data[snapshot2][fam])
            dom_chain_dif_list3 = stat_data[snapshot1][fam].difference(stat_data[snapshot3][fam])
            dom_chain_dif_list4 = stat_data[snapshot1][fam].difference(stat_data[snapshot4][fam])
            dom_chain_dif_list23 = dom_chain_dif_list2.intersection(dom_chain_dif_list3)
            dom_chain_dif_list = dom_chain_dif_list23.intersection(dom_chain_dif_list4)
            for dom_chain in dom_chain_dif_list:
               if (dom_chain,fam) not in counter:
                  counter[(dom_chain,fam)] = 1
               else:
                  counter[(dom_chain,fam)] += 1
               counter_all[(dom_chain, fam)]=counter[(dom_chain, fam)]
    val = set(sorted(list(counter_all.values())))
    number_fails={}
    for i in counter_all.keys():
       if counter_all[i] not in number_fails:
           number_fails[counter_all[i]] = 1
       else:
           number_fails[counter_all[i]] += 1
    print("domain chains which was missed in 3 snapshot and then added", number_fails)

    # domain x chains which was missed in 2 snapshot and added then
    counter_all = {}
    for fam in all_families:
        counter = {}
        for i in range(4, len(dates)):
            snapshot1, snapshot2, snapshot3, snapshot4, snapshot5 = dates[i-4], dates[i-3], dates[i-2], dates[i-1], dates[i]
            dom_chain_dif_list2 = stat_data[snapshot1][fam].difference(stat_data[snapshot2][fam])
            dom_chain_dif_list3 = stat_data[snapshot1][fam].difference(stat_data[snapshot3][fam])
            dom_chain_dif_list4 = stat_data[snapshot1][fam].difference(stat_data[snapshot4][fam])
            dom_chain_dif_list5 = stat_data[snapshot1][fam].difference(stat_data[snapshot5][fam])
            dom_chain_dif_list23 = dom_chain_dif_list2.intersection(dom_chain_dif_list3)
            dom_chain_dif_list234 = dom_chain_dif_list23.intersection(dom_chain_dif_list4)
            dom_chain_dif_list = dom_chain_dif_list234.intersection(dom_chain_dif_list5)
            for dom_chain in dom_chain_dif_list:
               if (dom_chain,fam) not in counter:
                  counter[(dom_chain,fam)] = 1
               else:
                  counter[(dom_chain,fam)] += 1
               counter_all[(dom_chain, fam)]=counter[(dom_chain, fam)]
    val = set(sorted(list(counter_all.values())))
    number_fails={}
    for i in counter_all.keys():
       if counter_all[i] not in number_fails:
           number_fails[counter_all[i]] = 1
       else:
           number_fails[counter_all[i]] += 1
    print("domain chains which was missed in 4 snapshot and then added", number_fails)

    # domain x chains which was missed in 2 snapshot and added then
    counter_all = {}
    for fam in all_families:
        counter = {}
        for i in range(5, len(dates)):
            snapshot1, snapshot2, snapshot3, snapshot4, snapshot5, snapshot6 = dates[i-5], dates[i-4], dates[i-3], dates[i-2], dates[i-1], dates[i]
            dom_chain_dif_list2 = stat_data[snapshot1][fam].difference(stat_data[snapshot2][fam])
            dom_chain_dif_list3 = stat_data[snapshot1][fam].difference(stat_data[snapshot3][fam])
            dom_chain_dif_list4 = stat_data[snapshot1][fam].difference(stat_data[snapshot4][fam])
            dom_chain_dif_list5 = stat_data[snapshot1][fam].difference(stat_data[snapshot5][fam])
            dom_chain_dif_list6 = stat_data[snapshot1][fam].difference(stat_data[snapshot6][fam])
            dom_chain_dif_list23 = dom_chain_dif_list2.intersection(dom_chain_dif_list3)
            dom_chain_dif_list234 = dom_chain_dif_list23.intersection(dom_chain_dif_list4)
            dom_chain_dif_list2345 = dom_chain_dif_list234.intersection(dom_chain_dif_list5)
            dom_chain_dif_list = dom_chain_dif_list2345.intersection(dom_chain_dif_list6)
            for dom_chain in dom_chain_dif_list:
               if (dom_chain,fam) not in counter:
                  counter[(dom_chain,fam)] = 1
               else:
                  counter[(dom_chain,fam)] += 1
               counter_all[(dom_chain, fam)]=counter[(dom_chain, fam)]
    val = set(sorted(list(counter_all.values())))
    number_fails={}
    for i in counter_all.keys():
       if counter_all[i] not in number_fails:
           number_fails[counter_all[i]] = 1
       else:
           number_fails[counter_all[i]] += 1
    print("domain chains which was missed in 5 snapshot and then added", number_fails)
    # domain x chains which belonges to more families
    counter_all = {}
    for snapshot in dates:
        for fam in all_families:
            for dom_chain in stat_data[snapshot][fam]:
                if dom_chain not in counter_all:
                    counter_all[dom_chain] = {fam}
                else:
                    counter_all[dom_chain].add(fam)
    number_fails={}
    for i in counter_all.keys():
       if len(counter_all[i]) not in number_fails:
           number_fails[len(counter_all[i])] = 1
       else:
           number_fails[len(counter_all[i])] += 1
    print("number of families to which domain chain belonges", number_fails)
    print("====================================================================")

    # domain x chains which was missed till x. snapshot and added then
    ldates=len(dates)
    for i in range(2, ldates):
        counter_before = {}
        counter_after = {}
        for j in range(1, i):
           snapshot = dates[j]
           for fam in all_families:
               for dom_chain in stat_data[snapshot][fam]:
                   if (dom_chain, fam) not in counter_before:
                       counter_before[(dom_chain, fam)] = 1
                   else:
                       counter_before[(dom_chain, fam)] += 1

        for j in range(i, ldates):
           snapshot = dates[j]
           for fam in all_families:
               for dom_chain in stat_data[snapshot][fam]:
                   if (dom_chain, fam) not in counter_after:
                       counter_after[(dom_chain, fam)] = 1
                   else:
                       counter_after[(dom_chain, fam)] += 1

        added=set()
        removed=set()
        r=len(range(i+1, ldates))
        for j in counter_after.keys():
            if (counter_after[j] == r) and (j not in counter_before):
                added.add(j)
        for j in counter_before.keys():
            if (counter_before[j] == i-1) and (j not in counter_after):
                removed.add(j)
        res = len(number_fails)
        if (i >3 and i < 19) or (i ==19):
            print(f"{dates[i]} {i} added {len(added)} removed {len(removed)} { added}")
        if (i == 16):
           families = set()
           for k in removed:
               families.add(k[1])
           for k in families:
               print(k)





stat_data={}
for i in sys.argv[1:]:
    read_stats(i, stat_data)

week_stats(stat_data)
family_stats(stat_data)