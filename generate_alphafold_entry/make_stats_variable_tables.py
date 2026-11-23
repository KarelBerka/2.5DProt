#!/usr/bin/env python3
# tool which generate for given protein its directory with
# * protein images
# and image per chain

import argparse
import os
import json
import collections


def add_statistics_from_file(filename, name, stats):
    with open(filename, "r") as fd:
        data = json.load(fd)
        for var in data:
            if type(data[var]) == float or type(data[var]) == int:
                if var not in stats:
                    stats[var] = {}
                stats[var][name] = data[var]
            else:
                if var not in stats:
                    stats[var] = {}
                for val in data[var]:
                    if int(val) not in stats[var]:
                        stats[var][int(val)] = {}
                    stats[var][int(val)][name] = data[var][val]
    return True


def create_stats_from_dir(af_dir):
    # from the input directory, find out the list of proteins
    organism_list = []
    stats = {}
    stats["sses"] = {}
    stats["helices"] = {}
    stats["strands"] = {}
    stats["helices_length"] = {}
    stats["strands_length"] = {}
    stats["len"] = {}
    stats["helices_len"] = {}
    stats["strands_len"] = {}
    stats["proteins"] = {}
    stats["avg_helices_number"] = {}
    stats["avg_strands_number"] = {}
    stats["avg_ses_number"] = {}
    stats["avg_helices_length"] = {}
    stats["avg_strands_length"] = {}
    stats["avg_helices_len"] = {}
    stats["avg_strands_len"] = {}
    stats["med_helices_number"] = {}
    stats["med_strands_number"] = {}
    stats["med_ses_number"] = {}
    stats["med_helices_length"] = {}
    stats["med_strands_length"] = {}
    stats["med_helices_len"] = {}
    stats["med_strands_len"] = {}

    stats = {}
    for filename in os.listdir(af_dir):
        if filename.endswith(".json"):
            add_statistics_from_file(f"{af_dir}/{filename}", filename[:-5], stats)
            organism_list.append(filename[:-5])
    return stats, organism_list


def output(organism_list, stats, stra=True, heli=True, both=True):

  if stra:
    print("strands")
    for o in sorted(organism_list):
        print(o, end=" ")
    print()
    for o in sorted(organism_list):
        for oo in sorted(organism_list):
            diff = 0
            o_strands_sum = stats["proteins"][o]*stats["avg_strands_number"][o]
            oo_strands_sum = stats["proteins"][oo]*stats["avg_strands_number"][o]
            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                    if oo in stats["strands"][i]:
                        strands_diff = stats["strands"][i][o]/o_strands_sum*100 - stats["strands"][i][oo]/oo_strands_sum*100
                    else:
                        strands_diff = stats["strands"][i][o]/o_strands_sum*100
                elif oo in stats["strands"][i]:
                    strands_diff = stats["strands"][i][oo]/oo_strands_sum*100
                else:
                    strands_diff = 0
                diff += strands_diff*strands_diff*10

  if heli:
    print("helices")
    for o in sorted(organism_list):
        print(o, end=" ")
    print()
    for o in sorted(organism_list):
        print(o, end=" ")
        for oo in sorted(organism_list):
            diff = 0
            o_helices_sum = stats["proteins"][o]*stats["avg_helices_number"][o]
            oo_helices_sum = stats["proteins"][oo]*stats["avg_helices_number"][o]
            for i in sorted(stats["helices"]):
                if o in stats["helices"][i]:
                    if oo in stats["helices"][i]:
                        helices_diff = stats["helices"][i][o]/o_helices_sum*100 - stats["helices"][i][oo]/oo_helices_sum*100
                    else:
                        helices_diff = stats["helices"][i][o]/o_helices_sum*100
                elif oo in stats["helices"][i]:
                    helices_diff = stats["helices"][i][oo]/oo_helices_sum*100
                else:
                    helices_diff = 0
                diff += helices_diff*helices_diff*10

  if both:
    print("helices + strands")
    for o in sorted(organism_list):
        print(o, end=" ")
    print()

    for o in sorted(organism_list):
        print(o, end=" ")
        for oo in sorted(organism_list):
            diff = 0
            o_strands_sum = stats["proteins"][o]*stats["avg_strands_number"][o]
            o_helices_sum = stats["proteins"][o]*stats["avg_helices_number"][o]
            oo_helices_sum = stats["proteins"][oo]*stats["avg_helices_number"][oo]
            oo_strands_sum = stats["proteins"][oo]*stats["avg_strands_number"][o]
            for i in sorted(stats["helices"]):
                if o in stats["helices"][i]:
                    if oo in stats["helices"][i]:
                        helices_diff = stats["helices"][i][o]/o_helices_sum*100 - stats["helices"][i][oo]/oo_helices_sum*100
                    else:
                        helices_diff = stats["helices"][i][o]/o_helices_sum*100
                elif oo in stats["helices"][i]:
                    helices_diff = stats["helices"][i][oo]/oo_helices_sum*100
                else:
                    helices_diff = 0
                diff += helices_diff*helices_diff*10

            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                    if oo in stats["strands"][i]:
                        strands_diff = stats["strands"][i][o]/o_strands_sum*100 - stats["strands"][i][oo]/oo_strands_sum*100
                    else:
                        strands_diff = stats["strands"][i][o]/o_strands_sum*100
                elif oo in stats["strands"][i]:
                    strands_diff = stats["strands"][i][oo]/oo_strands_sum*100
                else:
                    strands_diff = 0
                diff += strands_diff*strands_diff*10

        print()


def output_area(organism_list, stats, stra=True, heli=True, both=True):
    print("helices + strands")
    print(" ", end="")
    ol = ['Methanocaldococcus_without', 'Nocardia_without', 'Myle', 'MytuMycobacterium_tuberculosis_without', 'Streptococcu_without', 'Staphylococcus_without', 'Helicobacter_without', 'Campylobacter_without', 'Neisseria_without', 'Pseudomonas_without', 'Haemophilus_without', 'Klebsiella_without', 'Salmonella_without', 'Shigella_without', 'Escherichia_without', 'Leishmania-2022-12-08T19_06_12_without', 'Trypanosoma_cruzi-2022-12-10T04_39_58_without', 'Trypanosoma_brucei-2022-12-09T15_18_47_without', 'Schizosaccharomyces_without', 'Saccharomyces_without', 'Candida_without', 'Sporothrix_without', 'Madurella', 'Fonsecaea_without', 'Cladophialophora_without', 'Paracoccidioides_without', 'Ajellomyces_without', 'Plasmodium-2022-12-09T08_15_25_without', 'Dictypstelium_without', 'Arabidopsis_without', 'Glycine_without', 'Zea_without', 'Oryza_without', 'Caenorhabditis_without', 'Dracunculus-2022-12-10T19_00_55_without', 'Onchocerca_without', 'Brugia_without', 'Drosophila_without', 'Danio_without', 'Mus_without', 'Rattus_without', 'Octodon_degus_without', 'Rhinolophus_ferrumequinum_without', 'Rousettus_aegyptiacus_without', 'Phyllostomus_discolor_without', 'Molossus_molossus_without', 'Myotis_myotis_without', 'Sapajus_apella_without', 'Macaca_mulatta_without', 'Cercocebus_atys_without', 'Macaca_nemestrina_without', 'Gorilla_gorilla_gorilla_without', 'Homo_sapiens_without', 'Pan_paniscus_without', 'Pan_troglodytes_Chimpanzee_without', 'Felis_catus_without', 'Acinonyx_jubatus_without', 'Vulpes_culpes_without', 'Nyctereutes_procyonoides_without', 'Ailuropoda_melanoleuca_without', 'Enhydra_lutris_kenyoni_without', 'Callorhinus_ursinus_without', 'Odobenus_rosmarus_divergens_without', 'Odocoileus_virginianus_texanus_without', 'Capra_hircus_without', 'Bos_indicus_without', 'Bison_bison_bison_without', 'Vicugna_pacos_without', 'Balaenoptera_acutorostrata_scammoni_without', 'Physeter_macrocephalus_without', 'Lipotes_vexillifer_without', 'Tursiops_truncatus_without', 'Delphinapterus_leucas_without', 'Enterococcus_without']
    organism_list = ol
    for o in organism_list:
        print(o, end=" ")
    print()

    for o in organism_list:
        print(o, end=" ")
        for oo in organism_list:
            diff = 0
            o2_strands_sum = 0
            o2_helices_sum = 0
            oo2_strands_sum = 0
            oo2_helices_sum = 0

            for i in sorted(stats["helices"]):
                if o in stats["helices"][i]:
                     o2_helices_sum += stats["helices"][i][o]
                if oo in stats["helices"][i]:
                     oo2_helices_sum += stats["helices"][i][oo]

            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                     o2_strands_sum += stats["strands"][i][o]
                if oo in stats["strands"][i]:
                     oo2_strands_sum += stats["strands"][i][oo]

            for i in sorted(stats["helices"]):
                if i > 20:
                    break 
                if o in stats["helices"][i]:
                    if oo in stats["helices"][i]:
                        helices_diff = abs(stats["helices"][i][o]/o2_helices_sum*100 - stats["helices"][i][oo]/oo2_helices_sum*100)
                    else:
                        helices_diff = abs(stats["helices"][i][o]/o2_helices_sum*100)
#                        print(f"{stats['helices'][i][o]:.1f} ", end="")
                elif oo in stats["helices"][i]:
                    helices_diff = abs(stats["helices"][i][oo]/oo2_helices_sum*100)
                else:
                    helices_diff = 0
#                    print(f"0 ", end="")
##                diff += helices_diff/2#*helices_diff#*10
#                print(f"{helices_diff:.1f}({stats['helices'][i][o]/o2_helices_sum*100:.1f}-{stats['helices'][i][oo]/oo2_helices_sum*100:.1f}) ", end="")

            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                    if oo in stats["strands"][i]:
                        strands_diff = abs(stats["strands"][i][o]/o2_strands_sum*100 - stats["strands"][i][oo]/oo2_strands_sum*100)
                    else:
                        strands_diff = abs(stats["strands"][i][o]/o2_strands_sum*100)
                elif oo in stats["strands"][i]:
                    strands_diff = abs(stats["strands"][i][oo]/oo2_strands_sum*100)
                else:
                    strands_diff = 0
                diff += strands_diff/2*2

            print(f"{int(diff)}", end=" ")
#            print(f"| {diff:.1f} - {o_helices_sum} o2 {o2_helices_sum }", end=" ")
        print()


def output_max(organism_list, stats, stra=True, heli=True, both=True):
    print(" ", end="")
    ol = ['Methanocaldococcus_without', 'Nocardia_without', 'Myle', 'MytuMycobacterium_tuberculosis_without', 'Streptococcu_without', 'Staphylococcus_without', 'Helicobacter_without', 'Campylobacter_without', 'Neisseria_without', 'Pseudomonas_without', 'Haemophilus_without', 'Klebsiella_without', 'Salmonella_without', 'Shigella_without', 'Escherichia_without', 'Leishmania-2022-12-08T19_06_12_without', 'Trypanosoma_cruzi-2022-12-10T04_39_58_without', 'Trypanosoma_brucei-2022-12-09T15_18_47_without', 'Schizosaccharomyces_without', 'Saccharomyces_without', 'Candida_without', 'Sporothrix_without', 'Madurella', 'Fonsecaea_without', 'Cladophialophora_without', 'Paracoccidioides_without', 'Ajellomyces_without', 'Plasmodium-2022-12-09T08_15_25_without', 'Dictypstelium_without', 'Arabidopsis_without', 'Glycine_without', 'Zea_without', 'Oryza_without', 'Caenorhabditis_without', 'Dracunculus-2022-12-10T19_00_55_without', 'Onchocerca_without', 'Brugia_without', 'Drosophila_without', 'Danio_without', 'Mus_without', 'Rattus_without', 'Octodon_degus_without', 'Rhinolophus_ferrumequinum_without', 'Rousettus_aegyptiacus_without', 'Phyllostomus_discolor_without', 'Molossus_molossus_without', 'Myotis_myotis_without', 'Sapajus_apella_without', 'Macaca_mulatta_without', 'Cercocebus_atys_without', 'Macaca_nemestrina_without', 'Gorilla_gorilla_gorilla_without', 'Homo_sapiens_without', 'Pan_paniscus_without', 'Pan_troglodytes_Chimpanzee_without', 'Felis_catus_without', 'Acinonyx_jubatus_without', 'Vulpes_culpes_without', 'Nyctereutes_procyonoides_without', 'Ailuropoda_melanoleuca_without', 'Enhydra_lutris_kenyoni_without', 'Callorhinus_ursinus_without', 'Odobenus_rosmarus_divergens_without', 'Odocoileus_virginianus_texanus_without', 'Capra_hircus_without', 'Bos_indicus_without', 'Bison_bison_bison_without', 'Vicugna_pacos_without', 'Balaenoptera_acutorostrata_scammoni_without', 'Physeter_macrocephalus_without', 'Lipotes_vexillifer_without', 'Tursiops_truncatus_without', 'Delphinapterus_leucas_without', 'Enterococcus_without']
    organism_list = ol
    for o in organism_list:
        print(o, end=" ")
    print()

    for o in organism_list:
        print(o, end=" ")
        for oo in organism_list:
#            print(o,oo, end = " ")
            diff = 0 
            o2_strands_max = 0
            o2_helices_max = 0
            oo2_strands_max = 0
            oo2_helices_max = 0
            o2_strands_sum = 0
            o2_helices_sum = 0
            oo2_strands_sum = 0
            oo2_helices_sum = 0

            for i in sorted(stats["helices"]):
                if o in stats["helices"][i]:
                    o2_helices_sum += stats["helices"][i][o]
                if oo in stats["helices"][i]:
                    oo2_helices_sum += stats["helices"][i][oo]

            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                    o2_strands_sum += stats["strands"][i][o]
                if oo in stats["strands"][i]:
                    oo2_strands_sum += stats["strands"][i][oo]

            for i in sorted(stats["helices"]):
                if o in stats["helices"][i]:
                    if stats["helices"][i][o]/o2_helices_sum*100 > o2_helices_max:
                        o2_helices_max = stats["helices"][i][o]/o2_helices_sum*100
                if oo in stats["helices"][i]:
                    if stats["helices"][i][oo] > oo2_helices_max/oo2_helices_sum*100:
                        oo2_helices_max = stats["helices"][i][oo]/oo2_helices_sum*100
            if o2_helices_max > oo2_helices_max:
                helices_max = o2_helices_max
            else:
                helices_max = oo2_helices_max

            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                    if stats["strands"][i][o] > o2_strands_max:
                        o2_strands_max = stats["strands"][i][o]
                if oo in stats["strands"][i]:
                    if stats["strands"][i][oo] > oo2_strands_max:
                        oo2_strands_max = stats["strands"][i][oo]
            if o2_strands_max > oo2_strands_max:
                strands_max = o2_strands_max
            else:
                strands_max = oo2_strands_max

            for i in sorted(stats["helices"]):
                if i > 20:
                    break 
                if o in stats["helices"][i]:
                    if oo in stats["helices"][i]:
                        helices_diff = abs(stats["helices"][i][o]/o2_helices_sum*100 - stats["helices"][i][oo]/oo2_helices_sum*100)
                    else:
                        helices_diff = abs(stats["helices"][i][o]/o2_helices_sum*100)
#                        print(f"{stats['helices'][i][o]:.1f} ", end="")
                elif oo in stats["helices"][i]:
                    helices_diff = abs(stats["helices"][i][oo]/oo2_helices_sum*100)
                else:
                    helices_diff = 0
#                if diff < helices_diff:
#                    diff = helices_diff

            for i in sorted(stats["strands"]):
                if o in stats["strands"][i]:
                    if oo in stats["strands"][i]:
                        strands_diff = abs(stats["strands"][i][o]/o2_strands_max*100 - stats["strands"][i][oo]/oo2_strands_max*100)
                    else:
                        strands_diff = abs(stats["strands"][i][o]/o2_strands_max*100)
                elif oo in stats["strands"][i]:
                    strands_diff = abs(stats["strands"][i][oo]/oo2_strands_max*100)
                else:
                    strands_diff = 0
                if diff < strands_diff:
                    diff = strands_diff

            print(f"{int(diff)}", end=" ")
#            print()

#            print(f"| {diff:.1f} - {o_helices_sum} o2 {o2_helices_sum }", end=" ")
        print()


def output_hist_diff(organism_list, stats):
    for var in stats:
        print(var)
        print("", end=" ")
        for o in organism_list:
            print(o, end=" ")
        print()
        if var == "proteins" or var == "avg_helices_number" or var == "avg_strands_number" or var == "avg_ses_number" or \
           var == "avg_helices_length" or var == "avg_strands_length" or var == "avg_helices_len" or var == "avg_strands_len" or \
           var == "med_helices_number" or var == "med_strands_number" or var == "med_ses_number" or var == "med_helices_length" or \
           var == "med_strands_length" or var == "med_helices_len" or var == "med_strands_len":
            print("", end=" ")
            for o in organism_list:
                print(stats[var][o], end=" ")
            print()
        else:
            val_list = collections.OrderedDict(sorted(stats[var].items()))
            for val in val_list:
                print(val, end=" ")
                for o in organism_list:
                    if var == "helices_len" or var == "helices_length":
                        if o in stats[var][val]:
                            print(stats[var][val][o]*100/stats["proteins"][o]/stats["avg_helices_number"][o], end=" ")
                        else:
                            print("0", end=" ")
                    elif var == "strands_len" or var == "strands_length":
                        if o in stats[var][val]:
                            print(stats[var][val][o]*100/stats["proteins"][o]/stats["avg_strands_number"][o], end=" ")
                        else:
                            print("0", end=" ")
                    else:
                        if o in stats[var][val]:
                            print(stats[var][val][o]*100/stats["proteins"][o], end=" ")
                        else:
                            print("0", end=" ")
                print()


def add_to_dict(input_dict, key, key_value):
    if key not in input_dict:
        input_dict[key] = key_value
    else:
        input_dict[key] += key_value


# read arguments
parser = argparse.ArgumentParser(
    description='generate diagram for a alphafold structure')
parser.add_argument('--organisms_dir', action="store", required=True,
                    help="generate database stats for alfa fold structures in the given directory")
args = parser.parse_args()

# now we have to go through all -detected- files and count stats
stats, organism_list = create_stats_from_dir(args.organisms_dir)

# output(organism_list, stats)
# output_area(organism_list, stats)
output_max(organism_list, stats)
# output_hist_diff(organism_list, stats)
