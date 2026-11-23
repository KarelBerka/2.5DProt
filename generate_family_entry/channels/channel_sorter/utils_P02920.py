#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import json
import utils_read_cif

MODEL = ["1pv6", "A", "00"]


def sses_model_mapping_1pv6(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 7 and sse["end"] == 38:
                mapping["1"] = sse["label"]
            if sse["start"] == 42 and sse["end"] == 70:
                mapping["2"] = sse["label"]
            if sse["start"] == 74 and sse["end"] == 86:
                mapping["3a"] = sse["label"]
            if sse["start"] == 87 and sse["end"] == 101:
                mapping["3b"] = sse["label"]
            if sse["start"] == 104 and sse["end"] == 136:
                mapping["4"] = sse["label"]
            if sse["start"] == 140 and sse["end"] == 164:
                mapping["5"] = sse["label"]
            if sse["start"] == 166 and sse["end"] == 187:
                mapping["6"] = sse["label"]
            if sse["start"] == 210 and sse["end"] == 219:
                mapping["7a"] = sse["label"]
            if sse["start"] == 220 and sse["end"] == 250:
                mapping["7b"] = sse["label"]
            if sse["start"] == 254 and sse["end"] == 286:
                mapping["8"] = sse["label"]
            if sse["start"] == 288 and sse["end"] == 308:
                mapping["9"] = sse["label"]
            if sse["start"] == 312 and sse["end"] == 340:
                mapping["10"] = sse["label"]
            if sse["start"] == 346 and sse["end"] == 376:
                mapping["11"] = sse["label"]
            if sse["start"] == 378 and sse["end"] == 399:
                mapping["12a"] = sse["label"]
            if sse["start"] == 408 and sse["end"] == 416:
                mapping["12b"] = sse["label"]
    return mapping


def predict_P02920_tunnel_start_end(annotation_dir, pcd):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_1pv6(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    s1, s2, s3, e1, e2 = None, None, None, None, None
    print("mapping", mapping["2"], mapping["5"], mapping["7b"])
    print("mapping", mapping["4"], mapping["11"])
    for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
        if sse['label'] == mapping["2"]:   # 79
            s1 = sse["start"]
        if sse['label'] == mapping["5"]:   # 195
            s2 = sse["end"]
        if sse['label'] == mapping["7b"]:   # 36
            s3 = sse["end"]
#        if sse['label'] == mapping["12"]:   # 36 # necessaty for Q8EHE6 - don't have 7
#            s3 = sse["start"]
        if sse['label'] == mapping["4"]:   # 165
            e1 = sse["end"]
        if sse['label'] == mapping["11"]:   # 295
            e2 = sse["start"]
        print("sse", sse['label'], )

    if e1 != None and e2 != None and s1 != None and s2 != None and s3!= None:
        # transform start and end residue to coordinate in the original (not rotated) cif file
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        print("======================")
        print(s1, s2, s3, e1, e2)
        start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
        start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
        start3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s3)
        start = [(start1[0]+start2[0]+start3[0])/3,
                 (start1[1]+start2[1]+start3[1])/3,
                 (start1[2]+start2[2]+start3[2])/3]

        end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
        end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
        print("end1", end1)
        print("end2", end2)
        print("start1", start1)
        print("start2", start2)
        print("start3", start3)
        
        end = [(end1[0]+end2[0])/2,
               (end1[1]+end2[1])/2,
               (end1[2]+end2[2])/2]

        return start, end
#        return None, end
    else:
        return None, None




