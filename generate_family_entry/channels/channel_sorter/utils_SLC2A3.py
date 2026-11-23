#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import json
import utils_read_cif

MODEL = ["4zwb", "A", "00"]


def sses_model_mapping_4zwb(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 30 and sse["end"] == 56:
                mapping["1"] = sse["label"]
            if sse["start"] == 57 and sse["end"] == 71:
                mapping["1e"] = sse["label"]
            if sse["start"] == 78 and sse["end"] == 110:
                mapping["2"] = sse["label"]
            if sse["start"] == 112 and sse["end"] == 136:
                mapping["3"] = sse["label"]
            if sse["start"] == 139 and sse["end"] == 167:
                mapping["4"] = sse["label"]
            if sse["start"] == 173 and sse["end"] == 195:
                mapping["5"] = sse["label"]
            if sse["start"] == 206 and sse["end"] == 226:
                mapping["6"] = sse["label"]
            if sse["start"] == 286 and sse["end"] == 305:
                mapping["7a"] = sse["label"]
            if sse["start"] == 307 and sse["end"] == 321:
                mapping["7b"] = sse["label"]
            if sse["start"] == 326 and sse["end"] == 351:
                mapping["8"] = sse["label"]
            if sse["start"] == 353 and sse["end"] == 377:
                mapping["9"] = sse["label"]
            if sse["start"] == 384 and sse["end"] == 401:
                mapping["10"] = sse["label"]
            if sse["start"] == 417 and sse["end"] == 449:
                mapping["11"] = sse["label"]
            if sse["start"] == 450 and sse["end"] == 471:
                mapping["12"] = sse["label"]
    return mapping


def find_helix_start_end(annotation_dir, pcd, helix):
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_4zwb(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
        if sse['label'] == mapping[helix]:   # 79
            return sse['start'], sse['end']
    return None, None



def predict_SLC2A3_tunnel_start_end(annotation_dir, pcd):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_4zwb(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    s1, s2, s3, e1, e2 = None, None, None, None, None
    for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
        if sse['label'] == mapping["2"]:   # 79
            s1 = sse["start"]
        if sse['label'] == mapping["11"]:   # 195
            s2 = sse["end"]
        if sse['label'] == mapping["7b"]:   #
            s3 = sse["end"]
        if sse['label'] == mapping["12"]:   # 36 # necessaty for Q8EHE6 - don't have 7
            s3 = sse["start"]
        if sse['label'] == mapping["4"]:   # 165
            e1 = sse["end"]
        if sse['label'] == mapping["11"]:   # 295
            e2 = sse["start"]
        if sse['label'] == mapping["7a"]:   # hack for uncomplete proteins (e.g. P14672-2)
            e2 = sse["start"]

    if e1 != None and e2 != None and s1 != None and s2 != None and s3!= None:
        # transform start and end residue to coordinate in the original (not rotated) cif file
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
        start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
        start3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s3)
        start = [(start1[0]+start2[0]+start3[0])/3,
                 (start1[1]+start2[1]+start3[1])/3,
                 (start1[2]+start2[2]+start3[2])/3]
        end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
        end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
        end = [(end1[0]+end2[0])/2,
               (end1[1]+end2[1])/2,
               (end1[2]+end2[2])/2]

        return start, end
    else:
        return None, None




