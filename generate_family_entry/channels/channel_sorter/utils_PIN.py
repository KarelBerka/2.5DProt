#!/usr/bin/env python3
# go through all channels in given directory and find out group which goes through the given polygon
# defined by SSEs positions
# e.g. for cytochrpmes channel3 tunnels go through polygon defined by SSEs F and G

import json
import utils_read_cif

MODEL = ["V7B1K1", "A", "00"]


def sses_model_mapping_v7b1k1a(model_annotation_file):
    # maping based on  https://sestra.ncbr.muni.cz/job/cyp-1akd
    mapping = {}
    with open(model_annotation_file, "r") as fr:
        data = json.load(fr)
        pcd = f"{MODEL[0]}{MODEL[1]}{MODEL[2]}"
        for sse in data[pcd]["secondary_structure_elements"]:
            if sse["start"] == 4 and sse["end"] == 31:
                mapping["A1"] = sse["label"]
            if sse["start"] == 36 and sse["end"] == 61:
                mapping["A2"] = sse["label"]
            if sse["start"] == 69 and sse["end"] == 92:
                mapping["A3"] = sse["label"]
            if sse["start"] == 98 and sse["end"] == 110:
                mapping["A4a"] = sse["label"]
            if sse["start"] == 113 and sse["end"] == 125:
                mapping["A4b"] = sse["label"]
            if sse["start"] == 127 and sse["end"] == 165:
                mapping["A5"] = sse["label"]
            if sse["start"] == 252 and sse["end"] == 269:
                mapping["A6a"] = sse["label"]
            if sse["start"] == 271 and sse["end"] == 289:
                mapping["A6b"] = sse["label"]
            if sse["start"] == 295 and sse["end"] == 322:
                mapping["A7"] = sse["label"]
            if sse["start"] == 331 and sse["end"] == 356:
                mapping["A8"] = sse["label"]
            if sse["start"] == 360 and sse["end"] == 371:
                mapping["A9a"] = sse["label"]
            if sse["start"] == 376 and sse["end"] == 385:
                mapping["A9b"] = sse["label"]
            if sse["start"] == 389 and sse["end"] == 415:
                mapping["A10"] = sse["label"]
    return mapping


def predict_PIN_tunnel_start_end(annotation_dir, pcd):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_v7b1k1a(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    s1, s2, e1, e2 = None, None, None, None
    for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
        if sse['label'] == mapping["A5"]:   # 127
            e1 = sse["start"]
        if sse['label'] == mapping["A4b"]:   # 125
            e2 = sse["end"]
        if sse['label'] == mapping["A2"]:   # 36
            s1 = sse["start"]
        if sse['label'] == mapping["A5"]:   # 165
            s2 = sse["end"]
        if sse['label'] == mapping["A7"]:   # 295
            e4 = sse["start"]
        if sse['label'] == mapping["A2"]:   # 61
            e5 = sse["end"]

    if e1 != None and e2 != None and s1 != None and s2 != None:
        # transform start and end residue to coordinate in the original (not rotated) cif file
        cif_file = f"{annotation_dir}/{pcd[0]}.cif"
        start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
        start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
        start = [(start1[0]+start2[0])/2, (start1[1]+start2[1])/2, (start1[2]+start2[2])/2]

        end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
        end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)

        end4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e4)
        end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
        end = [(end1[0]+end2[0]+end5[0]+end4[0])/4,
               (end1[1]+end2[1]+end5[1]+end4[1])/4,
               (end1[2]+end2[2]+end5[2]+end4[2])/4]

        return start, end
    else:
        return None, None


def predict_PIN_tunnel_start(annotation_dir, pcd, number):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_v7b1k1a(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    s1, s2, s3, s4 = None, None, None, None

    if number == 1:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A2"]:   # 36
                s1 = sse["start"]
            if sse['label'] == mapping["A5"]:   # 165
                s2 = sse["end"]
        if s1 != None and s2 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
            start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
            start = [(start1[0]+start2[0])/2,
                     (start1[1]+start2[1])/2,
                     (start1[2]+start2[2])/2]
            return start

    elif number == 2:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A6a"]:    #
                s1 = sse["start"]
            elif sse['label'] == mapping["A5"]:   #
                s2 = sse["end"]
            elif sse['label'] == mapping["A4a"]:  #
                s3 = sse["start"]
            elif sse['label'] == mapping["A3"]:   #
                s4 = sse["end"]
        if s1 != None and s2 != None and s3 != None and s4 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
            start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
            start3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s3)
            start4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s4)
            start = [(start1[0]+start2[0]+start3[0]+start4[0])/4,
                     (start1[1]+start2[1]+start3[1]+start4[1])/4,
                     (start1[2]+start2[2]+start3[2]+start4[2])/4]
            return start

    elif number == 3:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A10"]:     #
                s1 = sse["start"]
            elif sse['label'] == mapping["A5"]:    #
                s2 = sse["end"]
            elif sse['label'] == mapping["A4a"]:   #
                s3 = sse["start"]
            elif sse['label'] == mapping["A9b"]:   #
                s4 = sse["end"]
        if s1 != None and s2 != None and s3 != None and s4 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
            start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
            start3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s3)
            start4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s4)
            start = [(start1[0]+start2[0]+start3[0]+start4[0])/4,
                     (start1[1]+start2[1]+start3[1]+start4[1])/4,
                     (start1[2]+start2[2]+start3[2]+start4[2])/4]
            return start

    elif number == 4:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A10"]:     #
                s1 = sse["start"]
            elif sse['label'] == mapping["A4a"]:   #
                s2 = sse["start"]
            elif sse['label'] == mapping["A8"]:    #
                s3 = sse["start"]
        if s1 != None and s2 != None and s3 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
            start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
            start3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s3)
            start = [(start1[0]+start2[0]+start3[0])/3,
                     (start1[1]+start2[1]+start3[1])/3,
                     (start1[2]+start2[2]+start3[2])/3]
            return start

    elif number == 5:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A6a"]:     #
                s1 = sse["start"]
            elif sse['label'] == mapping["A5"]:    #
                s2 = sse["end"]
            elif sse['label'] == mapping["A2"]:   #
                s3 = sse["start"]
        if s1 != None and s2 != None and s3 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            start1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s1)
            start2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s2)
            start3 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], s3)
            start = [(start1[0]+start2[0]+start3[0])/3,
                     (start1[1]+start2[1]+start3[1])/3,
                     (start1[2]+start2[2]+start3[2])/3]
            return start

    return None


def predict_PIN_tunnel_end(annotation_dir, pcd, number):
    # find start and end of membrane protein
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_v7b1k1a(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    e1, e2, e4, e5 = None, None, None, None
    if number == 1:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A5"]:   # 127
                e1 = sse["start"]
            if sse['label'] == mapping["A4b"]:   # 125
                e2 = sse["end"]
            if sse['label'] == mapping["A7"]:   # 295
                e4 = sse["start"]
            if sse['label'] == mapping["A2"]:   # 61
                e5 = sse["end"]

        if e1 != None and e2 != None and e4 != None and e5 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
            end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
            end4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e4)
            end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
            end = [(end1[0]+end2[0]+end5[0]+end4[0])/4,
                   (end1[1]+end2[1]+end5[1]+end4[1])/4,
                   (end1[2]+end2[2]+end5[2]+end4[2])/4]
            return end

    elif number == 2:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A1"]:    #
                e1 = sse["start"]
            if sse['label'] == mapping["A10"]:   #
                e2 = sse["end"]
            if sse['label'] == mapping["A9a"]:   #
                e4 = sse["start"]
            if sse['label'] == mapping["A4b"]:   #
                e5 = sse["end"]

        if e1 != None and e2 != None and e4 != None and e5 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
            end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
            end4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e4)
            end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
            end = [(end1[0]+end2[0]+end5[0]+end4[0])/4,
                   (end1[1]+end2[1]+end5[1]+end4[1])/4,
                   (end1[2]+end2[2]+end5[2]+end4[2])/4]
            return end

    elif number == 3:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A3"]:    #
                e1 = sse["start"]
            if sse['label'] == mapping["A5"]:   #
                e2 = sse["start"]
            if sse['label'] == mapping["A9a"]:   #
                e4 = sse["start"]
            if sse['label'] == mapping["A4b"]:   #
                e5 = sse["end"]

        if e1 != None and e2 != None and e4 != None and e5 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
            end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
            end4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e4)
            end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
            end = [(end1[0]+end2[0]+end5[0]+end4[0])/4,
                   (end1[1]+end2[1]+end5[1]+end4[1])/4,
                   (end1[2]+end2[2]+end5[2]+end4[2])/4]
            return end

    elif number == 4:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A3"]:    #
                e1 = sse["start"]
            if sse['label'] == mapping["A5"]:   #
                e2 = sse["start"]
            if sse['label'] == mapping["A6b"]:   #
                e5 = sse["end"]

        if e1 != None and e2 != None  and e5 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
            end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
            end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
            end = [(end1[0]+end2[0]+end5[0])/3,
                   (end1[1]+end2[1]+end5[1])/3,
                   (end1[2]+end2[2]+end5[2])/3]
            return end

    elif number == 5:
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A1"]:    #
                e1 = sse["start"]
            if sse['label'] == mapping["A7"]:    #
                e2 = sse["start"]
            if sse['label'] == mapping["A4b"]:   #
                e5 = sse["end"]

        if e1 != None and e2 != None and e5 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
            end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
            end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
            end = [(end1[0]+end2[0]+end5[0])/3,
                   (end1[1]+end2[1]+end5[1])/3,
                   (end1[2]+end2[2]+end5[2])/3]
            return end

    return None


def predict_PIN_ligand_position(annotation_dir, pcd, number):
    # at first get sses mapping to A-G (7 helices)
    model_annotation_file = f"{annotation_dir}/{MODEL[0]}{MODEL[1]}{MODEL[2]}-annotated.sses.json"
    mapping = sses_model_mapping_v7b1k1a(model_annotation_file)

    # read domain annotation
    annotation_file = f"{annotation_dir}/{pcd[0]}{pcd[1]}{pcd[2]}-annotated.sses.json"
    with open(annotation_file, "r") as fr:
        data = json.load(fr)

    # if there are helice C and use them to count start and end
    # read start and end residue number in annotation file
    e1, e2, e4, e5 = None, None, None, None
    if number == 1:
        # e.g.7wkw has a ligand between helice 2 and 9b and 4b
        for sse in data[f"{pcd[0]}{pcd[1]}{pcd[2]}"]["secondary_structure_elements"]:
            if sse['label'] == mapping["A9b"]:   #
                e1 = sse["start"]
            if sse['label'] == mapping["A4b"]:   #
                e2 = sse["start"]
            if sse['label'] == mapping["A2"]:    #
                e4 = sse["start"]
            if sse['label'] == mapping["A2"]:    #
                e5 = sse["end"]

        if e1 != None and e2 != None and e4 != None and e5 != None:
            # transform start and end residue to coordinate in the original (not rotated) cif file
            cif_file = f"{annotation_dir}/{pcd[0]}.cif"
            end1 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e1)
            end2 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e2)
            end4 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e4)
            end5 = utils_read_cif.find_res_avg_coord(cif_file, pcd[1], e5)
            end = [(end1[0]+end2[0]+end5[0]+end4[0])/4,
                   (end1[1]+end2[1]+end5[1]+end4[1])/4,
                   (end1[2]+end2[2]+end5[2]+end4[2])/4]
            return end

    return None



