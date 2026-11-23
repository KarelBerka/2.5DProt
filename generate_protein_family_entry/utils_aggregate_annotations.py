#!/usr/bin/env python3

import json


def aggregate_annotations(annotations_list, output_file_name, detect_dir, chain_id_mapping):
    # from the given annotations list which belonges to each chain of the input protein
    # and detect file {detect_dir}/{protein}-detected.sses.json
    # create aggregate annotaions
    final_annotation = {}
    if len(annotations_list) != 0:
        # start annotation is taken from the first chain annotation
        with open(annotations_list[0], "r") as fa:
            final_annotation = json.load(fa)
        protein = list(final_annotation.keys())[0]

        with open(f"{detect_dir}/{protein}-detected.sses.json", "r") as fd:
            data_detect = json.load(fd)

        transformation_table = {}

        # rename the annotation of the first af
        num = 0
        for sses in final_annotation[protein]["secondary_structure_elements"]:
            for sses_det in data_detect[protein]["secondary_structure_elements"]:
                if sses["chain_id"] == sses_det["chain_id"] and sses["start"] == sses_det["start"] and sses["end"] == sses_det["end"]:
                    sses["start_vector"] = sses_det["start_vector"]
                    sses["end_vector"] = sses_det["end_vector"]
                    sses["rainbow"] = "s"+str(int(100+800/len(annotations_list)*num))
                    sses["color"] = sses["rainbow"]
                    res = ""
                    for i in chain_id_mapping:
                        if chain_id_mapping[i][protein] == sses["chain_id"]:
                            res = i
                    transformation_table[sses["label"]] = sses["label"][0] + str(res) + sses["label"][1:]
                    sses["label"] = transformation_table[sses["label"]]

        # add beta_connactivity data as well
        new_beta_con = []
        for sses_pair in final_annotation[protein]["beta_connectivity"]:
            if sses_pair[0] in transformation_table and sses_pair[1] in transformation_table:
                new_beta_con.append([transformation_table[sses_pair[0]], transformation_table[sses_pair[1]], sses_pair[2]])

        final_annotation[protein]["beta_connectivity"] = new_beta_con

        # add another chain annotation to the af
        for annotation in annotations_list[1:]:
            num += 1
            with open(annotation, "r") as fa:
                data_add = json.load(fa)

            transformation_table = {}
            for sses in data_add[protein]["secondary_structure_elements"]:
                for sses_det in data_detect[protein]["secondary_structure_elements"]:
                    if sses["chain_id"] == sses_det["chain_id"] and sses["start"] == sses_det["start"] and sses["end"] == sses_det["end"]:
                        sses["start_vector"] = sses_det["start_vector"]
                        sses["end_vector"] = sses_det["end_vector"]
                        sses["rainbow"] = "s"+str(int(100+800/len(annotations_list)*num))
                        sses["color"] = sses["rainbow"]
                        res = ""
                        for i in chain_id_mapping:
                            if chain_id_mapping[i][protein] == sses["chain_id"]:
                                res = i
                        transformation_table[sses["label"]] = sses["label"][0] + str(res) + sses["label"][1:]

                        sses["label"] = sses["label"][0] + str(res)+sses["label"][1:]
                final_annotation[protein]["secondary_structure_elements"].append(sses)
                # add beta_connactivity data as well

            for sses_pair in data_add[protein]["beta_connectivity"]:
                if sses_pair[0] in transformation_table and sses_pair[1] in transformation_table:
                    final_annotation[protein]["beta_connectivity"].append([transformation_table[sses_pair[0]], transformation_table[sses_pair[1]], sses_pair[2]])

        data_detect[protein]["beta_connectivity"] = new_beta_con

    with open(output_file_name, "w") as fd:

        json.dump(final_annotation, fd, indent=4, separators=(",", ": "))

