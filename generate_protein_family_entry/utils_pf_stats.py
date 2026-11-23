#!/usr/bin/env python3

import logging
import json

# TODO: Code cleanup
# utils_stats should go through detected files not through pdb mappings api


def copy_and_rename_statistics(overprot_stats_list, family_sses_statistic):

    data = {}

    # open overprot data stats
    for overprot_stats in overprot_stats_list:
        try:
            with open(overprot_stats, "r") as fo:
                odata = json.load(fo)
        except FileNotFoundError as e:
            comment = f"Failed to load {overprot_stats} because of Exception ({e})"
            logging.error(comment)
            return False

        dir_split = overprot_stats.split("/")
        chain = dir_split[-2][-1]

        sses = odata["nodes"]
        for sse in sses:
            label = sse["label"][0] + str(chain) + sse["label"][1:]
            data[label] = {}
            data[label]["size"] = sse["avg_length"]
            data[label]["probability"] = sse["occurrence"]

    # save sses statistics
    try:
        with open(family_sses_statistic, "w") as fd:
            json.dump(data, fd, indent=4, separators=(',', ':'))
    except BaseException as e:
        comment = f"Failed to save {family_sses_statistic} because of Exception ({e})"
        logging.error(comment)
        return False

    logging.info("INFO: Statistics saved to %s" % family_sses_statistic)
    return True

