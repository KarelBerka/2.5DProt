#!/usr/bin/env python3

import logging
import json
import sys
sys.path.append('..')
import loader


def copy_statistic(directory, overprot_stats, family_sses_statistic):
    # copy overprot statistics to family_sses_statistic

    data = {}
    # open overprot data stats
    try:
        with open(overprot_stats, "r") as fo:
            odata = json.load(fo)
    except FileNotFoundError as e:
        comment = f"Failed to load {overprot_stats} because of Exception ({e})"
        logging.error(comment)
        return False

    # copy statistics to 2dprots structure
    sses = odata["nodes"]
    for sse in sses:
        data[sse["label"]] = {}
        data[sse["label"]]["size"] = sse["avg_length"]
        data[sse["label"]]["probability"] = sse["occurrence"]

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


def save_statistics(directory, d_list, dmax, family_sses_statistic, family_avg_protein):
    # variables to count average distances/sizes of sses
    # can be count from detected files
    sses_nontriv_occ_number = {}  # number of nontrivial helices
    size_sum = {}            # sum of all sizes per helic
    size_avg = {}            # average size per helic
    dist_count = 0           # number of all proteins in given family
    distances = {}
    occ = {}
    angles = {}

    err = []

    # go through all *detected* files
    for dom in d_list:
        f = f"{directory}/{dom[3]}.cif"
        prot = dom[0]
        chain = dom[1]
        domain = dom[2]
        logging.debug(f"loading {f} protein {prot} chain {chain} domain {domain}")
        try:
            protein = loader.Protein(f)
        except ValueError:
            err.append(dom)
            logging.error("Failed to load %s because of ValueError" % f)
            continue
        except Exception as e:
            err.append(dom)
            logging.error("Failed to load %s because of Exception (%s)" % (f, e))
            continue

        logging.debug("loading chain, domain")
        annotation_file = f"{directory}/{prot}-{chain}{domain}-annotated.sses.json"
        try:
            protein.add_chain_domain(chain, domain, annotation_file)
        except AssertionError as e:
            err.append(dom)
            logging.error("ERROR: Failed to load chain %s domain %s for %s" % (chain, domain, e))
            continue

        logging.debug("loading sses distances")
        try:
            dist = protein.get_distances(chain, domain)
        except AssertionError as e:
            err.append(dom)
            logging.error("ERROR: Unable to count distances for %s (%s)" % (f, e))
            continue
        except KeyError:
            err.append(dom)
            logging.error("Unable to find sec. structures for chain %s (file %s)"
                          % (chain, f))
            continue
        except TypeError as e:
            err.append(dom)
            logging.error("No secondary structures in chain %s (file %s): %s"
                          % (chain, f, e))
            continue

        try:
            sizes = protein.get_sizes_in_res(chain, domain)
        except AssertionError as e:
            err.append(protein)
            logging.error("Unable to count sizes for %s (%s)" % (f, e))
            continue
        except KeyError:
            err.append(protein)
            logging.error("Unable to find sec. structures for chain %s (file %s)"
                          % (chain, f))
            continue
        except TypeError as e:
            err.append(protein)
            logging.error("ERROR: No secondary structures in chain %s (file %s): %s"
                          % (chain, f, e))
            continue
        helices = protein.get_helices(chain, domain)
        ang = protein.get_all_helices_pairs_angles(chain, domain)

        # Count all secondary structures
        for h1 in helices:
            if h1 not in sses_nontriv_occ_number:
                size_sum[h1] = 0
                sses_nontriv_occ_number[h1] = 0
            size_sum[h1] += sizes[h1]
            sses_nontriv_occ_number[h1] += 1

        # Count average distances
        for h1 in dist:
            if h1 not in distances:
                distances[h1] = {}
                angles[h1] = {}
                occ[h1] = {}
            for h2 in dist[h1]:
                if h2 not in distances[h1]:
                    distances[h1][h2] = dist[h1][h2]
                    angles[h1][h2] = ang[h1, h2]
                    occ[h1][h2] = 1
                else:
                    distances[h1][h2] = distances[h1][h2]+dist[h1][h2]
                    angles[h1][h2] = angles[h1][h2]+ang[h1, h2]
                    occ[h1][h2] = occ[h1][h2]+1

        dist_count += 1
        if dist_count >= dmax and dmax != -1:
            break

    if sses_nontriv_occ_number == {} and err != []:
        # nothing to do no domain was sussesfully read
        return 0

    # find the sizes of helices
    for h in sses_nontriv_occ_number:
        size_avg[h] = float(size_sum[h]) / sses_nontriv_occ_number[h]

    # find average distance of sses
    for h1 in sses_nontriv_occ_number:
        if h1 not in sses_nontriv_occ_number:
            distances[h1] = {}
            angles[h1] = {}
        for h2 in sses_nontriv_occ_number:
            if (h1 in distances) and (h2 in distances[h1]):
                distances[h1][h2] = float(distances[h1][h2]) / occ[h1][h2]
                angles[h1][h2] = float(angles[h1][h2]) / occ[h1][h2]

    # Save average protein
    data = {}
    for i in sses_nontriv_occ_number:
        data[i] = {}
        data[i]["size"] = size_avg[i]
        data[i]["probability"] = float(sses_nontriv_occ_number[i])/dist_count
        data[i]["distances"] = distances[i]
        data[i]["angles"] = angles[i]
    fd = open(family_avg_protein, "w")
    json.dump(data, fd, indent=4, separators=(',', ':'))
    logging.info("INFO: Protein family saved to %s" % family_avg_protein)

    # Save sses statistics
    data = {}
    for i in sses_nontriv_occ_number:
        data[i] = {}
        data[i]["size"] = size_avg[i]
        data[i]["probability"] = float(sses_nontriv_occ_number[i])/dist_count
    fd = open(family_sses_statistic, "w")
    json.dump(data, fd, indent=4, separators=(',', ':'))
    logging.info("INFO: Statistics saved to %s" % family_sses_statistic)

    return err
