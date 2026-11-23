#!/usr/bin/env python3

import numpy as np
import json
import logging
import time
import utils_common


def rot_trans_sync(rot_trans, start_set, both, copy):
    data = {}
    # find necessary rotation
    new_inv = np.linalg.inv(np.array(rot_trans[both][0]))
    transf_m = np.dot(start_set[both][0], new_inv)
    # transform all input layouts
    for i in rot_trans:
        rti = rot_trans[i]
        trans = [rti[1][0][0], rti[1][1][0]]
        rt = [np.dot(transf_m, np.array(rti[0])),
              np.dot(transf_m, [trans[0], trans[1]])]
        data[i] = ([[rt[0][0][0], rt[0][0][1]],
                   [rt[0][1][0], rt[0][1][1]]],
                   ([[rt[1][0]], [rt[1][1]]]))
    if copy is True:
        for i in start_set:
            if i not in data:
                data[i] = start_set[i]

    return data


# input are two rotation files the second contains one layout from the first
# one (the layout name is in parameter both) the result is a file, which
# contains rotation file for all layouts (from the first and second file as
# well) which are rotated based on the first rotation file if the third
# argument is "copy", copy the layouts from the first rot. File as well if
# there is the common layout, then do not copy them
def synchronize_rotation_files(first_rot, sec_rot, both, out, copy=False):
    with open(first_rot, "r") as fd:
        frt = json.load(fd)
    with open(sec_rot, "r") as fd:
        srt = json.load(fd)
    if both not in frt:
        if frt == {} and srt == {}:
            # diagrams are empty
            data = {}
        else:
            return -1
    else:
        data = rot_trans_sync(srt, frt, both, copy)

    with open(out, "w") as fd:
        json.dump(data, fd, indent=4, separators=(',', ': '))
    logging.info("INFO: rotation files %s, %s aggregated to %s" % (first_rot, sec_rot, out))
    return 0


# input are two layout files. Both of them contains differently named sses
# this script sync the name of the second one with the first one.
def sync_layout_sses_names(orig_fn, rename_fn, output_fn):
    logging.debug("rename sses from file %s to sses names from %s" % (rename_fn, orig_fn))

    with open(orig_fn, "r") as fd:
        of_data = json.load(fd)

    with open(rename_fn, "r") as fd:
        rf_data = json.load(fd)

    data = {}
    data["layout"] = {}
    data["size"] = {}
    data['helices_residues'] = {}
    data['angles'] = {}
    data['h_type'] = {}
    data['color'] = {}
    data['3dlayout'] = {}
    data['beta_connectivity'] = []  # have to be dealed saparately: bcc,
    data['generated'] = rf_data["generated"]
    data['protein'] = of_data["protein"]
    data['chain_id'] = of_data['chain_id']
    data['domain'] = of_data['domain']
    if 'ligands' in of_data:
        data['ligands'] = of_data['ligands']

    switch = {}
    for i in rf_data["layout"]:
        if rf_data["helices_residues"][i][0] > rf_data["helices_residues"][i][1]:
            rf_data["helices_residues"][i] = [rf_data["helices_residues"][i][1], rf_data["helices_residues"][i][0]]
        for j in of_data["layout"]:
            if of_data["helices_residues"][j][0] > of_data["helices_residues"][j][1]:
                of_data["helices_residues"][j] = [of_data["helices_residues"][j][1], of_data["helices_residues"][j][0]]
            if of_data["helices_residues"][j] == rf_data["helices_residues"][i]:
                switch[i] = j
    for i in rf_data['layout']:
        if i in switch:
            j = switch[i]
        else:
            j = i
        data['layout'][j] = rf_data['layout'][i]
        data['size'][j] = rf_data['size'][i]
        data['helices_residues'][j] = rf_data['helices_residues'][i]
        data['angles'][j] = rf_data['angles'][i]
        if i in rf_data['h_type']:
            data['h_type'][j] = rf_data['h_type'][i]
        else:
            data['h_type'][j] = 3
        if i in rf_data['color']:
            data['color'][j] = rf_data['color'][i]
        data['3dlayout'][j] = rf_data['3dlayout'][i]

    for i in rf_data['size']:
        if i not in data['layout']:
            j = i
            data['size'][j] = rf_data['size'][i]
            data['helices_residues'][j] = rf_data['helices_residues'][i]
            if i in rf_data['h_type']:
                data['h_type'][j] = rf_data['h_type'][i]
            if i in rf_data['color']:
                data['color'][j] = rf_data['color'][i]

    for i in rf_data['beta_connectivity']:
        s = i.split(":")
        t = ""
        for j in s:
            if j in switch:
                k = switch[j]
            else:
                k = j
            if (t == ""):
                t = str(k)
            else:
                t = t+":"+str(k)
        data['beta_connectivity'].append(t)

    fd = open(output_fn, "w")
    json.dump(data, fd, indent=4, separators=(',', ': '))
    fd.close()
    logging.info("INFO: rotation files aggregated to %s" % (output_fn))


def sync_layout_n(directory, filename, orig_fn):
    rename_fn = f"{directory}/2dprots/{filename}"
    output_fn = f"{directory}/{filename}"
    try:
        sync_layout_sses_names(orig_fn, rename_fn, output_fn)
    except Exception as err:
        command = f"sync_layout_sses_names({orig_fn}, {rename_fn}, {output_fn})"
        log_fn = f"{directory}/sync_layout_sses_names-fail-{time.time()}.txt"
        utils_common.error_report(directory, command, None, err, log_fn)
        logging.error(f"sync_layout_sses_names failed, report save to {log_fn}")
        return False
    return True


# input is a rotation file and a name of layout file, which should be removed
# output file
def cleanup_rotation_file(rotation, remove, remove2, out, out2):
    # remove should be removed and then remove2 substitute by remove
    logging.debug("INFO: We remove %s from rotation file %s and put it to %s" % (remove, rotation, out2))
    with open(rotation, "r") as fd:
        data = json.load(fd)

    d = {}
    for i in data:
        if (i != remove):
            d[i] = data[i]

    with open(out2, "w") as fd:
        json.dump(d, fd, indent=4, separators=(',', ': '))

    logging.debug("INFO: We substitute %s to %s in file %s and put it to %s" % (remove2, remove, out2, out))
    with open(out2, "r") as fr, open(out, "w") as fw:
        s = fr.readline()
        while s:
            ind = s.find(remove2)
            if ind != -1:
                t = s[:ind]+remove+s[ind+len(remove2):]
            else:
                t = s
            fw.write(t)
            s = fr.readline()
    logging.debug("rotation file cleanup to %s" % (out))
