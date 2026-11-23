#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path
import math
import json


def len3d(a, b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


class Channels(object):
    # init ligand structure, the
    # input parameter is ligand json filename and
    # chain_id - could be either chain_id ot metacaracter "ALL"
    # which means ligand from all chains
    def __init__(self, filename, chain_id, mol_online=False, format="anicka"):
        self.filename = filename
        self.label = os.path.basename(self.filename).split('.')[0]
        self.ch_layout = None

        # read input data structure
        with open(self.filename, 'r') as fd:
            data = json.load(fd)
        self.a_id = []   # all annotations identifiers
        self.annotations = {}
        self.names = []     # all annotations of each part of Channels (each fiveth member)
        self.text_names = {}
        self.type = {}
        self.center = {}
        self.id = {}
        self.auth_chain = {}
        self.radius = {}
        self.color = {}
        if mol_online:
            for chann in data:
                counter = 0
                oldpos = None
                c = 0
                for i in chann["Profile"]:
                    if (oldpos == None) or (len3d(oldpos, [i['X'], i['Y'], i['Z']]) > 4) or (c > 20):
                        chan_inst_name = f'{chann["Id"]}A{counter}'
                        self.names.append(chan_inst_name)
                        self.center[chan_inst_name] = [i['X'], i['Y'], i['Z']]
                        self.type[chan_inst_name] = chann['Type']
                        self.radius[chan_inst_name] = i['Radius']
                        self.id[chan_inst_name] = chann["Id"]
                        # chann["Id"] do not have to be set
                        self.text_names[chan_inst_name] = str(chann["Id"])
                        self.auth_chain[chan_inst_name] = str(chann["Id"])
                        if "Colour" in chann:
                            self.color[chan_inst_name] = chann["Colour"]
                        else:
                            self.color[chan_inst_name] = None
                        oldpos = self.center[chan_inst_name]
                        c = 0
                        counter += 1
                    c += 1
        elif format == "anicka" or True:
            if "Channels" not in list(data.keys()):
                return

            if "Annotations" in data and isinstance(data["Annotations"], list):
                for ann in data["Annotations"]:
                    aid = str(ann.get("Id", ""))
                    if aid:
                        self.a_id.append(aid)
                        if aid not in self.annotations or not self.annotations[aid].get("Name"):
                            self.annotations[aid] = ann

            cat_priority = [
                "ReviewedChannels_MOLE",
                "TransmembranePores_MOLE",
                "CofactorTunnels_MOLE",
                "CSATunnels_MOLE",
                "ReviewedChannels_Caver",
                "TransmembranePores_Caver",
                "CofactorTunnels_Caver",
                "CSATunnels_Caver",
            ]
            all_cats = list(data["Channels"].keys())
            all_cats.sort(key=lambda k: cat_priority.index(k) if k in cat_priority else 99)

            has_mole = any(
                "MOLE" in cat and len(data["Channels"].get(cat, [])) > 0
                for cat in all_cats
            )

            counter = 0
            t_idx = 0
            kept_tunnels = []
            used_ids = set()

            for chann_l2 in all_cats:
                if has_mole and "Caver" in chann_l2:
                    continue
                for chann in data["Channels"][chann_l2]:
                    prof = chann.get("Profile", [])
                    if len(prof) < 2:
                        continue
                    p_start = [prof[0]["X"], prof[0]["Y"], prof[0]["Z"]]
                    p_end = [prof[-1]["X"], prof[-1]["Y"], prof[-1]["Z"]]
                    cid = str(chann.get("Id", ""))
                    ctype = str(chann.get("Type", ""))

                    if ctype and ctype not in ("Tunnel", "Path", "None"):
                        t_name = ctype
                    elif cid in self.annotations and self.annotations[cid].get("Name"):
                        t_name = self.annotations[cid]["Name"]
                    elif "TransmembranePores" in chann_l2 or ctype == "Path":
                        t_name = "TM pore"
                    else:
                        t_name = str(cid)

                    is_dup = False
                    for kt in kept_tunnels:
                        if t_name == kt["name"] and not t_name.replace(".", "").isdigit():
                            is_dup = True
                            break
                        if len3d(p_end, kt["end"]) < 6.0 and len3d(p_start, kt["start"]) < 6.5:
                            is_dup = True
                            break
                    if is_dup:
                        continue

                    t_idx += 1
                    kept_tunnels.append({"name": t_name, "start": p_start, "end": p_end})
                    uniq_id = cid if cid not in used_ids else f"{cid}.{t_idx}"
                    used_ids.add(uniq_id)
                    cav = int(chann.get("Cavity", 0))
                    prefix = f"{t_idx:02d}.{cid}.{cav}"

                    oldpos = None
                    c = 0
                    for idx_p, i in enumerate(prof):
                        pt = [i["X"], i["Y"], i["Z"]]
                        is_last = idx_p == len(prof) - 1
                        if (oldpos is None) or (len3d(oldpos, pt) > 3.5) or is_last or (c > 18):
                            chan_inst_name = f"{prefix}A{counter:04d}"
                            self.names.append(chan_inst_name)
                            self.center[chan_inst_name] = pt
                            self.type[chan_inst_name] = ctype
                            self.radius[chan_inst_name] = i["Radius"]
                            self.id[chan_inst_name] = uniq_id
                            self.text_names[chan_inst_name] = t_name
                            self.auth_chain[chan_inst_name] = ""
                            if cid in self.annotations and "Colour" in self.annotations[cid]:
                                self.color[chan_inst_name] = self.annotations[cid]["Colour"]
                            else:
                                self.color[chan_inst_name] = None
                            oldpos = self.center[chan_inst_name]
                            c = 0
                            counter += 1
                        c += 1

        return
