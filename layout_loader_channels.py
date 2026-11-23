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
        elif format == "anicka":
            print("data.keys()", data.keys())
            if "Channels" not in list(data.keys()):
                return

            # for wanted chain_id we read all present channels
            counter = 0
            for chann_l2 in data["Channels"].keys():
              for chann in data["Channels"][chann_l2]:
                oldpos = None
                c = 0
                for i in chann["Profile"]:
                    if (oldpos == None) or (len3d(oldpos, [i['X'], i['Y'], i['Z']]) > 4) or (c > 20):
                        chan_inst_name = f'{chann["Id"]}{int(chann["Cavity"])}A{counter:04}'
                        self.names.append(chan_inst_name)
                        self.center[chan_inst_name] = [i['X'], i['Y'], i['Z']]
                        self.type[chan_inst_name] = chann['Type']
                        self.radius[chan_inst_name] = i['Radius']
                        self.id[chan_inst_name] = chann["Id"]
                        # chann["Id"] do not have to be set
                        if chann["Id"] in self.annotations:
                            self.text_names[chan_inst_name] = self.annotations[chann["Id"]]["Name"]
                        else:
                            self.text_names[chan_inst_name] = str(chann["Id"])
                        self.auth_chain[chan_inst_name] = ""
                        oldpos = self.center[chan_inst_name]
                        c = 0
                        counter += 1
                    c += 1

        else:
            if "Annotations" not in list(data.keys()):
                return

            # for wanted chain_id we read all present channels
            for chann in data["Annotations"]:
                self.a_id.append(chann["Id"])
                self.annotations[chann["Id"]] = chann
                if "Colour" in chann:
                    self.color[chan_inst_name] = chann["Colour"]
                else:
                    self.color[chan_inst_name] = None

            counter = 0
#        for chann_l2 in ["CSATunnels_MOLE", "CSATunnels_Caver", "TransmembranePores_MOLE", "TransmembranePores_Caver", "CofactorTunnels_MOLE", "CofactorTunnels_Caver"]:
#        ["CSATunnels_MOLE", "CSATunnels_Caver", "TransmembranePores_MOLE", "TransmembranePores_Caver", "CofactorTunnels_MOLE", "CofactorTunnels_Caver"]:

            for chann_l2 in data["Channels"].keys():
              for chann in data["Channels"][chann_l2]:
                oldpos = None
                c = 0
                for i in chann["Profile"]:
                    if (oldpos == None) or (len3d(oldpos, [i['X'], i['Y'], i['Z']]) > 4) or (c > 20):
                        chan_inst_name = f'{chann["Id"]}{int(chann["Cavity"])}A{counter:04}'
                        self.names.append(chan_inst_name)
                        self.center[chan_inst_name] = [i['X'], i['Y'], i['Z']]
                        self.type[chan_inst_name] = chann['Type']
                        self.radius[chan_inst_name] = i['Radius']
                        self.id[chan_inst_name] = chann["Id"]
                        # chann["Id"] do not have to be set
                        if chann["Id"] in self.annotations:
                            self.text_names[chan_inst_name] = self.annotations[chann["Id"]]["Name"]
                        else:
                            self.text_names[chan_inst_name] = str(chann["Id"])
                        self.auth_chain[chan_inst_name] = ""
                        oldpos = self.center[chan_inst_name]
                        c = 0
                        counter += 1
                    c += 1

        return
