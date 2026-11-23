#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path
import logging
import math
import json
import datetime
from collections import OrderedDict
import utils
import utils_matrix2layout
import utils_ligand_matrix2layout
import utils_angleslayout
import utils_sses_types
import utils_cluster
import utils_palettes
import sys


class Ligands(object):
    # init ligand structure, the
    # input parameter is ligand json filename and
    # chain_id - could be either chain_id ot metacaracter "ALL"
    # which means ligand from all chains
    def __init__(self, filename, chain_id, mol_online=False):
        self.filename = filename
        self.label = os.path.basename(self.filename).split('.')[0]
        self.chain = chain_id
        self.l_layout = None

        # read input data structure
        with open(self.filename, 'r') as fd:
            data = json.load(fd)
        self.names = []
        self.atoms = {}
        self.center = {}
        self.e_type = {}
        self.ept = {}
        self.edsc = {}
        self.ecmp = {}
        self.auth_chain = {}
        self.residues = []
        self.res = {}

        if mol_online:
            # mol_online minimalistic style of input - show all 
            counter = 0
            for k in data:
                c = data[k]
                if (c["entity_comp"] == "HOH") or (c["entity_type"] == "polymer"):
                    continue
                self.names.append(counter)
                self.atoms[counter] = c["n_atoms"]
                self.center[counter] = c["center"]
                if "entity_type" not in c:
                    self.e_type[counter] = ""
                else:
                    self.e_type[counter] = c["entity_type"]
                self.edsc[counter] = c["entity_description"]
                self.ecmp[counter] = c["entity_comp"]
                if "res" not in c:
                    self.res[counter] = ""
                else:
                    self.res[counter] = c["res"]
#                self.chain[counter] = c["chain"]
                self.auth_chain[counter] = c["auth_chain"]
                counter += 2

        else:
        # for wanted chain_id we read all present ligands
          for c in data:
            if (self.chain == "ALL" or data[c]["auth_chain"] == chain_id):
                if (data[c]["entity_type"] != "water") and (data[c]["entity_polymer_type"] != "polypeptide(L)"):
                    if "residues" in data[c] and data[c]["residues"] is not None:
                        for res in data[c]["residues"]:
                            cc = c+str(res["residue"])
                            self.names.append(cc)
                            self.atoms[cc] = res["n_atoms"]
                            self.center[cc] = res["center"]
                            self.e_type[cc] = data[c]["entity_type"]
                            self.ept[cc] = data[c]["entity_polymer_type"]
                            self.edsc[cc] = data[c]["entity_description"]
                            self.ecmp[cc] = data[c]["entity_comp"]
                            self.auth_chain[cc] = data[c]["auth_chain"]
                    else:
                        self.names.append(c)
                        self.atoms[c] = data[c]["n_atoms"]
                        self.center[c] = data[c]["center"]
                        self.e_type[c] = data[c]["entity_type"]
                        self.ept[c] = data[c]["entity_polymer_type"]
                        self.edsc[c] = data[c]["entity_description"]
                        self.ecmp[c] = data[c]["entity_comp"]
                        self.auth_chain[c] = data[c]["auth_chain"]


class SSESStats(object):

    def __init__(self, filename):
        self.filename = filename
        self.label = os.path.basename(self.filename).split('.')[0]
        with open(self.filename, 'r') as fd:
            self._sses_stats = json.load(fd)
        self.label_to_type = self._determine_types()

    def most_probable(self, sses_set):
        prob_sse = 0
        prob = 0
        for sse in sses_set:
            if self._sses_stats[sse]["probability"] > prob:
                prob = self._sses_stats[sse]["probability"]
                prob_sse = sse
        return prob_sse

    def _determine_types(self):
        """
        Determine type (i.e. importance level) of secondary structure based on
        how often does it appear in given group for which this sses stats file
        was generated and how big is it in average in that group.
        """
        out = {}
        prime = 0
        for label in self._sses_stats:
           out[label] = utils_sses_types.determine_types(label, self._sses_stats[label]["probability"], self._sses_stats[label]["size"])
        return out

    def determine_types(self,label, size):
        if label in self.label_to_type:
            return self.label_to_type[label]
        else:
            if (utils_sses_types.is_nontrivial_sse(label,size[label])):
               return 2
            return 3

    def is_prime(self, label):
        if label in self.label_to_type:
            return self.label_to_type[label] == 1
        else:
            return False

    def is_sec(self, label, size):
        if label in self.label_to_type:
            return self.label_to_type[label] == 2
        else:
            if (utils_sses_types.is_nontrivial_sse(label,size[label])):
               return True
            return False

    def is_third(self, label):
        if label in self.label_to_type:
            return self.label_to_type[label] == 3
        else:
            return True

    def get_avg_size(self, label):
        return self._sses_stats[label]["size"]


class Annotation(object):

    def __init__(self, filename, chain, domain):
        self.filename = filename
        self.chain = chain
        self.domain = domain
        self.load()

    def load(self):
        with open(self.filename, 'r') as infile:
            self.data = json.load(infile)
        k = list(self.data.keys())
        if len(k) != 1:
            raise Exception("Annotation file has more than one structure")
        self.name = k[0]
        self.sses_data = self.data[self.name]["secondary_structure_elements"]
        self.beta_connectivity = self.data[self.name]["beta_connectivity"]
        self.beta_connectivity_next = {}
        self.beta_connectivity_prev = {}
        self.sses = []
        self.res_sizes = {}
        self.sizes = {}
        self.helices_residues = {}
        self.chain_id = {}
        self.vectors = {}
        self.color = {}
        self.m_distances = {}
        self.domains = {}
        for sse_number in range(len(self.sses_data)):
            self.sses.append(self.sses_data[sse_number]["label"])
            self.res_sizes[self.sses_data[sse_number]["label"]] = self.sses_data[sse_number]["end"] - self.sses_data[sse_number]["start"]
            self.helices_residues[self.sses_data[sse_number]["label"]] = [self.sses_data[sse_number]["end"], self.sses_data[sse_number]["start"]]
            self.vectors[self.sses_data[sse_number]["label"]] = [self.sses_data[sse_number]["start_vector"], self.sses_data[sse_number]["end_vector"]]
            self.sizes[self.sses_data[sse_number]["label"]] = math.dist(self.sses_data[sse_number]["start_vector"], self.sses_data[sse_number]["end_vector"])
            self.chain_id[self.sses_data[sse_number]["label"]] = self.sses_data[sse_number]["chain_id"]
            if "rainbow" in self.sses_data[sse_number]:
                self.color[self.sses_data[sse_number]["label"]] = self.sses_data[sse_number]["rainbow"]
            else:
                self.color[self.sses_data[sse_number]["label"]] = "grey"
        for item in self.beta_connectivity:
            if item[0] not in self.beta_connectivity_next:
                self.beta_connectivity_next[item[0]]=[]
            self.beta_connectivity_next[item[0]].append([item[1],item[2]])
            if item[1] not in self.beta_connectivity_next:
                self.beta_connectivity_next[item[1]]=[]
            self.beta_connectivity_next[item[1]].append([item[0],item[2]])

    def chain_sses(self, chain_set):
        ret = []
        for i in self.sses:
            if self.chain_id[i] in chain_set:
                ret.append(i)
        return ret

    
    def layout_3d(self):
        rl = {}
        for i in self.vectors:
                 rl_plain = self.vectors[i]
                 rl[i] = [list(rl_plain[0]), list(rl_plain[1])]
        return rl


    def not_in_chain_sses(self, chain_set):
        ret = []
        for i in self.sses:
            if self.chain_id[i] not in chain_set:
                ret.append(i)
        return ret

    def get_beta_connectivity_next(self, sheet):
        """
        return the sheets which are next in beta connections section
        """
        if sheet not in self.beta_connectivity_next.keys():
            return []
        else:
             return self.beta_connectivity_next[sheet]


    def get_beta_connectivity_prev(self, sheet):
        """
        return the sheets which are previous in beta connections section
        """
        if sheet not in self.beta_connectivity_prev.keys():
            return []
        else:
            return self.beta_connectivity_prev[sheet]


    def get_beta_connectivity_near(self, sheet):
        """
        return the sheets which are near (next or previous) in 
        beta connections section
        """
        if sheet not in self.beta_connectivity_next.keys():
            if sheet not in self.beta_connectivity_prev.keys():
                near = []
            else:
                near = self.beta_connectivity_prev[sheet]
        else:
            if sheet not in self.beta_connectivity_prev.keys():
                near = self.beta_connectivity_next[sheet]
            else:
               near = self.beta_connectivity_next[sheet] + self.beta_connectivity_prev[sheet]
        return near


    def get_beta_connectivity_cluster_orientation(self, start):
        """
        return cluster with decomposition to same/opposit orientation sets
        """
        # set with sheets with the same/opposit orientation
        sameset = [start]
        opposit = []

        # wariables which helps to track which sses have to be
        to_do = [start]                                           # parsed
        cluster = []                                              # are already parsed

        while to_do != []:
            sheet = to_do[0]
            nextsheet = self.get_beta_connectivity_next(sheet)
            # for all successors
            for i in nextsheet:
                if (i[0] not in cluster) and (i[0] not in to_do):
                    # if it was not processed already, add it to same/opposit list
                    if (i[1] == 1):
                        if (sheet in sameset):
                            sameset.append(i[0])
                        else:
                            opposit.append(i[0])
                    else:
                        if (sheet in sameset):
                            opposit.append(i[0])
                        else:
                            sameset.append(i[0])
                    to_do.append(i[0])

            to_do.remove(sheet)
            cluster.append(sheet)
        return (cluster, sameset, opposit)

    def get_beta_connectivity_clusters(self):
        """
        return the list of sheets which have the opposit orientation
        """
        output = []
        for h in self.sses:
            isnot = 0
            for i in output:
                if h in i:
                    isnot = 1
            # if h is not in any added cluster, add its cluster now
            if isnot == 0:
                set = self.get_beta_connectivity_cluster_orientation(h)
                output.append(":".join(set[0]))
        return output


    def get_distances(self):
        return utils.get_distance_matrix4n(self.vectors)

    def get_all_helices_pairs_angles(self):
        # get angles between all helic pairs

        # at first approximate helices with lines
        lines = {}
        for sse in self.sses:
            l = self.vectors[sse]
            lines[sse] = [l[0][0]-l[1][0],l[0][1]-l[1][1],l[0][2]-l[1][2]]

        # count angles between approximated lines
        angle = {}
        for h in self.sses:
            for g in self.sses:
                #  this is angle between h and g in degrees values 0..179
                if h != g:
                    a = (lines[h][0]*lines[g][0] + lines[h][1]*lines[g][1] + lines[h][2]*lines[g][2])/ (utils.vector_lenght(lines[h]) * utils.vector_lenght(lines[g]))
                    if a > 1:
                       a = 1
                    angle[h,g] = round(math.acos(a) * 180 / math.pi)
                else:
                    angle[h,g] = 0
        return angle


class Layout(object):

    COUNT_LAYOUT_MAX_GEN = 20000
    ###HELIC_TYPE_PRIM = 'prim'
    ###HELIC_TYPE_DIGIT = 'digit'
    ###HELIC_TYPE_DOT = 'dot'

    def __init__(self, sses_stats=None, ignore_ordering=None, annotation=None):
        self.ignore_ordering = ignore_ordering
        self.sses_stats = sses_stats
        self.annotation = annotation

    def read_sl_layout_from_json(self, filename):
        with open(filename, 'r') as fd:
            data = json.load(fd)
            self.layout=data["layout"]
            self.a_layout=data["layout"]
            self.h_type=data["h_type"]
            self.layout_3d=data["3dlayout"]

    def get_most_probable_sse(self, set):
       return self.sses_stats.most_probable(set)

    def get_prime_sses(self):
        prime = []
        for sse in self.annotation.sses:
            if self.sses_stats and self.sses_stats.is_prime(sse):
               prime.append(sse)
        return prime

    def get_prime_clusters(self,chain,domain, size):
        prime = []
        for cluster in self.annotation.get_beta_connectivity_clusters():
            out = 3
            for i in cluster.split(":"):
                if self.sses_stats and self.sses_stats.determine_types(i, size) < out:
                   out = self.sses_stats.determine_types(i,size)
            if (out == 1):
                prime.append(cluster)
        return prime

    def get_sec_sses(self,size):
        sec = []
        for sse in self.annotation.sses:
            if self.sses_stats and self.sses_stats.is_sec(sse, size):
               sec.append(sse)
        return sec

    def get_sec_clusters(self,chain,domain, size):
        sec = []
        for cluster in self.annotation.get_beta_connectivity_clusters():
            out = 3
            for i in cluster.split(":"):
                if self.sses_stats and self.sses_stats.determine_types(i, size) < out:
                   out = self.sses_stats.determine_types(i, size)
            if (out == 2):
                sec.append(cluster)
        return sec


    def get_prime2_clusters(self, chain, domain, size):
    # used for spike proteins - A, B, C are parts of trimer
        prime = []
        for cluster in self.annotation.get_beta_connectivity_clusters():
            out = 3
            inabc = True
            for i in cluster.split(":"):
                if self.annotation.chain_id[i] not in ("A", "B", "C"):
                   inabc = False
                if self.sses_stats and self.sses_stats.determine_types(i, size) < out:
                   out = self.sses_stats.determine_types(i,size)
            if (out == 1) and inabc:
                prime.append(cluster)
        return prime


    def get_sec2_clusters(self, chain, domain, size):
    # used for spike proteins - A, B, C are parts of trimer
        sec = []
        for cluster in self.annotation.get_beta_connectivity_clusters():
            out = 3
            inabc = False
            for i in cluster.split(":"):
                if self.annotation.chain_id[i] in ("A", "B", "C"):
                   inabc = True

                if self.sses_stats and self.sses_stats.determine_types(i, size) < out:
                   out = self.sses_stats.determine_types(i, size)
            if (out == 2) and inabc:
                sec.append(cluster)
        return sec


    def get_th2_clusters(self, chain, domain, size):
    # used for spike proteins - A, B, C are parts of trimer
        sec = []
        for cluster in self.annotation.get_beta_connectivity_clusters():
            inabc = False
            for i in cluster.split(":"):
                if self.annotation.chain_id[i] in ("A", "B", "C"):
                   inabc = True
            if not inabc:
                sec.append(cluster)
        return sec


    def get_third_clusters(self,chain,domain):
        third = []
        for cluster in self.annotation.get_beta_connectivity_clusters():
            out = 3
            for i in cluster.split(":"):
                if self.sses_stats and self.sses_stats.determine_types(i) < out:
                   out = self.sses_stats.determine_types(i)
            if (out == 3):
                third.append(cluster)
        return third


    def evolution_computation(self, prime_s, sec_s, prime_c, sec_c, all_c, \
        sses_sizes, sses_distances, sses_angles, start_layout, start_angles, \
        simple_clusters, near, orientation, sses_positions, chain, domain, tmp_dir):
        # We want to count layout for both prime and sec sses (third are not used)
        relevant_sses = prime_s
        relevant_sses.update(sec_s)

        # At first we have to compute clusters_layout
        if simple_clusters == True:
           # for start layout - not implemented
           self.c_layout, self.ca_layout = utils_cluster.simple_clusters_layout(all_c)
        else:
           self.c_layout, self.ca_layout = utils_cluster.compute_clusters_layout(
                            all_c, sses_distances, sses_sizes, near,
                            relevant_sses, orientation, start_layout,
                            start_angles, prime_s, chain, domain)

        # Now we are ready to do start layout if it is possible apply coords from starting layout
        self.layout = utils_matrix2layout.get_fake_layout(all_c, start_layout, relevant_sses, {})
        # Set layout angles to start version
        self.a_layout = utils_angleslayout.get_fake_angles_layout(all_c, start_angles, relevant_sses)

        # First round - optimize prime SSEs coords if we have enough of them
        relevant_clusters = prime_c   # primary SSEs only in this round
        relevant_sses = prime_s
        logging.debug(f"first round fixed are {set()} delam {relevant_clusters}")
        if len(relevant_clusters) > 0:
            self.layout, self.a_layout, self.c_layout, self.ca_layout = \
                    utils_matrix2layout.modify_layout_via_evolution(\
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, \
                    sses_distances,  sses_sizes, relevant_clusters, set(),
                    start_layout, tmp_dir)
            self.a_layout, self.ca_layout = \
                    utils_angleslayout.modify_angles_layout_via_evolution(\
                    self.a_layout, self.ca_layout, sses_angles,  sses_sizes, \
                    relevant_clusters, set(), start_angles, self.layout, \
                    self.c_layout, sses_positions)

            # Setup for a second round in which we do not want to move with what we have optimised now
            fixed_clusters = relevant_clusters
        else:
            fixed_clusters = set()

        # Lets set variables for next round
        relevant_clusters =  sec_c  # primary SSEs only in this round
        relevant_sses = sec_s
        logging.debug (f"second round fixed are {fixed_clusters} delam {relevant_clusters} layout {self.layout} angles {self.a_layout}")
        # Do second round
        if len(relevant_clusters) > 0:
            self.layout, self.a_layout, self.c_layout, self.ca_layout = \
                    utils_matrix2layout.modify_layout_via_evolution(\
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, \
                    sses_distances, sses_sizes, relevant_clusters, fixed_clusters, \
                    start_layout, tmp_dir)
            self.a_layout, self.ca_layout = \
                    utils_angleslayout.modify_angles_layout_via_evolution(\
                    self.a_layout, self.ca_layout,  sses_angles, sses_sizes, \
                    relevant_clusters, fixed_clusters, start_angles, \
                    self.layout, self.c_layout, sses_positions)

        # compute layout from cluster sfuff to per sse stuff
        self.layout, self.a_layout = \
                    utils_matrix2layout.count_sse_coords_based_on_cluster_coords( \
                    self.layout, self.a_layout, self.c_layout, self.ca_layout)


    def compute(self, chain, domain, tmp_dir, start_layout={}, start_angles={}, alpha_fold=False):
        self.chain = chain

        # Find out sses sizes, distances and angles in 3d (set of sses which we deal with in the newxt round)
        sses_sizes = self.annotation.res_sizes
        sses_distances = self.annotation.get_distances()
        sses_angles = self.annotation.get_all_helices_pairs_angles()

        # Get list of all/prime/sec  SSEs/clusters we are going to show in this protei
        if (self.sses_stats):
            prime_s = set(self.get_prime_sses())
            sec_s = set(self.get_sec_sses(sses_sizes))
            all_s = set(self.annotation.sses)
        else:
            prime_s = set()
            sec_s = set(self.annotation.sses)
            all_s = set(self.annotation.sses)

        self.h_type = {}
        for s in all_s:
           if s in prime_s:
              self.h_type[s] = 1
           else:
              if s in sec_s:
                 self.h_type[s] = 2
              else:
                 self.h_type[s] = 3

        if self.sses_stats:
            prime_c = set(self.get_prime_clusters(chain, domain, sses_sizes))
            sec_c = set(self.get_sec_clusters(chain, domain, sses_sizes))
            all_c = self.annotation.get_beta_connectivity_clusters()
        else:
            prime_c = set()
            sec_c = set(self.annotation.get_beta_connectivity_clusters())
            all_c = self.annotation.get_beta_connectivity_clusters()

        # Get list of 3d positions of all sses
        sses_positions = self.annotation.vectors

        self.evolution_computation(prime_s, sec_s, prime_c, sec_c, all_c, sses_sizes, sses_distances, sses_angles, \
                        start_layout, start_angles, False, self.annotation.get_beta_connectivity_near, \
                        self.annotation.get_beta_connectivity_cluster_orientation, sses_positions, chain, domain, \
                        tmp_dir)


    def three_steps_evolution_computation(self, prime_s, sec_s, th_s, prime_c, sec_c, all_c, th_c, \
        sses_sizes, sses_distances, sses_angles, start_layout, start_angles, \
        simple_clusters, near, orientation, sses_positions, chain, domain):
        # We want to count layout for both prime and sec sses (third are not used)

        relevant_sses = set()
        relevant_sses.update(prime_s)
        relevant_sses.update(sec_s)
        relevant_sses.update(th_s)

        # At first we have to compute clusters_layout
        if simple_clusters == True:
           # for start layout - not implemented
           self.c_layout, self.ca_layout = utils_cluster.simple_clusters_layout(all_c)
        else:
           self.c_layout, self.ca_layout = utils_cluster.compute_clusters_layout(
                            all_c, sses_distances, sses_sizes, near,
                            relevant_sses, orientation, start_layout,
                            start_angles, prime_s, chain, domain)

        # Now we are ready to do start layout if it is possible apply coords from starting layout
        self.layout = utils_matrix2layout.get_fake_layout(all_c, start_layout, relevant_sses, OrderedDict())
        # Set layout angles to start version
        self.a_layout = utils_angleslayout.get_fake_angles_layout(all_c, start_angles, relevant_sses)

        # First round - optimize prime SSEs coords if we have enough of them
        fixed_clusters = set()
        relevant_clusters = prime_c   # primary SSEs only in this round

        logging.debug(f"First round fixed are {set()} delam {relevant_clusters}")
        if len(relevant_clusters) > 0:
            self.layout, self.a_layout, self.c_layout, self.ca_layout = \
                    utils_matrix2layout.modify_layout_via_evolution(\
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, \
                    sses_distances,  sses_sizes, relevant_clusters, set(), \
                    start_layout, flags={"step_size" : 10, "dbs" : "spike-protein"})
            self.a_layout, self.ca_layout = \
                    utils_angleslayout.modify_angles_layout_via_evolution(\
                    self.a_layout, self.ca_layout, sses_angles,  sses_sizes, \
                    relevant_clusters, fixed_clusters, start_angles, self.layout, \
                    self.c_layout, sses_positions)

        # Setup for a second round in which we do not want to move with what we have optimised now
        fixed_clusters = relevant_clusters
        relevant_clusters = sec_c   # primary SSEs only in this round
        # set relevant clusters to the position of the nearest fixed sse
        if start_layout == {}:
            for cl in relevant_clusters:
                sse1 = cl.split(":")[0]
                dist = 10000
                for cl2 in fixed_clusters:
                    for sse2 in cl2.split(":"):
                        if sses_distances[sse1][sse2] < dist:
                            self.layout[cl] = self.layout[cl2]
                            dist = sses_distances[sse1][sse2]

        # Second round - optimize prime SSEs coords if we have enough of them
        logging.debug(f"Second round fixed are {fixed_clusters} delam {relevant_clusters}")

        if len(relevant_clusters) > 0:
            self.layout, self.a_layout, self.c_layout, self.ca_layout = \
                    utils_matrix2layout.modify_layout_via_evolution(\
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, \
                    sses_distances,  sses_sizes, relevant_clusters, set(), \
                    start_layout, flags={"step_size" : 10, "dbs" : "spike-protein"})
            self.a_layout, self.ca_layout = \
                    utils_angleslayout.modify_angles_layout_via_evolution(\
                    self.a_layout, self.ca_layout, sses_angles,  sses_sizes, \
                    relevant_clusters, set(), start_angles, self.layout, \
                    self.c_layout, sses_positions)

        # Lets set variables for next round
        fixed_clusters = fixed_clusters.union(relevant_clusters)
        relevant_clusters =  th_c  # primary SSEs only in this round
        if start_layout == {}:
            for cl in relevant_clusters:
                sse1 = cl.split(":")[0]
                dist = 10000
                for cl2 in fixed_clusters:
                    for sse2 in cl2.split(":"):
                        if sses_distances[sse1][sse2] < dist:
                            self.layout[cl] = self.layout[cl2]
                            dist = sses_distances[sse1][sse2]

        logging.debug (f"Third round fixed are {fixed_clusters} delam {relevant_clusters}")
        # Do third round
        if len(relevant_clusters) > 0:
            self.layout, self.a_layout, self.c_layout, self.ca_layout = \
                    utils_matrix2layout.modify_layout_via_evolution(\
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, \
                    sses_distances, sses_sizes, relevant_clusters, fixed_clusters,
                    start_layout, flags={"step_size" : 20, "dbs" : "spike-protein"})
            self.a_layout, self.ca_layout = \
                    utils_angleslayout.modify_angles_layout_via_evolution(\
                    self.a_layout, self.ca_layout,  sses_angles, sses_sizes, \
                    relevant_clusters, fixed_clusters, start_angles, \
                    self.layout, self.c_layout, sses_positions)

        # compute layout from cluster sfuff to per sse stuff
        self.layout, self.a_layout = \
                    utils_matrix2layout.count_sse_coords_based_on_cluster_coords( \
                    self.layout, self.a_layout, self.c_layout, self.ca_layout)


    def compute_whole_spike_protein(self, chain, domain, start_layout={}, start_angles={}, alpha_fold=False):
        self.chain = chain

        # Find out sses sizes, distances and angles in 3d (set of sses which we deal with in the newxt round)
        sses_sizes = self.annotation.res_sizes
        sses_distances = self.annotation.get_distances()
        sses_angles = self.annotation.get_all_helices_pairs_angles()

        # Get list of all/prime/sec  SSEs/clusters we are going to show in this protei
        trimer = set(self.annotation.chain_sses(("A", "B", "C")))
        prime_s = trimer.intersection(set(self.get_prime_sses()))
        sec_s = trimer.intersection(set(self.get_sec_sses(sses_sizes)))
        th_s = set(self.annotation.not_in_chain_sses(("A", "B", "C")))
        all_s = prime_s.union(sec_s,th_s)

        self.h_type = {}
        for s in all_s:
           if s in prime_s:
              self.h_type[s] = 1
           else:
              if (s in sec_s) or (s in th_s):
                 self.h_type[s] = 2
              else:
                 self.h_type[s] = 3

        prime_c = set(self.get_prime2_clusters(chain, domain, sses_sizes))
        sec_c = set(self.get_sec2_clusters(chain, domain, sses_sizes))
        th_c = set(self.get_th2_clusters(chain, domain, sses_sizes))
        all_c = self.annotation.get_beta_connectivity_clusters()

        # Get list of 3d positions of all sses
        sses_positions = self.annotation.vectors

        self.three_steps_evolution_computation(prime_s, sec_s, th_s, prime_c, sec_c, all_c, th_c, sses_sizes, sses_distances, sses_angles, \
                        start_layout, start_angles, False, self.annotation.get_beta_connectivity_near, \
                        self.annotation.get_beta_connectivity_cluster_orientation, sses_positions, chain, domain)


    def compute_whole_protein(self, layouts, chains, domains, tmp_dir):
        # TODO - propagate h_type to inner fuctions
        # Find out sses, sizes, distances in 3d
        sses_sizes = {}      # sizes of all sses in all domains, clusters
        sses_rsizes = {}
        sses_distances = {}  # distances of all sses in all domains, cluster pairs
        number = 0

        for c1 in chains:
           number += len(domains[c1])

        if number > 70:
            short = True
        else:
            short = False

        d_layout = {}
        for c1 in chains:
           sses_sizes[c1] = {}
           sses_rsizes[c1] = {}
           sses_distances[c1] = {}
           d_layout[c1] = {}
           for d1 in domains[c1]:
              sses_sizes[c1][d1] = {}
              sses_rsizes[c1][d1] = {}
              d_layout[c1][d1] = {}

              for i in layouts["layout"][c1][d1]:
                  if i in layouts["hres"][c1][d1]:
                      sses_rsizes[c1][d1][i] = layouts["hres"][c1][d1][i][1] - layouts["hres"][c1][d1][i][0]
                      for j in self.annotation.helices_residues:
                           if self.annotation.helices_residues[j][0] == layouts["hres"][c1][d1][i][0] and self.annotation.helices_residues[j][0] == layouts["hres"][c1][d1][i][0] and layouts["chain_id"][c1][d1][i] == self.annotation.chain_id[j]:
                               d_layout[c1][d1][i] = self.annotation.vectors[j]
                               layouts["layout_3d"][c1][d1][i] = self.annotation.vectors[j]
                      sses_sizes[c1][d1][i] = math.dist(layouts["layout_3d"][c1][d1][i][0], layouts["layout_3d"][c1][d1][i][1])
              sses_distances[c1][d1] = {}

        for c1 in chains:
           for d1 in domains[c1]:
              for c2 in chains:
                 sses_distances[c1][d1][c2] = {}
                 for d2 in domains[c2]:
                     sses_distances[c1][d1][c2][d2] = utils.get_distance_matrix4n_2(d_layout[c1][d1], d_layout[c2][d2])
 
        self.sizes = sses_sizes
        self.rsizes = sses_rsizes

        # get fake layout and angles layout
        # Now we are ready to do start layout if it is possible apply coords from starting layout
        self.c_layout = layouts["layout"]
        self.ca_layout = layouts["a_layouts"]
        self.layout = {}
        self.a_layout = {}
        for c1 in chains:
           self.layout[c1] = {}
           self.a_layout[c1] = {}
           for d1 in domains[c1]:
              self.layout[c1][d1] = [0,0]
              self.a_layout[c1][d1] = 0

        if len(chains) > 0:
            self.layout, self.a_layout, self.c_layout, self.ca_layout = \
                    utils_matrix2layout.modify_protein_layout_via_evolution( \
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, \
                    sses_distances, sses_sizes, tmp_dir=tmp_dir)

        # compute layout from cluster sfuff to per sse stuff
        self.layout, self.a_layout = \
                    utils_matrix2layout.count_sse_protein_coords_based_on_cluster_coords( \
                    self.layout, self.a_layout, self.c_layout, self.ca_layout, 1)


    def compute_avg(self, p_data, c_data,  sse_stats):
#        generations = 40000
        start_layout = {}
        start_angles = {}

        # Find out sses sizes, distances and angles in 3d (set of sses which we deal with in the newxt round)
        sses_sizes = {}
        sses_distances = {}
        sses_angles = {}
        for i in p_data:
            sses_sizes[i] = p_data["size"]
            sses_distances[i] = {}
            sses_angles[i] = {}
            for j in p_data[i]["distances"]:
                sses_distances[i][j] = p_data[i]["distances"][j]
                sses_angles[i][j] = p_data[i]["angles"][j]

        sses_s = SSESStats(sse_stats)
        # Get list of all/prime/sec  SSEs/clusters we are going to show in this protei
        all_c = []
        for i in c_data:
           if c_data[i] > 0.5:
              all_c.append(i)

        prime_c = []
        sec_c = []
        for i in all_c:
           state = 3
           for j in i.split(":"):
              if sses_s.determine_types(i) < state:
                 state = sses_s.determine_types(i)
           if state == 1:
              prime_c.append(i)
           if state == 2:
              sec_c.append(i)

        all_s = []
        for i in all_c:
           for j in i.split(":"):
              all_s.append(j)

        prime_s = []
        sec_s = []
        for i in all_s:
           stype = sses_s.determine_types(i)
           if stype == 1:
              prime_s.append(i)
           if stype == 2:
              sec_s.append(i)

        self.evolution_computation(prime_s, sec_s, prime_c, sec_c, all_c, sses_sizes, sses_distances, sses_angles,
                                   start_layout, start_angles, True, False, False)

    def compute_ligands(self, ligands, gcl=False, layout_3d=False, near=False):
        ligands_positions = ligands.center             # Get list of ligands 3d=positions

        # Get list of all SSEs we are going to show in this protein
        if gcl:
            # whole protein
            sses_positions = {}
            sses_layout = {}
            sses_sizes = {}
            sses_chain_id = {}
            for c in layout_3d:
                for d in layout_3d[c]:
                    for sse in layout_3d[c][d]:
                        if sse in self.layout[c][d]:
                            sses_positions[sse+c+d] = layout_3d[c][d][sse]
                            sses_layout[sse+c+d] = self.layout[c][d][sse]
                            sses_sizes[sse+c+d] = self.sizes[c][d][sse]
                            sses_chain_id[sse+c+d] = c
            all_s = sses_positions.keys()
            prime_s = set()                          # TODO - should be set as well
        else:

            all_s = set()
            sses_positions = self.annotation.vectors  # Get list of sses 3d positions
            sses_layout = self.layout
            sses_sizes = self.annotation.res_sizes
            sses_chain_id = {}
            for sse in self.layout:
                 sses_chain_id[sse] = ligands.chain
                 all_s.add(sse)
            if (self.sses_stats):
                 prime_s = set(self.get_prime_sses())
            else:
                 prime_s = set()

        if near:
            # if we create per domain diagram, we have to choose ligands which "belonges"
            # to given domain manualy - we set the minimal distance from domain sses
            # count the minimal distance from each ligand to sses TODO
            near_ligands = []
            for l1 in ligands.names:
                for sse in all_s:
                     line_start = sses_positions[sse][0]
                     line_diff = [sses_positions[sse][1][0] - sses_positions[sse][0][0],
                                  sses_positions[sse][1][1] - sses_positions[sse][0][1],
                                  sses_positions[sse][1][2] - sses_positions[sse][0][2]]
                     l1_sse_dist = int(utils.point_to_line_dist(line_start, line_diff, ligands_positions[l1]))
                     if l1_sse_dist < 15:                                              # distance which denotes ligand is near
                         near_ligands.append(l1)
                         break
            ligands.names = near_ligands                                               # Get lis of all ligands
            ligands.atoms = {key: ligands.atoms[key] for key in ligands.names}         # Get list of ligands sizes
            ligands_positions = {key: ligands.center[key] for key in ligands.names}    # Get list of ligands 3d positions
        # if the minimal distance is > 10 the ligand will not be shown - move it from the ligand lists TODO

        # count distances
        ligands_sses_distances = {}
        for l1 in ligands.names:
            ligands_sses_distances[l1] = {}
            for l2 in ligands.names:
                 ligands_sses_distances[l1][l2] = int(math.dist(ligands_positions[l1], ligands_positions[l2]))
            for sse in all_s:
                 line_start = sses_positions[sse][0]
                 line_diff = [sses_positions[sse][1][0] - sses_positions[sse][0][0],
                              sses_positions[sse][1][1] - sses_positions[sse][0][1],
                              sses_positions[sse][1][2] - sses_positions[sse][0][2]]
                 ligands_sses_distances[l1][sse] = int(utils.point_to_line_dist(line_start, line_diff, ligands_positions[l1]))

        # Now we are ready to create start layout
        ligands.layout = {}
        for l in ligands.names:
             ligands.layout[l] = [0, 0]

        # All ligands are computed in one bunch
        logging.debug(f"finding ligands layout for protein ligands {ligands.names}")

        if len(ligands.names) > 0:
            ligands_layout = utils_ligand_matrix2layout.modify_ligand_layout_via_evolution(
                ligands, all_s, prime_s, sses_sizes,
                sses_layout, sses_chain_id, ligands_sses_distances,
                layout_3d, flag="ligands", tmp_dir="")
            self.l_layout = ligands_layout
        else:
            self.l_layout = {}


    def compute_channels(self, channels, tmp_dir, layout_3d=False, input_sizes=None, flag=""):
        channels_positions = channels.center             # Get list of ligands 3d=positions

        # Get list of all SSEs we are going to show in this protein
        # whole protein
        sses_positions = {}
        sses_layout = {}
        sses_sizes = {}
        sses_chain_id = {}

        if (flag == "domain") or (flag =="add_channels"):
            # there is no need to have domain and chain
            for sse in layout_3d:
                if sse in self.layout:
                    sses_positions[sse] = layout_3d[sse]
                    sses_layout[sse] = self.layout[sse]
                    if input_sizes != None: # used in generate_scl_layout
                        sses_sizes[sse] = input_sizes[sse]
                    else:
                        sses_sizes[sse] = self.annotation.sizes[sse]
                    sses_chain_id[sse] = "ALL"
            layout_3d = {'ALL':{'ALL': layout_3d}}
        else:
            # protein uses domain and chain
            for c in layout_3d:
                for d in layout_3d[c]:
                    for sse in layout_3d[c][d]:
                        if sse in self.layout[c][d]:
                            sses_positions[sse+c+d] = layout_3d[c][d][sse]
                            sses_layout[sse+c+d] = self.layout[c][d][sse]
                            sses_sizes[sse+c+d] = self.sizes[c][d][sse]
                            sses_chain_id[sse+c+d] = c
        all_s = sses_positions.keys()
        prime_s = set()                          # TODO - should be set as well

        # count distances
        channels_sses_distances = {}
        for ch1 in channels.names:
            channels_sses_distances[ch1] = {}
            for sse in all_s:
                 line_start = sses_positions[sse][0]
                 line_diff = [sses_positions[sse][1][0] - sses_positions[sse][0][0],
                              sses_positions[sse][1][1] - sses_positions[sse][0][1],
                              sses_positions[sse][1][2] - sses_positions[sse][0][2]]
                 channels_sses_distances[ch1][sse] = int(utils.point_to_line_dist(line_start, line_diff, channels_positions[ch1]))

        aux_names = []
        aux_channels_sses_distances = {}
        if flag == "domain":
            # filter out all chain members which are more then 7A from any sses
            for ch1 in channels.names:
                min_dist = 100
                for sse in all_s:
                    if channels_sses_distances[ch1][sse] < min_dist:
                        min_dist = channels_sses_distances[ch1][sse]
                if min_dist < 10 or flag == "protein":
                    aux_names.append(ch1)
            channels.names = aux_names

        # count distances
        channels_sses_distances = {}
        for ch1 in channels.names:
            channels_sses_distances[ch1] = {}
            for ch2 in channels.names:
                 channels_sses_distances[ch1][ch2] = int(math.dist(channels_positions[ch1], channels_positions[ch2]))
            for sse in all_s:
                 line_start = sses_positions[sse][0]
                 line_diff = [sses_positions[sse][1][0] - sses_positions[sse][0][0],
                              sses_positions[sse][1][1] - sses_positions[sse][0][1],
                              sses_positions[sse][1][2] - sses_positions[sse][0][2]]
                 channels_sses_distances[ch1][sse] = int(utils.point_to_line_dist(line_start, line_diff, channels_positions[ch1]))

        # Now we are ready to create start layout
        channels.layout = {}
        for ch in channels.names:
             channels.layout[ch] = [0, 0]

        # All ligands are computed in one bunch
        logging.debug(f"finding ligands layout for protein ligands {channels.names}")
        if len(channels.names) > 0:
            channels_layout = utils_ligand_matrix2layout.modify_ligand_layout_via_evolution(
                channels, all_s, prime_s, sses_sizes,
                sses_layout, sses_chain_id, channels_sses_distances,
                layout_3d, flag="channels", tmp_dir=tmp_dir,
                channels_positions=channels_positions,
                sses_positions=sses_positions)
            self.ch_layout = channels_layout
        else:
            self.ch_layout = {}


    def save_stats(self, chain, domain, statsfn=None):
        if statsfn is None:
            return False

        name = "layout-%s-%s%s.json" % (self.annotation.name, chain, domain)
        data = {}
        sizese = self.annotation.helices_residues
        size = {}
        for i in sizese:
            size[i] = sizese[i][1] - sizese[i][0] + 1
        data[name] = {"helices": self.annotation.helices_residues.keys(),
                      "prime": self.get_prime_sses(),
                      "sec": self.get_sec_sses(size),
                      "size": size, }
        with open(statsfn, 'a') as fd:
            json.dump(data, fd, indent=4, separators=(',', ': '))
        return name

    def save(self, chain, domain, filename, ligands=None, channels=None):
        if len(self.layout) == 0:
            print("WARN: Saving layout which is empty!")

        rl = {}
        for i in self.layout:
            if i in self.annotation.vectors:
                 rl_plain = self.annotation.vectors[i]
                 rl[i] = [list(rl_plain[0]), list(rl_plain[1])]

        data = {
            'layout': self.layout,
            'helices_residues': self.annotation.helices_residues,
            'size': self.annotation.sizes,
            'protein': self.annotation.name,
            'chain_id': self.annotation.chain_id,
            'chain': chain,
            'domain': domain,
            'filename_protein': self.annotation.filename,
            'generated': datetime.datetime.now().isoformat(),
            'angles': self.a_layout,
            'h_type': self.h_type,
            'beta_connectivity': self.annotation.get_beta_connectivity_clusters(),
            'color': self.annotation.color,
            '3dlayout': rl,
        }

        if ligands:
            l_data = {}
            if "l_layout" in dir(self):
                for l in self.l_layout.keys():
                    l_data[l] = {
                        'layout': self.l_layout[l],
                        'atoms': ligands.atoms[l],
                        'e_type': ligands.e_type[l],
                        'ept': ligands.ept[l],
                        'edsc': ligands.edsc[l],
                        'ecmp': ligands.ecmp[l], }
            data['ligands'] = l_data

        if channels:
            ch_data = {}
            for ch in channels.names:
                ch_data[ch] = {
                    'layout': self.ch_layout[ch],
                    'radius': channels.radius[ch],
                    'id': channels.id[ch],
                    'text_names': channels.text_names[ch],
                    '3d_layout': channels.center[ch],
                }
            data['channels'] = ch_data

        with open(filename, 'w') as fd:
            json.dump(data, fd, indent=4, separators=(',', ': '))
        return filename

    def save_protein(self, chains, domains, h_type, bcco, hres, filename, prot_name, layout_3d, ligands, channels, mol_online=False):
        # count the number of cdomains in whole protein (necessary to code sses names)
        number = 0
        maxhr = 0
        for c in chains:
           for d in domains[c]:
              number = number + 1
              for sses in hres[c][d]:
                 if hres[c][d][sses][1] > maxhr:
                    maxhr = hres[c][d][sses][1]
        maxhr = maxhr+1
        layout = {}
        a_layout = {}
        size = {}
        bcc = []
        hr = {}
        color = {}
        layout3d = {}
        ch = {}
        dom = {}
        ht = {}
        n = 0
        nn = 0
        for c in chains:
            ncolor = utils_palettes.set_domain_color(c)
            for d in domains[c]:
                for sses in self.layout[c][d]:
                    if sses not in self.sizes[c][d]:
                        continue
                    if (number < 11):
                        p = "%01d" % (n)
                    elif (number < 100):
                        p = "%02d" % (n)
                    elif (number < 1000):
                        p = "%03d" % (n)
                    elif (number < 10000):
                        p = "%04d" % (n)
                    else:
                        sys.exit()
                    if sses[0] == "_":                       # minor SSEs starts _E or _H - thus skip first two chars
                        name = str(sses[:2]+p+sses[2:])
                    else:                                  # major SSEs starts with E or H - thus skip only one character
                        name = str(sses[0]+p+sses[1:])
                    layout[name] = self.layout[c][d][sses]
                    a_layout[name] = self.a_layout[c][d][sses]
                    size[name] = self.sizes[c][d][sses]
                    hr[name] = (hres[c][d][sses][0], hres[c][d][sses][1])
                    ch[name] = c
                    dom[name] = d
                    color[name] = ncolor
                    layout3d[name] = layout_3d[c][d][sses]
                    if sses in h_type[c][d]:
                        ht[name] = h_type[c][d][sses]
                    else:
                        ht[name] = 3
                for i in bcco[c][d]:
                    bccj = {}
                    isplit = i.split(":")
                    for sses in isplit:
                        if (number < 11):
                            p = "%01d" % (n)
                        elif (number < 100):
                            p = "%02d" % (n)
                        elif (number < 1000):
                            p = "%03d" % (n)
                        elif (number < 10000):
                            p = "%04d" % (n)
                        else:
                            # TODO: fix this
                            sys.exit()

                        name = str(sses[0]+str(n)+sses[1:])
                        if bccj == {}:
                            bccj = name
                        else:
                            bccj = bccj+":"+name
                    bcc.append(bccj)
                n = n + 1
            nn = nn + 1
        if len(self.layout) == 0:
            print("WARN: Saving layout which is empty!")
        data = {
            'layout': layout,
            'size': size,
            'helices_residues': hr,
            'generated': datetime.datetime.now().isoformat(),
            'angles': a_layout,
            'h_type': ht,
            'beta_connectivity': bcc,
            'protein': "protein " + prot_name,
            'chain_id': ch,
            'color': color,
            'domain': dom,
            '3dlayout': layout3d,
        }
        if ligands:
            l_data = {}
            if mol_online:
                for l in ligands.names:
                    l_data[l] = {
                       'layout': self.l_layout[l],
                        'atoms': ligands.atoms[l],
                        'e_type': ligands.e_type[l],
#                        'ept': ligands.ept[l],
                        'chain': ligands.auth_chain[l],
                        'auth_chain': ligands.auth_chain[l],
                        'edsc': ligands.edsc[l],
                        'res': ligands.res[l],
                        'ecmp': ligands.ecmp[l], }
            else:
                for l in ligands.names:
                    l_data[l] = {
                       'layout': self.l_layout[l],
                        'atoms': ligands.atoms[l],
                        'e_type': ligands.e_type[l],
                        'ept': ligands.ept[l],
                        'edsc': ligands.edsc[l],
                        'ecmp': ligands.ecmp[l], }
            data['ligands'] = l_data

        if channels:
            ch_data = {}
            for ch in channels.names:
                if ch in channels.color:
                   c = channels.color[ch]
                else:
                   c = None
                ch_data[ch] = {
                    'layout': self.ch_layout[ch],
                    'radius': channels.radius[ch],
                    'id': channels.id[ch],
                    'text_names': channels.text_names[ch],
                    '3d_layout': channels.center[ch],
                    'color': c,
                }
            data['channels'] = ch_data

        with open(filename, 'w') as fd:
            json.dump(data, fd, indent=4, separators=(',', ': '))
        return filename

    def add_channels_to_protein(self, orig_filename, dest_filename, channels):
        with open(orig_filename, 'r') as infile:
            data = json.load(infile)

        if channels:
            ch_data = {}
            for ch in channels.names:
                if ch in channels.color:
                   c = channels.color[ch]
                else:
                   c = None
                ch_data[ch] = {
                    'layout': self.ch_layout[ch],
                    'radius': channels.radius[ch],
                    'id': channels.id[ch],
                    'text_names': channels.text_names[ch],
                    '3d_layout': channels.center[ch],
                    'color': c,
                }
            data['channels'] = ch_data

        with open(dest_filename, 'w') as fd:
            json.dump(data, fd, indent=4, separators=(',', ': '))
        return dest_filename

    def load(self, filename):
        with open(filename, 'r') as infile:
            self.data = json.load(infile)
        if self.ignore_ordering:
            # This is not needed e.g. when we are loading average protein which
            # as of now does not have order of SSEs refined ('helices_residues'
            # is empty)
            self.data['layout'] = OrderedDict(self.data['layout'])
        else:
            # self.data['layout'] have to be ordered as per residues
            helices_order = OrderedDict(sorted(self.data['layout'].items(), key=lambda t: t[1][0]))
            layout_raw = self.data['layout']
            self.data['layout'] = OrderedDict()
            for helic in helices_order.keys():
                self.data['layout'][helic] = layout_raw[helic]
