#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path
import csv
import logging
import math
import numpy
import json
from collections import OrderedDict
try:
    import Bio.PDB
    _AtomBase = Bio.PDB.Atom.Atom
except ImportError:
    Bio = None
    _AtomBase = object
import utils
import utils_image
try:
    import two_d_prot.pdbe_api
    import two_d_prot.pdbe_downloader
    import two_d_prot.parser_workaround
except ImportError:
    two_d_prot = None


# used in utils_rmsd.py
class AtomLike(_AtomBase):

    def __init__(self, coord, name=None):
        self.coord = coord
        if name is not None:
            self.id = name


# used in Family and Domain classes
class Item(object):

    pdbe_api_mapping = two_d_prot.pdbe_api.PDBeApiMapping() if two_d_prot else None


class Protein(object):

    QUIET = True   # do not print warnings on PDB structure when loading
    DSSP_PATH = os.path.join(os.getcwd(), 'dssp')

    def __init__(self, filename, protein_name=False):
        if filename.endswith('.pdb') or filename.endswith('.cif'):
            filename = filename[:-4]
        self.domain_name = os.path.basename(filename)
        if protein_name:
            self.name = protein_name
        else:
            self.name = self.domain_name[:4]
        self.filename_protein = filename+".cif"
        assert os.path.exists(self.filename_protein)
        self.domains = {}   # loaded later from CathDB
        self._load_structure()
        self.heme = None   # to be loaded upon first call to self.get_heme()
        self.m_distances = {}
        self.distances = {}
        self.sizes = {}
        self._pdbe_api_mapping_object = None   # use self._pdbe_api_mapping to access this to allow lazy loading
        self.secstr = {}
        self.secstr_atoms = {}
        self.beta_connectivity_next = {}
        self.beta_connectivity_prev = {}
        self.whole_protein = 0
        self.sizes_in_res = {}

    def _load_structure(self):
        parser = two_d_prot.parser_workaround.MMCIFParser(QUIET=self.QUIET)
        self.structure = parser.get_structure(self.name, self.filename_protein)
        self.header = self.structure.header
        if len(list(self.structure.get_models())) > 1:
            print("WARNING: There are more models in %s file, using only first one" % self.filename_protein)

    #
    def add_chain_domain(self, chain, domain, filename=None):
        """
        Any method that accesses any data related to annotations should call
        this method at the beginning. It will make sure that we do not have
        outdated data in various caches (scleanup is in the `_load_secstrs()`)
        and that we have correct annotations loaded.

        E.g. protein "3s99" have A01 and A02 both in one family. If we would
        use anontation created for A02 when accessing A01, we would get
        incorrect data (and vice versa).
        """
        if filename is None or filename == "None":
            filename_anotation = "input/annotations/%s-%s%s.json" \
                                 % (self.name, chain, domain)
        else:
            filename_anotation = filename

        if os.path.isfile(filename_anotation) is False:
            filename_anotation = "input/annotations/%s.json" % (self.name)
        self._load_secstrs(filename_anotation, chain, domain)

    def add_chain_domain_from_file(self, chain, domain, filename, test=True):
        """    """
        self._load_secstrs(filename, chain, domain, test)

    def _load_secstrs(self, filename_anotation, chain, domain, test=True):
        # add secondary structures for chain, domain to directory, use annotation file
        assert os.path.isfile(filename_anotation)
        with open(filename_anotation, 'r') as fd:
            data = json.load(fd)
        if chain not in self.secstr:
            self.secstr[chain] = {}
            self.secstr_atoms[chain] = {}
            self.beta_connectivity_next[chain] = {}
            self.beta_connectivity_prev[chain] = {}
        if domain not in self.secstr[chain]:
            self.secstr[chain][domain] = {}
            self.secstr_atoms[chain][domain] = {}
            self.beta_connectivity_next[chain][domain] = {}
            self.beta_connectivity_prev[chain][domain] = {}
        if self.name.lower() in data:
            name = self.name.lower()
        elif self.name in data:
            name = self.name
        elif self.domain_name in data:
            name = self.domain_name
        for item in data[name]['secondary_structure_elements']:
            # there is the possibility to have concreate chain or metacharacter ALL
            if item['chain_id'] == chain or chain == "ALL":
                if not test or \
                    (self.is_res_in_domain(chain, domain, (' ', item['start'], ' '))
                        or self.is_res_in_domain(chain, domain, (' ', item['end'], ' '))
                        or domain == "ALL"):
                        self.secstr[chain][domain][item['label']] = (item['start'], item['end'])
        # we have to load beta connectivity data as well
        for item in data[name]['beta_connectivity']:
            if (chain in self.secstr) and (domain in self.secstr[chain]) and (item[0] in self.secstr[chain][domain]) and (item[1] in self.secstr[chain][domain]):
                if item[0] not in self.beta_connectivity_next[chain][domain]:
                    self.beta_connectivity_next[chain][domain][item[0]] = []
                self.beta_connectivity_next[chain][domain][item[0]].append([item[1], item[2]])
                if item[1] not in self.beta_connectivity_next[chain][domain]:
                    self.beta_connectivity_next[chain][domain][item[1]] = []
                self.beta_connectivity_next[chain][domain][item[1]].append([item[0], item[2]])

    #  FIXME: Maybe we should have some checks per superfamily or something? This works for CYP, but not for GPCR
    # for chain in self.secstr:
    #    secstr_prim = [k for k,v in self.secstr[chain].iteritems() if v[0] == self.HELIC_TYPE_PRIM]
    #    assert set(secstr_prim) == set(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']), \
    #        "Annotation for protein %s, chain %s defines unexpected set of primary helices: %s" % (self.name, chain, ",".join(secstr_prim))
    # Make sure it is sorted by residue ID
    #    for chain in self.secstr[chain][domain]:
    #        self.secstr[chain][domain] = OrderedDict(sorted(self.secstr[chain][domain].items(), key=lambda t: t[1][0]))

    def add_chain_domain_from_layout(self, chain, domain, layout3d):
        self.whole_protein = 1
        if chain not in self.secstr:
            self.secstr[chain] = {}
        if domain not in self.secstr[chain]:
            self.secstr[chain][domain] = {}
        # read secstr data
        for sses in layout3d:
            self.secstr[chain][domain][sses] = (layout3d[sses])
    # we have to load beta connectivity data as well
#    for item in data[self.name.lower()]['beta_connectivity']:
#        if item[0] not in self.beta_connectivity_next[chain][domain]:
#            self.beta_connectivity_next[chain][domain][item[0]]=[]
#        self.beta_connectivity_next[chain][domain][item[0]].append([item[1],item[2]])
#        if item[1] not in self.beta_connectivity_next[chain][domain]:
#            self.beta_connectivity_next[chain][domain][item[1]]=[]
#        self.beta_connectivity_next[chain][domain][item[1]].append([item[0],item[2]])

    def get_beta_connectivity_next(self, chain, domain, sheet):
        """
        return the sheets which are next in beta connections section
        """
        if sheet not in self.beta_connectivity_next[chain][domain].keys():
            return []
        else:
            return self.beta_connectivity_next[chain][domain][sheet]

    def get_beta_connectivity_prev(self, chain, domain, sheet):
        """
        return the sheets which are previous in beta connections section
        """
        if sheet not in self.beta_connectivity_prev[chain][domain].keys():
            return []
        else:
            return self.beta_connectivity_prev[chain][domain][sheet]

    def get_beta_connectivity_near(self, chain, domain, sheet):
        """
        return the sheets which are near (next or previous) in
        beta connections section
        """
        if sheet not in self.beta_connectivity_next[chain][domain].keys():
            if sheet not in self.beta_connectivity_prev[chain][domain].keys():
                near = []
            else:
                near = self.beta_connectivity_prev[chain][domain][sheet]
        else:
            if sheet not in self.beta_connectivity_prev[chain][domain].keys():
                near = self.beta_connectivity_next[chain][domain][sheet]
            else:
                near = self.beta_connectivity_next[chain][domain][sheet] + self.beta_connectivity_prev[sheet]
        return near

    def get_beta_connectivity_cluster_orientation(self, chain, domain, start):
        """
        return cluster with decomposition to same/opposit orientation sets
        """
        # set with sheets with the same/opposit orientation
        sameset = [start]
        opposit = []

        # wariables which helps to track which sses have to be
        to_do = [start]                                             # parsed
        cluster = []                                              # are already parsed

        while to_do != []:
            sheet = to_do[0]
            nextsheet = self.get_beta_connectivity_next(chain, domain, sheet)
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

    def get_beta_connectivity_clusters(self, chain, domain):
        """
        return the list of sheets which have the opposit orientation
        """
        output = []
        for h in self.get_helices(chain, domain):
            isnot = 0
            for i in output:
                if h in i:
                    isnot = 1
            # if h is not in any added cluster, ad it cluster now
            set = self.get_beta_connectivity_cluster_orientation(chain, domain, h)
            if isnot == 0:
                output.append(":".join(set[0]))
        return output

    def _get_residue_atoms(self, chain, residue):
        return list(self.structure[0][chain][residue])

    def _get_secstr_atoms(self, chain, domain, secstr):
        # If we already have it, just return it
        if secstr not in self.secstr_atoms[chain][domain]:
            # Load atoms for given secstr
            self.secstr_atoms[chain][domain][secstr] = []
            if secstr in self.secstr[chain][domain]:
                res_range = range(self.secstr[chain][domain][secstr][0], self.secstr[chain][domain][secstr][1]+1)
                for r in self.structure[0][chain].get_residues():
                    if r.get_id()[1] in res_range:
                        for a in r:
                            self.secstr_atoms[chain][domain][secstr].append(a)
            else:
                raise Exception("Atoms for unknown helic %s, chain %s requested " % (secstr, chain))
        return self.secstr_atoms[chain][domain][secstr]

    @property
    def _pdbe_api_mapping(self):
        if self._pdbe_api_mapping_object is None:
            self._pdbe_api_mapping_object = two_d_prot.pdbe_api.PDBeApiMapping()
        return self._pdbe_api_mapping_object

    def apply_super_imposer(self, chain, super_imposer):
        super_imposer.apply(self.structure[0][chain].get_atoms())

    def info(self):
        out = ''
        out += "Name: %s\n" % self.name
        out += "Files: %s, %s\n" % (self.filename_protein, self.filename_anotation)
        out += "Structure: %s\n" % self.structure
        out += "First model: %s\n" % self.structure[0]
        out += "Chains: %s\n" % (','.join(self.structure[0].get_chains()))
        first = sorted(self.secstr.keys())[0]
        out += "Secondary structures in chain %s: %s\n" % (self.secstr[first], first)
        out += "Heme in chain %s: %s (%s) -> %s" % (first, self.get_heme(first).get_resname(), self.get_heme_center(first), str(self.get_heme(first).get_full_id()))
        return out

    def dump(self):
        # Heme
        csvfile = 'data/%s/heme.dat' % self.name
        print("Writing heme coords into %s" % csvfile)
        csvout = open(csvfile, 'wb')
        csvwriter = csv.writer(csvout, delimiter='\t')
        csvwriter.writerow(self.get_heme_center())
        csvout.close()
        # Atoms per helices
        for h in self.get_helices():
            csvfile = 'data/%s/atoms-%s.dat' % (self.name, h)
            print("Writing helic %s atom coords into %s" % (h, csvfile))
            csvout = open(csvfile, 'wb')
            csvwriter = csv.writer(csvout, delimiter='\t')
            for a in self.get_secstr_atoms(h):
                csvwriter.writerow(a.get_coord())
            csvout.close()

    def get_dssp(self, chain, domain):
        if self.dssp == {}:
            self.dssp = Bio.PDB.DSSP(self.structure[0], self.filename_protein, dssp=self.DSSP_PATH)
        filtered = {k[1]: v for k, v in dict(self.dssp).iteritems() if k[0] == chain}   # k is e.g. ('A', (' ', 106, ' ')), filter only these with our chain
        filtered = {k: v for k, v in filtered.iteritems() if self.is_res_in_domain(chain, domain, k)}   # k is e.g. (' ', 106, ' '), filter only what belonfs to our domain
        ordered = OrderedDict(sorted(filtered.items(), key=lambda t: t[1][0]))   # order by residuum ID
        return ordered

    #
    def get_distances(self, chain, domain):
        if chain not in self.distances:
            self.distances[chain] = {}
        if domain not in self.distances[chain]:
            self.distances[chain][domain] = utils.get_distance_matrix4(self, chain, domain)
        return self.distances[chain][domain]

    def get_m_distances(self, chain1, domain1, chain2, domain2, short=False):
        if chain1 not in self.m_distances:
            self.m_distances[chain1] = {}
        if domain1 not in self.m_distances[chain1]:
            self.m_distances[chain1][domain1] = {}
        if chain2 not in self.m_distances[chain1][domain1]:
            self.m_distances[chain1][domain1][chain2] = {}
        if domain2 not in self.m_distances[chain1][domain1][chain2]:
            pp = utils.get_distance_matrix4_2(self, chain1, domain1, chain2, domain2, short=short)
            self.m_distances[chain1][domain1][chain2][domain2] = pp
        return self.m_distances[chain1][domain1][chain2][domain2]

    def get_cluster_distances(self, chain, domain):
        sse_dist = self.get_distances(chain, domain)
        if chain not in self.cl_distances:
            self.cl_distances[chain] = {}
        if domain not in self.cl_distances[chain]:
            self.cl_distances[chain][domain] = {}
            for i in self.get_beta_connectivity_clusters(chain, domain):
                self.cl_distances[chain][domain][i] = {}
                for j in self.get_beta_connectivity_clusters(chain, domain):
                    ii = i.split(":")
                    jj = j.split(":")
                    dist_sum = 0
                    dist_number = 0
                    for h1 in ii:
                        for h2 in jj:
                            dist_sum = dist_sum + sse_dist[h1][h2]
                            dist_number = dist_number + 1
                    self.cl_distances[chain][domain][i][j] = float(dist_sum)/dist_number
        return self.cl_distances[chain][domain]

    def get_sizes(self, chain, domain):
        if chain not in self.sizes:
            self.sizes[chain] = {}
        if domain not in self.sizes[chain]:
            self.sizes[chain][domain] = {}
            for helic in self.get_helices(chain, domain):
                if self.whole_protein == 1:
                    bestline = self.secstr[chain][domain][helic]
                    self.sizes[chain][domain][helic] = math.dist(numpy.array(bestline[0]), numpy.array(bestline[1]))
                else:
                    bestline = self.get_helic_aprox_bestline(chain, domain, helic)
                    self.sizes[chain][domain][helic] = math.dist(bestline[0], bestline[1])
        return self.sizes[chain][domain]

    def get_cluster_sizes(self, chain, domain):
        if chain not in self.cl_sizes:
            self.cl_sizes[chain] = {}
        if domain not in self.cl_sizes[chain]:
            self.cl_sizes[chain][domain] = {}
            for cluster in self.get_beta_connectivity_clusters(chain, domain):
                self.cl_sizes[chain][domain][cluster] = 0
                clusters = cluster.split(":")
                for sse in clusters:
                    bestline = self.get_helic_aprox_bestline(chain, domain, sse)
                    self.cl_sizes[chain][domain][cluster] = self.cl_sizes[chain][domain][cluster] + math.dist(bestline[0], bestline[1])
        return self.cl_sizes[chain][domain]

    #
    def get_sizes_in_res(self, chain, domain):
        if chain not in self.sizes_in_res:
            self.sizes_in_res[chain] = {}
        if domain not in self.sizes_in_res[chain]:
            self.sizes_in_res[chain][domain] = {}
        for helic in self.get_helices(chain, domain):
            self.sizes_in_res[chain][domain][helic] = self.secstr[chain][domain][helic][1] - self.secstr[chain][domain][helic][0]
        return self.sizes_in_res[chain][domain]

    #
    def get_all_helices_pairs_angles(self, chain, domain):
        # get angles between all helic pairs
        helices = self.get_helices(chain, domain)

        # at first approximate helices with lines
        lines = {}
        for h in helices:
            l = self.get_helic_aprox_bestline(chain, domain, h)
            lines[h] = [l[0][0]-l[1][0], l[0][1]-l[1][1], l[0][2]-l[1][2]]

        # count angles between approximated lines
        angle = {}
        for h in helices:
            for g in helices:
                #  this is angle between h and g in degrees values 0..179
                if h != g:
                    angle[h, g] = round(math.acos(
                      (lines[h][0] * lines[g][0] + lines[h][1]*lines[g][1] + lines[h][2] * lines[g][2]) /
                      (utils.vector_lenght(lines[h]) * utils.vector_lenght(lines[g]))) * 180 / math.pi)
                else:
                    angle[h, g] = 0
        return angle

    def get_chains(self):
        return sorted([k.id for k in self.structure[0].get_chains()])

    def get_domains(self, chain):
        if chain not in self.domains:
            self.domains[chain] = {}
            for domain in self._pdbe_api_mapping.get_chain_domains(self.name, chain):
                self.domains[chain][domain] = self._pdbe_api_mapping.get_domain_ranges(self.name, chain, domain)
        return self.domains[chain]

    def is_res_in_domain(self, chain, domain, res):
        # res ... something like "(' ', 2, ' ')"
        borders = self._pdbe_api_mapping.get_domain_ranges(self.name, chain, domain)
        for border_pair in borders:
            if border_pair[0] <= res[1] <= border_pair[1]:
                return True
        return False

    #
    def get_helices(self, chain, domain):
        out = []
        if (chain not in self.secstr):
            return out
        elif (domain not in self.secstr[chain]):
            return out
        if self.whole_protein == 1:
            return self.secstr[chain][domain]
        elif (domain == "ALL") or (domain == "re"):
            return self.secstr[chain][domain]
        else:
            for sse in self.secstr[chain][domain].keys():
                # FIXME: We might want to check if whole range is in domain
                if self.is_res_in_domain(chain, domain, (' ', self.secstr[chain][domain][sse][0], ' ')) \
                    or self.is_res_in_domain(chain, domain, (' ', self.secstr[chain][domain][sse][1], ' ')):
                    out.append(sse)
        return out

    def get_secstr_residues(self, chain, domain=None, secstr=None):
        if secstr is None:
            assert domain is not None
            out = OrderedDict()
            if (chain not in self.secstr):
                return out
            for k, v in self.secstr[chain][domain].items():
                if (domain == "ALL") or \
                    (self.is_res_in_domain(chain, domain, (' ', v[0], ' ')) or self.is_res_in_domain(chain, domain, (' ', v[1], ' '))):
                    out[k] = (v[0], v[1])
            return out
        else:
            return self.secstr[chain][secstr]

    def get_residue_atoms(self, chain, residue):
        out = self._get_residue_atoms(chain, residue)
        assert len(out) > 0, "No atoms in residue %s chain %s? Maybe this is error in %s.sses file?" % (residue, chain, self.name)
        return out

    def get_secstr_atoms(self, chain, domain, secstr):
        out = self._get_secstr_atoms(chain, domain, secstr)
        assert len(out) > 0, "No atoms in helic %s chain %s domain %s? Maybe this is error in %s.sses file?" % (secstr, chain, domain, self.name)
        return out

    def get_residue_atoms_coords(self, chain, residue):
        return [a.get_coord() for a in self.get_residue_atoms(chain, residue)]

    def get_secstr_atoms_coords(self, chain, domain, helic):
        return [a.get_coord() for a in self.get_secstr_atoms(chain, domain, helic)]

    def _get_coords_aprox_bestline(self, coords):
        # data = [(x1, y1, z1), (x2, y2, z2),...]
        # Thaks http://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
        data = numpy.array(coords)
        data_mean = data.mean(axis=0)
        data_normalized = data - data_mean
        # Using 3 * standard deviation, we should catch 99.7% of all points (for normal distribution)
        # https://en.wikipedia.org/wiki/Standard_deviation#Rules_for_normally_distributed_data
        std_dev_3 = data_normalized.std() * 3
        uu, dd, vv = numpy.linalg.svd(data_normalized)
        linepts = vv[0] * numpy.mgrid[-std_dev_3:std_dev_3:2j][:, numpy.newaxis]
        linepts += data_mean
        # Verify that line have same orientation (i.e. start <---> end matches) as underlying atoms
        oriented = 0
        disoriented = 0
        for i in range(3):
            oriented += math.dist(linepts[0], data[i])
            disoriented += math.dist(linepts[0], data[-1*(i+1)])
            oriented += math.dist(linepts[1], data[-1*(i+1)])
            disoriented += math.dist(linepts[1], data[i])
        if oriented < disoriented:
            return linepts
        else:
            return [linepts[1], linepts[0]]

    def get_residues_aprox_bestline(self, chain, residues):
        coords = []
        for r in residues:
            try:
                # Going to ignore KeyError here, because e.g. for "3tmm" residues
                # 85, 92 and 101 do not have "ATOM" entries, only "HETATM" entries
                # and that causes KeyError here:
                #
                #   [...]
                #   File "/home/pok/Checkouts/bilkoviny/loader.py", line 118, in _get_residue_atoms
                #     return list(self.structure[0][chain][residue])
                #   File "/usr/lib64/python2.7/site-packages/Bio/PDB/Chain.py", line 74, in __getitem__
                #     return Entity.__getitem__(self, id)
                #   File "/usr/lib64/python2.7/site-packages/Bio/PDB/Entity.py", line 40, in __getitem__
                #     return self.child_dict[id]
                # KeyError: (' ', 85, ' ')
                #
                # So that KeyError sounds reasonable. Sent email asking this.
                coords += self.get_residue_atoms_coords(chain, r)
            except KeyError as e:
                print("WARNING: Failed to get atoms for res %s, ignoring. Error was: %s" % (r, e))
        return self._get_coords_aprox_bestline(coords)

    def get_helic_read_line22(self, chain, domain, helic):
        return self.secstr[chain][domain][helic]

    def get_helic_aprox_bestline(self, chain, domain, helic):
        coords = self.get_secstr_atoms_coords(chain, domain, helic)
        return self._get_coords_aprox_bestline(coords)

    def get_helic_aprox_startend_all(self, chain, domain, set):
        se = {}
        for i in set:
            se[i] = [self.get_secstr_atoms(chain, domain, i)[0].get_coord(), self.get_secstr_atoms(chain, domain, i)[-1].get_coord()]
        return se

    def get_helic_aprox_startend(self, chain, domain, helic):
        return [self.get_secstr_atoms(chain, domain, helic)[0].get_coord(), self.get_secstr_atoms(chain, domain, helic)[-1].get_coord()]


# used in utils_generate_svg.py
class Layout(object):

    COUNT_LAYOUT_MAX_GEN = 20000

    def __init__(self, filename=None, sses_stats=None, protein=None, ignore_ordering=False):
        # If we are loading already generated layout, filename have to be set
        self.ignore_ordering = ignore_ordering
        self.filename = filename
        if self.filename is not None and protein is None:
            self.load(self.filename)
            self.protein = None
            self.sses_stats = None
        elif self.filename is None and protein is not None:
            self.data = {}
            self.protein = protein
            self.sses_stats = sses_stats
        else:
            raise Exception("Unknown set of options when creating layout")

    def get_canvas_size(self, sses_filter, avg=1, min_canvas_width=None, filename_out=None):
        # Find how far to the up/right/down/left does the protein go and determine
        # dicmentions and center of the image
        xmin, xmax = "p", "p"
        ymin, ymax = "p", "p"
        for helic, coord in self.data['layout'].items():
            if (sses_filter(helic, self.data['h_type'], self.data['size']) == 0):
                if (xmax == "p") or (coord[0] > xmax):
                    xmax = coord[0]
                if (xmin == "p") or (coord[0] < xmin):
                    xmin = coord[0]
                if (ymax == "p") or (coord[1] > ymax):
                    ymax = coord[1]
                if (ymin == "p") or (coord[1] < ymin):
                    ymin = coord[1]

        if 'ligands' in self.data:
            for l in self.data['ligands']:
                if (xmax == "p") or (self.data['ligands'][l]['layout'][0] > xmax):
                    xmax = self.data['ligands'][l]['layout'][0]
                if (xmin == "p") or (self.data['ligands'][l]['layout'][0] < xmin):
                    xmin = self.data['ligands'][l]['layout'][0]
                if (ymax == "p") or (self.data['ligands'][l]['layout'][1] > ymax):
                    ymax = self.data['ligands'][l]['layout'][1]
                if (ymin == "p") or (self.data['ligands'][l]['layout'][1] < ymin):
                    ymin = self.data['ligands'][l]['layout'][1]

        if xmin == "p":
            xmin, xmax, ymin, ymax = 0, 0, 0, 0

        width = xmax - xmin
        height = ymax - ymin
        dx = xmin
        dy = ymin

        # Find biggest helic and add its radius to determined sizes
        r_max = 0
        for h, r in self.data['size'].items():
            if r > r_max:
                r_max = r
        width += 2 * r_max * utils_image.Image.CIRCLE_MULTIPLIER
        height += 2 * r_max * utils_image.Image.CIRCLE_MULTIPLIER
        dx -= r_max * utils_image.Image.CIRCLE_MULTIPLIER
        dy -= r_max * utils_image.Image.CIRCLE_MULTIPLIER

        # Add 10% of space from each side to frame the image nicely
        dx = dx - width * 0.05
        width = width * 1.1
        dy = dy - height * 0.005 - 5
        height = height * 1.01 + 5

        return (width, height), (dx, dy)

    def get_canvas(self, dimensions=None, min_dimensions=None, filename_out=None, fill=True):
        # Determine requested canvas dimensions
        if dimensions is None:
            width_height, center = self.get_canvas_size(min_canvas_width=min_dimensions[0])
            # Add extra 10 space on the upper end for all the text we are going to put in
            width_height = (width_height[0], width_height[1]+10)
        else:
            width_height, center = dimensions
            if min_dimensions is not None:
                if width_height[0] < min_dimensions[0]:
                    width_height = (min_dimensions[0], width_height[1])
                if width_height[1] < min_dimensions[1]:
                    width_height = (width_height[0], min_dimensions[1])

        # Instantiate image object
        if filename_out is None:
            filename_out = "generated/imageCircle-%s-%s-%s%s.svg" % (self.data['protein_group_safe'], self.data['protein'], self.data['chain_id'], self.data['domain'])
        return utils_image.Image(filename_out, width_height, center=center, fill=fill)

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
#            assert len(self.data['helices_residues']) > 0
            helices_order = OrderedDict(sorted(self.data['layout'].items(), key=lambda t: t[1][0]))
            layout_raw = self.data['layout']
            self.data['layout'] = OrderedDict()
            for helic in helices_order.keys():
                self.data['layout'][helic] = layout_raw[helic]
