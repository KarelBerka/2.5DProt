#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path
import datetime
import atexit
import Bio.PDB
import loader
import pickle
import utils_vector
import logging
import utils_pdbe_api

class Rmsd3DCache(object):

    NAME = 'input/rmsd_3d_cache.pickle'

    def __init__(self):
        # Check cache exists, if not, create it
        if not os.path.isfile(self.NAME):
            with open(self.NAME, 'wb') as fp:
                pickle.dump({}, fp)
        # Load cache
        self.cache_changed = False
        self.cache_mtime = os.path.getmtime(self.NAME)
        with open(self.NAME, 'rb') as fp:
            self.cache = pickle.load(fp)
        # Register cleanup function
        atexit.register(self._save)

    def __contains__(self, item):
        logging.debug("DEBUG: Checking if %s is in Rmsd3DCache" % (item,))
        return item in self.cache

    def __getitem__(self, key):
        logging.debug("DEBUG: Gettinging %s from Rmsd3DCache" % (key,))
        return self.cache[key]

    def __setitem__(self, key, value):
        logging.debug("DEBUG: Settinging %s to %s in Rmsd3DCache" % (key, value))
        self.cache_changed = True
        self.cache[key] = value

    def __delitem__(self, key):
        logging.debug("DEBUG: Deleting %s from Rmsd3DCache" % (key,))
        self.cache_changed = True
        del(self.cache[key])

    def _save(self):
        if self.cache_changed:
            new_mtime = os.path.getmtime(self.NAME)
            print(">>>", self.cache_mtime, new_mtime, self.cache_mtime == new_mtime)
            if self.cache_mtime != new_mtime:
                logging.warning("WARNING: Not saving new items into Rmsd3DCache cache because timestamp changed")
            else:
                logging.debug("DEBUG: Saving Rmsd3DCache")
                with open(self.NAME, 'wb') as fp:
                    pickle.dump(self.cache, fp)


_rmsd_3d_cache = Rmsd3DCache()


def get_rmsd_3d(protein1, chain1, domain1, protein2, chain2, domain2, only_for=None):
    """
    Return RMSD of these two helices and number of matched SSEs. RMSD
    is conted as per centers of SSEs aproximations in 3D.
    """
    sses1 = set(protein1.get_helices(chain1, domain1))
    sses2 = set(protein2.get_helices(chain2, domain2))
    sses_common = list(sses1.intersection(sses2))   # we need this ordered
    if only_for is not None:
        sses_common = set(sses_common).intersection(only_for)
    assert len(sses_common) >= 3, "At least 3 common SSEs needed"

    # Count SSE centers and fake them as atoms
    points1 = []
    points2 = []
    for sse in sses_common:
        line1 = protein1.get_helic_aprox_bestline(chain1, sse)
        line2 = protein2.get_helic_aprox_bestline(chain2, sse)
        point1 = loader.AtomLike(utils_vector.center_of_line(line1))
        point2 = loader.AtomLike(utils_vector.center_of_line(line2))
        # Add points that many times as is size
        size = int(min(protein1.get_sizes(chain1, domain1)[sse], protein2.get_sizes(chain2, domain2)[sse]))
        i = 0
        while i < size:
           points1.append(point1)
           points2.append(point2)
           i = i + 1

    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(points1, points2)

    # Rotate second protein as needed
    protein2.apply_super_imposer(chain2, super_imposer)

    # Get RMSD so we do not have to count it ourselves
    rmsd = super_imposer.rms

    return (rmsd, len(sses_common))

def get_rmsd_2d(layout1, layout2, only_for=None):
    """
    Return RMSD of these two helices and number of matched SSEs. RMSD
    is conted as per centers in 2D layout.
    """
    sses1 = set(layout1.data['layout'].keys())
    sses2 = set(layout2.data['layout'].keys())
    sses_common = sses1.intersection(sses2)
    if only_for is not None:
        sses_common = set(sses_common).intersection(only_for)
    if (len(sses_common) <3):
       return (0,0)
    #assert len(sses_common) >= 3, "At least 3 common SSEs needed"

    # Fake our 2D points as 3D atoms
    points1 = []
    points2 = []
    for sse in sses_common:
        point1 = layout1.data['layout'][sse] + [0]
        point2 = layout2.data['layout'][sse] + [0]
        atom1 = loader.AtomLike(point1)
        atom2 = loader.AtomLike(point2)
        # Add points that many times as is size
        size = int(min(layout1.data['size'][sse], layout2.data['size'][sse]))
        i = 0
        while i < size:
           points1.append(atom1)
           points2.append(atom2)
           i = i + 1

    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(points1, points2)

    # Rotate points
    super_imposer.apply(points2)

    # Get RMSD so we do not have to count it ourselves
    rmsd = super_imposer.rms

    return (rmsd, len(sses_common))

def _pymot_init():
    """
    Helper function to init pymol interface
    """
    import pymol
    pymol.finish_launching(['pymol', '-cq'])
    return pymol

def _pymol_load(pymol, protein, chain, domain):
    """
    Helper function to load protein into Pymol if it was not loaded
    before and create variable per chain&domain if it was not created already
    """
    #print("DEBUG: Loading %s:%s" % (protein, chain))
    pymol.cmd.load(protein.filename_protein)
    assert chain in pymol.cmd.get_chains(protein.name)

    domain_borders = utils_pdbe_api.get_domain_ranges(protein.name, chain, domain)
    domain_borders_str = ','.join(['%s:%s' % borders for borders in domain_borders])
    assert len(domain_borders) > 0

    fraction = "%s_%s%s" % (protein.name, chain, domain)
    pymol.cmd.create(fraction, "%s and chain %s and resi %s" \
        % (protein.name, chain, domain_borders_str))
    return fraction

def get_rmsd_3d_pymol(pymol,fraction1, fraction2):
    # If we have the value in the cache, return it
    cache1_key = (fraction1, fraction2)
    cache2_key = (fraction2, fraction1)
    if cache1_key in _rmsd_3d_cache:
        return _rmsd_3d_cache[cache1_key][1:3]
    if cache2_key in _rmsd_3d_cache:
        return _rmsd_3d_cache[cache2_key][1:3]

    try:
        cealign_out = pymol.cmd.cealign(fraction1, fraction2, quiet=0)
    except pymol.CmdException as e:
        logging.error("ERROR: cealign %s %s failed with: %s" % (fraction1, fraction2, e))
        return 'Err'
    logging.debug("DEBUG: Pymol RMSD: %s vs. %s: %s" % (fraction1, fraction2, cealign_out))

    # Add newly couned value into cache
    _rmsd_3d_cache[cache1_key] = (datetime.datetime.now(), cealign_out['RMSD'], cealign_out['alignment_length'])

    return (cealign_out['RMSD'], cealign_out['alignment_length'])
