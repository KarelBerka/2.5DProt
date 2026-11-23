#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import requests
import logging
import time
import json
import os

class PDBeApiError(Exception):
    """Generic PDBe API error"""
    pass

class PDBeApiReturnCodeError(PDBeApiError):
    """PDBe API return code is not correct"""
    pass

class PDBeApi(object):

    API_URL = 'https://www.ebi.ac.uk/pdbe/api'
    DB = 'input/pdbe_api.json'

    def __init__(self):
        self.cache = {}
        self.loaded = False
        self.use_cache = False

    def load(self):
        if not os.path.isfile(self.DB):
            return
        with open(self.DB, 'r') as fp:
            try:
                self.cache = json.load(fp)
            except ValueError as e:
                raise ValueError("Broken %s: %s" % (self.DB, e.message))
        logging.debug("Loaded cache from %s with %s records" % (self.DB, len(self.cache)))

    def save(self):
        with open(self.DB, 'w') as fp:
            json.dump(self.cache, fp)
        logging.debug("Saved cache to %s with %s records" % (self.DB, len(self.cache)))

    def get(self, query):
        # Load cache from the disk if we have not done so yet
        if self.use_cache and not self.loaded:
            self.load()
            self.loaded = True

        if not self.use_cache and query in self.cache:
            return self.cache[query]
        else:
            url = "%s/%s" % (self.API_URL, query)

            attempt = 0
            attempt_max = 10
            while True:
                try:
                    response = requests.get(url, timeout=30)
                except requests.exceptions.Timeout as e:
                    logging.warning("Request to %s timed for a %s time: %s" % (url, attempt, e))
                    if attempt >= attempt_max:
                        logging.error("Failed the request to %s because it is timeouting" % url)
                        raise
                    attempt += 1
                    time.sleep(1)
                else:
                    break

            if response.ok:
                self.cache[query] = json.loads(response.text)
            else:
                logging.error('HTTP request %s failed, status code %s' % (url, response.status_code))
                raise PDBeApiReturnCodeError(response.status_code)
            if self.use_cache:
                self.save()
            return self.cache[query]


class PDBeApiMapping(PDBeApi):

    API_URI = 'mappings'
    CATH_OR_B = 'CATH-B'
    CATH = 'CATH'

    def get_domain_ranges(self, pdbid, chain, domain):
        """
        Return residue ranges for given domain. Note that chain we feed
        in does not necessarily equal to the chain we print out, because:

        auth_* numbering (old numbering we used before with *.pdb files)
        +-----------------------------------------------------------------------------------------------------------------------------------------------------------+
        |                 | PDBe API                                                  | mmCIF file                      | PDB file                                  |
        +-----------------+-----------------------------------------------------------+---------------------------------+-------------------------------------------+
        | Chain           | chain_id                                                  | auth_asym_id                    | chainID (column 22)                       |
        | First rezidue # | start:author_residue_number + start:author_insertion_code | auth_seq_id + pdbx_PDB_ins_code | resSeq (column 23-26) + iCode (column 27) |
        | Last rezidue #  | end:author_residue_number + end:author_insertion_code     | auth_seq_id + pdbx_PDB_ins_code | resSeq (column 23-26) + iCode (column 27) |
        +-----------------+-----------------------------------------------------------+---------------------------------+-------------------------------------------+

        label_* numbering (new numbering we use now with *.cif files)
        +-----------------+----------------------+---------------+----------+
        |                 | PDBe API             | mmCIF file    | PDB file |
        +-----------------+----------------------+---------------+----------+
        | Chain           | struct_asym_id       | label_asym_id | -        |
        | First residue # | start:residue_number | label_seq_id  | -        |
        | Last residue #  | end:residue_number   | label_seq_id  | -        |
        +-----------------+----------------------+---------------+----------+
        """
        data = self.get(self.API_URI + "/" + pdbid)
        out = {}
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if v2['domain'][0:4] == pdbid \
                   and v2['struct_asym_id'] == chain \
                   and v2['domain'][-2:] == domain:
                    out[v2['segment_id']] = (v2['start']['residue_number'], v2['end']['residue_number'])

        if out == {}:
            for k1, v1 in data[pdbid][self.CATH].items():
                for v2 in v1['mappings']:
                    if v2['domain'][0:4] == pdbid \
                       and v2['struct_asym_id'] == chain \
                       and v2['domain'][-2:] == domain:
                        out[v2['segment_id']] = (v2['start']['residue_number'], v2['end']['residue_number'])
        out = [out[i] for i in sorted(out.keys())]
        return out


    def find_domain_name(self, pdbid, chain, domain):
        data = self.get(self.API_URI + "/" + pdbid)
        out = None
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if v2['domain'][0:4] == pdbid \
                   and v2['struct_asym_id'] == chain \
                   and v2['domain'][-2:] == domain:
                    out = v2['domain']

        if out == None:
            for k1, v1 in data[pdbid][self.CATH].items():
                for v2 in v1['mappings']:
                    if v2['domain'][0:4] == pdbid \
                       and v2['struct_asym_id'] == chain \
                       and v2['domain'][-2:] == domain:
                        out = v2['domain']
        return out


    def transform_segment_id_to_chain_id(self, pdbid, chain, domain):
        data = self.get(self.API_URI + "/" + pdbid)
        out = None
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if v2['domain'][0:4] == pdbid \
                   and v2['struct_asym_id'] == chain \
                   and ((v2['domain'][-2:] == domain) or (domain == "ALL")):
                    out = v2['chain_id']

        if out == None:
            for k1, v1 in data[pdbid][self.CATH].items():
                for v2 in v1['mappings']:
                    if v2['domain'][0:4] == pdbid \
                       and v2['struct_asym_id'] == chain \
                       and ((v2['domain'][-2:] == domain) or (domain == "ALL")):
                        out = v2['chain_id']
        return out


    def transform_chain_id_to_segment_id(self, pdbid, chain):
        data = self.get(self.API_URI + "/" + pdbid)
        out = None
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if v2['domain'][0:4] == pdbid \
                   and v2['chain_id'] == chain:
                    out = v2['struct_asym_id']

        if out == None:
            for k1, v1 in data[pdbid][self.CATH].items():
                for v2 in v1['mappings']:
                    if v2['domain'][0:4] == pdbid \
                       and v2['chain_id'] == chain \
                       and v2['domain'][-2:] == domain:
                        out = v2['struct_asym_id']
        return out


    def get_domain_family(self, pdbid, chain, domain):
        """
        For given domain (e.g. "1tqnA00"), return its family (e.g. "1.10.630.10")
        """
        data = self.get(self.API_URI + "/" + pdbid)
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if v2['domain'][0:4] == pdbid \
                   and v2['chain_id'] == chain \
                   and v2['domain'][-2:] == domain:
                    return str(k1)
        domain_str = "%s%s%s" % (pdbid, chain, domain)
        raise KeyError("Fammily for domain %s not found" % domain_str)

    def get_family_domains(self, familyid):
        """
        For given family (e.g. "1.10.10.10") returns list of domains
        (e.g. ['1tqnA00', ...]
        """
        data = self.get(self.API_URI + "/" + familyid)
        out = set()
        for protein, mapping in data[familyid]['PDB'].items():
            for domain in mapping['mappings']:
                newd = str(domain['domain'])[:4]+ str(domain['struct_asym_id'])+ str(domain['domain'])[-2:]
                out.add(newd)
        return out

    def get_family_cath_domains(self, familyid):
        """
        For given family (e.g. "1.10.10.10") returns list of domains
        (e.g. ['1tqnA00', ...]
        """
        data = self.get(self.API_URI + "/" + familyid)
        out = set()
        for protein, mapping in data[familyid]['PDB'].items():
            for domain in mapping['mappings']:
                newd = str(domain['domain'])[:4]+ str(domain['chain_id'])+ str(domain['domain'])[-2:]
                out.add(newd)
        return out

    def _get_protein_chains_from(self, proteinid, out, fromkey):
        """
        For given protein (e.g. "1tqn") returns list of chains (e.g. ["A"])
        """
        data = self.get(self.API_URI + "/" + proteinid)
        for k1, v1 in data[proteinid][fromkey].items():
            for v2 in v1['mappings']:
                if v2['domain'].startswith(proteinid):
#                    chain = v2['domain'][4:][:-2]
                    chain = v2['struct_asym_id']
                    if chain not in out:
                        out.append(str(chain))
        return out

    def get_protein_chains(self, proteinid):
        """
        TODO
        """
        out = []
        self._get_protein_chains_from(proteinid, out, self.CATH_OR_B)
        self._get_protein_chains_from(proteinid, out, self.CATH)
        return list(set(out))

    def get_protein_chains_cath_only(self, proteinid):
        """
        TODO
        """
        out = []
        self._get_protein_chains_from(proteinid, out, self.CATH)
        return out

    def get_protein_chains_cath_beta_only(self, proteinid):
        """
        TODO
        """
        out = []
        self._get_protein_chains_from(proteinid, out, self.CATH_OR_B)
        return out

    def get_chain_domains(self, pdbid, chain):
        """
        For given protein+chain (e.g. "1tqnA"), return list of domains (e.g. ["00"])
        """
        data = self.get(self.API_URI + "/" + pdbid)
        out = []
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if v2['domain'].startswith(pdbid) \
                   and v2['struct_asym_id'] == chain:
                    domain = v2['domain'][-2:]
                    if domain not in out:
                        out.append(str(domain))
        for k1, v1 in data[pdbid][self.CATH].items():
            for v2 in v1['mappings']:
                if v2['domain'].startswith(pdbid) \
                   and v2['struct_asym_id'] == chain:
                    domain = v2['domain'][-2:]
                    if domain not in out:
                        out.append(str(domain))

        return out


    def get_chain_domain_family(self, pdbid, chain, domain):
        """
        For given protein+chain+domain (e.g. "1tqnA00"), return to which family it belonges (e.g. "1.10.10.630")
        """
        data = self.get(self.API_URI + "/" + pdbid)
        out = None
        for k1, v1 in data[pdbid][self.CATH_OR_B].items():
            for v2 in v1['mappings']:
                if (v2['domain'].startswith(pdbid)) \
                   and (v2['struct_asym_id'] == chain) \
                   and (v2['domain'][-2:] == domain):
                    return k1
        for k1, v1 in data[pdbid][self.CATH].items():
            for v2 in v1['mappings']:
                if (v2['domain'].startswith(pdbid)) \
                   and (v2['struct_asym_id'] == chain) \
                   and (v2['domain'][-2:] == domain):
                    return k1
        return out
