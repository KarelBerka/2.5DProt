#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import utils_common
import copy


PDBE_URL = "https://www.ebi.ac.uk/pdbe"
MAPPINGS_URL = PDBE_URL + "/api/mappings"

CATH_B = 'CATH-B'
CATH = 'CATH'

# --------------------
# per family functions
# --------------------


def get_family_domains(familyid):
    """
    For given family (e.g. "1.10.10.10") returns list of domains
    (e.g {('1yfl', 'E', '02'), ('4gbe', 'D', '02'), ('1q0t', 'A', '02'), }
    """
    url = f"{MAPPINGS_URL}/{familyid}"
    data = utils_common.get_url_to_json(url)

    out = set()
    for protein, mapping in data[familyid]['PDB'].items():
        for domain in mapping['mappings']:
            # e.g. for 5w97h00 - (5w97, h, 00)
            newd = (domain['domain'][:4],
                    domain['domain'][4:-2],
                    domain['domain'][-2:])
            out.add(newd)
    return out


# ---------------------
# per protein functions
# ---------------------

def protein_uniprot_mapping(proteinid):
    """
    For given protein (e.g. "1m1n") returns list of uniprots
    (e.g ['P068782', 'P783946']
    """
    url = f"{MAPPINGS_URL}/{proteinid}"
    data = utils_common.get_url_to_json(url)

    uniprot_list = set()
    for uniprot, mapping in data[proteinid]["UniProt"].items():
        uniprot_list.add(uniprot)
    return uniprot_list



def get_protein_chains(proteinid):
    """
    For given protein (e.g. "1m1n") returns list of chains
    (e.g ('A', 'B', 'C') )
    """
    url = f"{MAPPINGS_URL}/{proteinid}"
    data = utils_common.get_url_to_json(url)
    chain_list = set()
    for family, mapping in data[proteinid][CATH_B].items():
        for domain in mapping['mappings']:
            if domain['domain'].startswith(proteinid):
                chain_list.add(domain['struct_asym_id'])
    return chain_list


def get_protein_domain_families(proteinid):
    """
    For given protein (e.g. "1m1n") returns list of domains
    (e.g ['2.140.10.30', '2.140.10.30', '3.40.50.1820', '3.40.50.1820']
    """
    url = f"{PDBE_URL}/{proteinid}"
    data = utils_common.get_url_to_json(url)

    out = []
    domain_list = set()
    for family, mapping in data[proteinid][CATH_B].items():
        for domain in mapping['mappings']:
            if domain['domain'] not in domain_list:
                out.append(family)
                domain_list.add(domain['domain'])
    return out


def get_protein_chain_domains_and_trans(proteinid, sses):
    """
    For given protein (e.g. "3f0e")
    and list of sses per chain (e.g. 
    {'A': [[23, 37], [41, 43], [46, 48], [52, 54], [62, 75], [80, 84], [96, 110], [113, 122], [128, 143], [147, 150], [150, 155], [162, 165], [168, 179]], \
    'B': [[24, 36], [41, 43], [46, 48], [53, 54], [62, 75], [80, 83], [96, 110], [113, 122], [128, 143], [147, 150], [153, 155], [162, 165], [168, 178]], \
    'C': [[24, 37], [41, 43], [46, 48], [52, 54], [62, 75], [80, 83], [96, 110], [113, 122], [128, 143], [150, 155], [162, 165], [168, 178]]} )
    returns list of families per domain ( and "Rest" if there is a free part)
    (e.g {'A': ['3.30.1330.50'], 'B': ['3.30.1330.50'], 'C': ['3.30.1330.50']}
    and transition table for mapping the domain names
    (e.g. {'A': ['A'], 'B': ['B'], 'C': ['C']})
    """

    url = f"{MAPPINGS_URL}/{proteinid}"
    data = utils_common.get_url_to_json(url)

    out = {}
    trans = {}
    for c in sses:
        out[c] = []
        trans[c] = []
    sses_output = copy.deepcopy(sses)
    for family, mapping in data[proteinid][CATH_B].items():
        for domain in mapping['mappings']:
            out[domain["struct_asym_id"]].append(family)
            trans[domain["struct_asym_id"]].append(domain["chain_id"])
            for sse in sses[domain["struct_asym_id"]]:
                if (((domain["start"]["residue_number"] < sse[0]) and (sse[0] < domain["end"]["residue_number"])) or
                    ((domain["start"]["residue_number"] < sse[1]) and (sse[1] < domain["end"]["residue_number"]))) \
                    and (sse in sses_output[domain["struct_asym_id"]]):
                        sses_output[domain["struct_asym_id"]].remove(sse)
    sses = copy.deepcopy(sses_output)
    for d in sses:
        if len(sses[d]) > 1:
            out[d].append("Rest")
    return out, trans


def get_protein_domains(proteinid):
    """
    For given protein (e.g. "1m1n")
    (e.g {"1m1nD04", "1m1nF04"} }
    """

    url = f"{MAPPINGS_URL}/{proteinid}"
    data = utils_common.get_url_to_json(url)

    out = set()
    for family, mapping in data[proteinid][CATH_B].items():
        for domain in mapping['mappings']:
            out.add(domain["domain"])
    return out


def download_protein_cif(directory, protein):
    pdbe_protein_url = f"{PDBE_URL}/entry-files/download/{protein}.cif"
    dest_file = f"{directory}/{protein}.cif"
    ret = utils_common.get_url_to_file(pdbe_protein_url, dest_file)
    return ret


def download_protein_pdb(directory, protein):
    pdbe_protein_url = f"{PDBE_URL}/entry-files/download/{protein}.pdb"
    dest_file = f"{directory}/{protein}.pdb"
    ret = utils_common.get_url_to_file(pdbe_protein_url, dest_file)
    return ret



def download_protein_list_cif(directory, protein_list):
    for protein in protein_list:
        out_file = "%s/%s.cif" \
            % (directory, protein)
        if not os.path.isfile(out_file):
            download_protein_cif(directory, protein)


def download_protein_list_pdb(directory, protein_list):
    for protein in protein_list:
        out_file = "%s/%s.pdb" \
            % (directory, protein)
        if not os.path.isfile(out_file):
            download_protein_pdb(directory, protein)


def transform_uniprot_chain_id_to_struct_asym_id(protein_id, uniprot_id, chain_id):
    """
    For given proteinid and chainid found for uniprot protein
    e.g. in https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/
    transform chain_id to struct_asym_id
    using https://www.ebi.ac.uk/pdbe/api/mappings/7pzy
    """
    url = f"{MAPPINGS_URL}/{protein_id}"
    data = utils_common.get_url_to_json(url)

    for domain in data[protein_id]['UniProt'][uniprot_id]["mappings"]:
        if domain['chain_id'] == chain_id:
            return domain["struct_asym_id"]
    return ""


def transform_cath_chain_id_to_struct_asym_id(protein_id, chain_id):
    """
    For given proteinid and chainid found for uniprot protein
    e.g. in https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/
    transform chain_id to struct_asym_id
    using https://www.ebi.ac.uk/pdbe/api/mappings/7pzy
    """
    url = f"{MAPPINGS_URL}/{protein_id}"
    data = utils_common.get_url_to_json(url)
    for fam_dom in data[protein_id][CATH_B]:
        for domain in data[protein_id][CATH_B][fam_dom]['mappings']:
            if domain['chain_id'] == chain_id:
                return domain["struct_asym_id"]
    return ""


def transform_cath_struct_asym_id_to_chain_id(protein_id, chain_id):
    """
    For given proteinid and chainid found for uniprot protein
    e.g. in https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/
    transform chain_id to struct_asym_id
    using https://www.ebi.ac.uk/pdbe/api/mappings/7pzy
    """
    url = f"{MAPPINGS_URL}/{protein_id}"
    data = utils_common.get_url_to_json(url)
    for fam_dom in data[protein_id][CATH_B]:
        for domain in data[protein_id][CATH_B][fam_dom]['mappings']:
            if domain['struct_asym_id'] == chain_id:
                return domain["chain_id"]
    return ""

# --------------------
# per chain functions
# --------------------


def get_protein_chain_domains(proteinid, chain_id):
    """
    For the given protein and chain this function returns
    the set of its chains
    E.G. "1m1n", "A" returns {"01", "02"}
    """

    url = f"{MAPPINGS_URL}/{proteinid}"
    data = utils_common.get_url_to_json(url)

    dom = set()
    for family, mapping in data[proteinid][CATH_B].items():
        for domain in mapping['mappings']:
            if domain["struct_asym_id"] == chain_id:
                dom.add(domain['domain'][-2:])

    return dom

# --------------------
# per domain functions
# --------------------


def get_domain_ranges(pdbid, chain, domain, family = ""):
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
    if domain is "" read all domains
    if family is set use only domains from this family
    """
    data = utils_common.get_url_to_json(MAPPINGS_URL + '/' + pdbid)
    out = {}
    for k1, v1 in data[pdbid][CATH_B].items():
        if (family == "") or (family == k1):
            # return only domains from wanted family,
            # if there is no wanted family return alldomains
            for v2 in v1['mappings']:
                if domain == "":
                    # we want all domains in the chain
                    if v2["chain_id"] == chain:
                        out[v2['segment_id']] = (v2['start']['residue_number'], v2['end']['residue_number'])
                else:
                    # we want only 'domain' domain
                    if v2['domain'] == pdbid + chain + domain:
                        out[v2['segment_id']] = (v2['start']['residue_number'], v2['end']['residue_number'])
    if out == {}:
        # if there was no domain in CATH-B dbs try CATH dbs
        if (family == "") or (family == k1):
            # return only domains from wanted family,
            # if there is no wanted family return alldomains
            for k1, v1 in data[pdbid][CATH].items():
                for v2 in v1['mappings']:
                  if domain == "":
                    # we want all domains in the chain
                    if v2["chain_id"] == chain:
                        out[v2['segment_id']] = (v2['start']['residue_number'], v2['end']['residue_number'])
                  else:
                    if v2['domain'] == pdbid + chain + domain:
                        out[v2['segment_id']] = (v2['start']['residue_number'], v2['end']['residue_number'])
    out = [out[i] for i in sorted(out.keys())]
    return out


def get_domain_boundaries(protein, chain, domain, family=""):
    # find domain boundaries of given protein chain domain tuple
    sequence = get_domain_ranges(protein, chain, domain, family)
    ret = ','.join(['%s:%s' % borders for borders in sequence])
    return ret


