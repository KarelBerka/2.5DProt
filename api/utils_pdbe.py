#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import requests
import logging
import time
import json
import os
import utils_common

API_URL = 'https://www.ebi.ac.uk/pdbe/api/'
PROTEIN_URL = API_URL + '/mappings/'
CATH = 'CATH'

def get_protein_pdbe_domain_numbers(protein):
    # return list/number (if number_only is set)
    # of all cath domains of the given protein

    protein_filename = PROTEIN_URL + protein
    data = utils_common.get_url_to_json(protein_filename)
    str_data=""
    str_number = 0
    if (data != None) and (protein in data) and \
        CATH in data[protein]:
        for i in data[protein][CATH]:
            for j in data[protein][CATH][i]['mappings']:
                str_data += f"{j['domain']},"
                str_number += 1
        return str_number
    else:
        return ""
