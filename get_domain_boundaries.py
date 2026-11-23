#!/usr/bin/env python3

import sys

import utils_pdbe_api

protein = sys.argv[1]
chain = sys.argv[2]
domain = sys.argv[3]

sequence = utils_pdbe_api.get_domain_ranges(protein, chain, domain)
print(','.join(['%s:%s' % borders for borders in sequence]))

if len(sequence) == 0:
    sys.exit(2)
else:
    sys.exit(0)
