#!/usr/bin/env python3
# example:
# ./generate_protein_comparison_web.py --testset testsets/prot16 ../generate_protein_entry
# create html page conteining for sample testset of proteins
# 2dprot diagrams from given directory and relevant 3d models from
# pdbe and cath databases

import logging
import argparse
import os
import shutil
import requests
import time
import random
from jinja2 import Template


INDEX = "index.html"
HEADER_TEMPLATE = """
<html>
    <head>
        <title>2DProt status page</title>
        <style>
            table, td, th {
                border: 1px solid black;
            }

            table {
                border-collapse: collapse;
                width: 90%;
            }

            td {
                padding: .5em;
            }
        </style>
    </head>
    <body>
        <table>
"""

RECORD_TEMPLATE = """
<tr>
    <td>
        <p> <small> <small> Protein: {{ protein }} </small> </small> </p>
        <p> <small> <small> Pocet chainu: {{ chains }} </small> </small> </p>
        <p> <small> <small> Timestamp: {{ timestamp }} </small> </small> </p>
    </td>
    <td>
        <p> <small> <small> 2DProt </small> </small> </p>
        <img src='{{ image_path }}' width='400'/>
    </td>
    <td>
        <p> <small> <small> www.ebi.ac.uk/pdbe/ front: </small> </small> </p>
        <img src='http://www.ebi.ac.uk/pdbe/static/entry/{{ protein }}_deposited_chain_front_image-800x800.png' width='240'/>
    </td>
    <td>
        <p> <small> <small> www.ebi.ac.uk/pdbe/ side: </small> </small> </p>
        <img src='http://www.ebi.ac.uk/pdbe/static/entry/{{ protein }}_deposited_chain_side_image-800x800.png' width='240'/>
    </td>
    <td>
        <p> <small> <small> www.cathdb.info: </small> </small> </p>
        <img src='http://www.cathdb.info/version/v4_2_0/molscript/{{ protein }}/L?colour=ss' width='240'/>
    </td>
</tr>
"""


def get_hera_url(protein=None, chain=None, domain=None):
    if protein is not None and chain is not None and domain is not None:
        url = "http://www.ebi.ac.uk/thornton-srv/databases/PDBsum/%s/%s//dom%s%s.gif" % (protein[1:3], protein, chain, domain)
        r = requests.head(url)
        if r.status_code == 404:
            url = "http://www.ebi.ac.uk/thornton-srv/databases/PDBsum/%s/%s//dom%s%s.gif" % (protein[1:3], protein, 'A', domain)
        return url
    else:
        return ''


def add_dir(status_web, web_dir, directory, protein):
    """Add info from directory into status web"""
    for d in os.listdir(f"{directory}/{protein[0]}/{protein[1]}"):
        if d.startswith(f'generated-{protein}-'):
            old_image_path = f"{directory}/{protein[0]}/{protein[1]}/{d}/{protein}_l.svg"
            new_image_path = f"{web_dir}/images/{protein}l.svg"
            chains = 0
            try:
                shutil.copyfile(old_image_path, new_image_path)
                ts = d[15:]
                chains = ""
            except:
                ts = "----"

    data = {
        'protein': protein,
        'chains': chains,
        'timestamp': ts,
        'failures': False,
        'image_path': f"images/{protein}l.svg"
    }
    logging.debug('Adding directory record with %s' % data)
    with open(status_web, 'a') as fp:
        fp.write(Template(RECORD_TEMPLATE).render(data))


def main():
    parser = argparse.ArgumentParser(
        description="Generate web page which comparises input 2d pictures and 3d pictures of proteins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('directory', action='store',
                        help='directory which contains test protein 2d diagrams')
    parser.add_argument('--output', default='results',
                        help='path to the directory with output results')
    parser.add_argument('--test-set', action='store',
                        help='set of output proteins')
    parser.add_argument("--debug", action="store_true",
                        help="show debug output")
    parser.add_argument("--random-sample", action="store",
                        help="show given number of random proteins (default 10), could be 'all' as well")
    args = parser.parse_args()

    # By default, set logging to INFO
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logging.debug(args)

    #Read all input proteins
    testset = []
    if args.test_set is not None:
        name = args.test_set
        with open(args.test_set,"r") as fp:
            protein = fp.readline()
            while protein:
                testset.append(protein)
                protein = fp.readline()
    else:
        name = f"rand{args.random_sample}"
        f_list = []
        for d1 in os.listdir(args.directory):
            if len(d1) == 1:
                for d2 in os.listdir(f"{args.directory}/{d1}"):
                    if len(d2) == 1:
                        for d3 in os.listdir(f"{args.directory}/{d1}/{d2}"):
                            if d3.startswith("generated-"):
                                s = d3.split("-")
                                f_list.append(s[1])
        if args.random_sample == "all":
            testset = f_list
        else:
            testset = random.sample(f_list, int(args.random_sample))

    logging.debug("testset", testset)

    # If directory does not exist, create it
    timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
    f = name.rfind("/")
    if f != -1:
        fname = name[f+1:]
    else:
        fname = name
    web_dir = f"{args.output}/testset-{fname}-{timestamp}"

    os.mkdir(web_dir)
    os.mkdir(f"{web_dir}/images")
    logging.debug("Created directory %s" % web_dir)

    # If index page is empty, initiate it
    status_web_index = f"{web_dir}/index.html"
    with open(status_web_index, 'w') as fp:
        fp.write(Template(HEADER_TEMPLATE).render())
        logging.debug("Created file %s" % status_web_index)

    # Add one row to the table per directory
    for protein in testset:
        add_dir(status_web_index, web_dir, args.directory, protein[0:4])


if __name__ == '__main__':
    main()
