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
        <p> <small> <small> Comment: {{ comment }} </small> </small> </p>
        <p> <small> <small> Number of chains: {{ chains }} </small> </small> </p>
        <p> <small> <small> Timestamp: {{ timestamp }} </small> </small> </p>
    </td>
    <td>
        <p> <small> <small> pLDDT score diagram: </small> </small> </p>
        <img src='{{ image_path }}' width='240'/>
    </td>
    <td>
        <p> <small> <small> AlphaMissense: </small> </small> </p>
        <img src='{{ image_am }}' width='340'/>
    </td>
    <td>
        <p> <small> <small> charges: </small> </small> </p>
        <img src='{{ image_charges }}' width='340'/>
    </td>
    <td>
        <p> <small> <small> 2DProts version: </small> </small> </p>
        <img src='https://2dprots.ncbr.muni.cz/files/family/{{ protein }}/multiple4.svg' width='240'/>
    </td>
    <td>
        <p> <small> <small> www.cathdb.info 3D: </small> </small> </p>
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


def add_family(status_web, web_dir, directory, family):
    """Add info from directory into status web"""
    comm = family.find("#")
    if comm != -1:
        family_name = family[:comm].strip()
    else:
        family_name = family.strip()

    chains = 0
    ts = "--"
    for d in os.listdir(directory):
        if d.startswith(f'generated-{family_name}-'):
            # the right family dir
            old_image_path1 = f"{directory}/{d}/{family_name}_sses.png"
            new_image_path1 = f"{web_dir}/images/{family_name}_sses.png"
            try:
                shutil.copyfile(old_image_path1, new_image_path1)
            except:
                logging.error(f"file {old_image_path1} does not exist")
            old_image_path2 = f"{directory}/{d}/{family_name}_am.png"
            new_image_path2 = f"{web_dir}/images/{family_name}_am.png"
            try:
                shutil.copyfile(old_image_path2, new_image_path2)
            except:
                logging.error(f"file {old_image_path2} does not exist")
            old_image_path3 = f"{directory}/{d}/{family_name}_charges.png"
            new_image_path3 = f"{web_dir}/images/{family_name}_charges.png"
            try:
                shutil.copyfile(old_image_path3, new_image_path3)
            except:
                logging.error(f"file {old_image_path3} does not exist")
            ts_index = d[12:].find("-")
            ts = d[ts_index+12+1:]
            chains = 0
            for f in os.listdir(f"{directory}/{d}"):
                if f.endswith(f'l.svg'):
                    chains += 1
            break

    data = {
        'protein': family_name,
        'chains': int(chains/2),
        'timestamp': ts,
        'failures': False,
        'image_path': f"images/{family_name}_sses.png",
        'image_am': f"images/{family_name}_am.png",
        'image_charges': f"images/{family_name}_charges.png",
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
                        help='set of output families')
    parser.add_argument('--debug', action='store_true',
                        help='show debug output')
    parser.add_argument('--random-sample', default=10,
                        help="show <given number> of random samples (default 10) could be 'all' as well")
    args = parser.parse_args()

    # By default, set logging to INFO
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logging.debug(args)

    # Read all input proteins
    testset = []
    if args.test_set is not None:
        name = args.test_set
        with open(args.test_set, "r") as fp:
            family = fp.readline()
            while family:
                testset.append(family)
                family = fp.readline()
    else:
        name = f"rand{args.random_sample}"
        f_list = []
        for d in os.listdir(args.directory):
            if d.startswith(f'generated-'):
                s = d.split("-")
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
    for family in testset:
        print("family", family)
        add_family(status_web_index, web_dir, args.directory, family)


if __name__ == '__main__':
    main()
