#!/usr/bin/env python3
# example:
# ./s35_comparison_in_one_family.py ../generate_family_entry
# TODO: add comments

import logging
import argparse
import os
import shutil
import requests
import time
from jinja2 import Template
import json


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
        <p> <small> <small> Protein: {{ s35_family }} </small> </small> </p>
        <p> <small> <small> Number of chains: {{ chains }} </small> </small> </p>
        <p> <small> <small> Percentage of chains with channels: {{ channels }} </small> </small> </p>
        <p> <small> <small> Average number of channels in chins with channels: {{ channels_num }} </small> </small> </p>
        <p> <small> <small> Channels type: {{ channels_type }} </small> </small> </p>
    </td>
    <td>
        <p> <small> <small> 2DProt </small> </small> </p>
        <img src='{{ image_path4 }}' width='240'/>
    </td>
    <td>
        <p> <small> <small> 2DProt </small> </small> </p>
        <img src='{{ image_path2 }}' width='340'/>
    </td>
"""


DOM_TEMPLATE = """
    <td>
        <p> <small> <small> domain {{f}} </small> </small> </p>
        <p> <small> <small> Channels type: {{ channels_type }} </small> </small> </p>
        <img src='{{ image_path }}' width='340'/>
    </td>
"""


END_TEMPLATE = """
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


def add_s35_family(status_web, web_dir, directory, s35_family):
    """Add info from directory into status web"""

    chains = 0
    occ = s35_family.rfind(".")
    fam = s35_family[:occ]
    print("fam", fam, s35_family)
    # the right family dir
    old_image_path1 = f"{directory}/only_all.png"
    new_image_path1 = f"{web_dir}/images/{s35_family}_2.png"
    try:
        shutil.copyfile(old_image_path1, new_image_path1)
    except:
        logging.error(f"file {old_image_path1} does not exist")
    old_image_path2 = f"{directory}/only_channels.png"
    new_image_path2 = f"{web_dir}/images/{s35_family}_4.png"
    try:
        shutil.copyfile(old_image_path2, new_image_path2)
    except:
        logging.error(f"file {old_image_path1} does not exist")
    chains = 0
    channels = 0
    channels_num = 0
    channels_type = set()

    for f in os.listdir(f"{directory}"):
        if f.endswith('.json') and f.startswith("layout-"):
            # for each layout
            chains += 1
            with open(f"{directory}/{f}", 'r') as fp:
                data = json.load(fp)
                print("name ", f"{directory}/{f}")
                chan = set()
                if "channels" in data and (len(data["channels"]) != 0):
                    channels += 1
                    for ch in data["channels"]:
                        ind = ch.rfind(".")
                        if ch[:ind] not in chan:
                            chan.add(ch[:ind])
                            channels_num += 1
                        if "text_names" in data["channels"][ch] and ("." not in data["channels"][ch]["text_names"]):
                            channels_type.add(data["channels"][ch]["text_names"])

    if channels_type == set() or channels_type == {"Tunnel"}:
        return

    if chains > 0:
        ratio = channels/chains
    else:
        ratio = 0

    if channels > 0:
        ch_ratio = channels_num/channels
    else:
        ch_ratio = 0

    data = {
        's35_family': s35_family,
        'chains': chains,
        'channels': ratio,
        'channels_num': ch_ratio,
        'failures': False,
        'image_path2': f"images/{s35_family}_2.png",
        'image_path4': f"images/{s35_family}_4.png",
        'channels_type': channels_type,
    }
    logging.debug('Adding directory record with %s' % data)
    with open(status_web, 'a') as fp:
        fp.write(Template(RECORD_TEMPLATE).render(data))

    for f in os.listdir(f"{directory}"):
        if f.endswith('.svg') and f.startswith("layout"):
            with open(f"{directory}/{f[:-4]}.json", 'r') as fp:
                channels_type = set()
                data = json.load(fp)
                if "channels" in data and (len(data["channels"]) != 0):
                    for ch in data["channels"]:
                        ind = ch.rfind(".")
                        if "text_names" in data["channels"][ch] and ("." not in data["channels"][ch]["text_names"]):
                            channels_type.add(data["channels"][ch]["text_names"])

            if channels_type == set() or channels_type == {"Tunnel"}:
                continue

            old_image_path1 = f"{directory}/{f}"
            new_image_path1 = f"{web_dir}/images/{f}"
            try:
                shutil.copyfile(old_image_path1, new_image_path1)
            except:
                logging.error(f"file {old_image_path1} does not exist")
            data = {'fname': f,
                    'image_path': f"images/{f}",
                    'channels_type': channels_type, }
            with open(status_web, 'a') as fp:
                fp.write(Template(DOM_TEMPLATE).render(data))

    with open(status_web, 'a') as fp:
        fp.write(Template(END_TEMPLATE).render(data))


def main():
    parser = argparse.ArgumentParser(
        description="Generate web page which comparises input 2d pictures and 3d pictures of proteins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('directory', action='store',
                        help='family directory which contains test protein 2d diagrams')
    parser.add_argument('--output', default='results',
                        help='path to the directory with output results')
    parser.add_argument('--debug', action='store_true',
                        help='show debug output')
    args = parser.parse_args()

    # By default, set logging to INFO
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logging.debug(args)

    # Create output directory
    timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
    name = args.directory
    f = name.rfind("/")
    if f != -1:
        fname = name[f+1:]
    else:
        fname = name
    web_dir = f"{args.output}/testset-{fname}-{timestamp}"
    os.mkdir(web_dir)
    os.mkdir(f"{web_dir}/images")
    logging.debug("Created directory %s" % web_dir)

    # Create initiate index page
    status_web_index = f"{web_dir}/index.html"
    with open(status_web_index, 'w') as fp:
        fp.write(Template(HEADER_TEMPLATE).render())
        logging.debug("Created file %s" % status_web_index)

    # Read all input proteins
    testset = []
    name = "all"
    for d in os.listdir(args.directory):
        full_d = f"{args.directory}/{d}"
        if os.path.isdir(full_d) and (d != "annotation") and (d != "2dprots") and (d != "layout2"):
            add_s35_family(status_web_index, web_dir, full_d, d)

    logging.debug("testset", testset)


if __name__ == '__main__':
    main()
