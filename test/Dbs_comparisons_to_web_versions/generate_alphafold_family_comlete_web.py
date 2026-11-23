#!/usr/bin/env python3
# example:
# ./generate_alphafold_family_comlete_web.py --dir <dirname>
# create html page conteining alphafolf family 2dprot multi diagrams 
# and 2dprotsingle diagrams for all domains from given directory and
# relevant pdb  model from 2dprot and 3D CATH databases

import logging
import argparse
import os
import shutil
import requests
import time
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
    <td>
        <p> <small> <small> {{ protein }} </small> </small> </p>
        <img src='{{ image_path }}' width='400'/>
    </td>
"""

RECORD_MULTIPLE_TEMPLATE = """
<tr>
    <td>
        <p> <small> <small> 2DProt - {{ family}} </small> </small> </p>
        <img src='{{ image_path1 }}' width='400'/>
    </td>
        <td>
        <p> <small> <small> 2DProt - {{ family}}</small> </small> </p>
        <img src='{{ image_path2 }}' width='400'/>
    </td>
    <td>
        <p> <small> <small> 2DProts cath - {{family}} </small> </small> </p>
        <img src='https://2dprots.ncbr.muni.cz/files/family/{{ family }}/multiple4.svg' width='240'/>
    </td>
    <td>
        <p> <small> <small> www.cathdb.info: {{family}} </small> </small> </p>
        <img src='http://www.cathdb.info/version/v4_2_0/molscript/{{ family }}/L?colour=ss' width='240'/>
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


def add_protein_family(status_web, web_dir, directory, family):
    """Add info from directory into status web"""
    single_list = []
    chains = 0
    multi = None

    print("directory", directory)
    for f in os.listdir(f"{directory}"):
        if f.endswith(".svg") and f.startswith("image"):
            old_image_path = f"{directory}/{f}"
            new_image_path = f"{web_dir}/images/{f}"
            shutil.copyfile(old_image_path, new_image_path)
            chains += 1
            single_list.append(f)
            print("single", f)
        if f.endswith(".png") and not f.endswith("_2.png") and not f.endswith("_3.png") and not f.endswith("_4.png"):
            old_image_path = f"{directory}/{f}"
            new_image_path = f"{web_dir}/images/{f}"
            shutil.copyfile(old_image_path, new_image_path)
            old_image_path = f"{directory}/{f[:-4]}_2.png"
            new_image_path = f"{web_dir}/images/{f[:-4]}_2.png"
            shutil.copyfile(old_image_path, new_image_path)
            old_image_path = f"{directory}/{f[:-4]}_3.png"
            new_image_path = f"{web_dir}/images/{f[:-4]}_3.png"
            shutil.copyfile(old_image_path, new_image_path)
            old_image_path = f"{directory}/{f[:-4]}_4.png"
            new_image_path = f"{web_dir}/images/{f[:-4]}_4.png"
            shutil.copyfile(old_image_path, new_image_path)
            multi = f
            print("multi", multi)

    if not multi:
        return
    with open(status_web, 'a') as fp:
            data = {
                'image_path1': f"images/{multi}",
                'image_path2': f"images/{multi[:-4]}_2.png",
                'image_path3': f"images/{multi[:-4]}_3.png",
                'image_path4': f"images/{multi[:-4]}_4.png",
                'family': multi[:-4],
            }
            fp.write(Template(RECORD_MULTIPLE_TEMPLATE).render(data))

            fp.write("<tr><td>-</td><td>-</td><td>-</td><td>-</td></tr><tr>")
            counter = 0
            for f in single_list:
                data = {
                    'protein': f[6:],
                    'image_path': f"images/{f}"
                }
                fp.write(Template(RECORD_TEMPLATE).render(data))
                counter += 1
                if counter % 4 == 0:
                    fp.write("</tr><tr>")
                logging.debug('Adding directory record with %s' % data)
            fp.write("</tr>")


def main():
    parser = argparse.ArgumentParser(
        description="Generate web page which comparises input 2d pictures and 3d pictures of proteins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('directory', action='store',
                        help='directory which contains test protein 2d diagrams')
    parser.add_argument('--output', default='results',
                        help='path to the directory with output results')
    parser.add_argument("--debug", action="store_true",
                        help="show debug output")
    args = parser.parse_args()

    # By default, set logging to INFO
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logging.debug(args)

    # If directory does not exist, create it
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

    # If index page is empty, initiate it
    status_web_index = f"{web_dir}/index.html"
    with open(status_web_index, 'w') as fp:
        fp.write(Template(HEADER_TEMPLATE).render())
        logging.debug("Created file %s" % status_web_index)
    add_protein_family(status_web_index, web_dir, args.directory, name)


if __name__ == '__main__':
    main()
