#!/usr/bin/env python3
# example:
# ./generate_protein_family_inside_comparison_web.py  ../generate_protein_entry
# create html page conteining family 2dprot diagrams from given directory and
# relevant 3d models from pdbe and cath databases

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
    <center>
    <h2> 2DProt visualisation of CATH families generated from AlphaFold database</h2>
    <h3>(October 24, 2022 snapshot)</h3>
    </center>
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
         <small>CATH family:</small>
         <break>
         <a href="{{family}}/index.html">{{family}}</a>
    </td>
    <td>
        <p> <small> <small> AlphaFold dbs - opaque </small> </small> </p>
        <img src='{{ image_path1 }}' width='400'/>
    </td>
    <td>
        <p> <small> <small> AlphaFold dbs - transparent with template {{ family}}</small> </small> </p>
        <img src='{{ image_path2 }}' width='400'/>
    </td>
    <td>
        <p> <small> <small> pdbe dbs - transparent with template </small> </small> </p>
        <img src='https://2dprots.ncbr.muni.cz/files/family/{{ family }}/multiple2.png' width='240'/>
    </td>
    <td>
        <p> <small> <small> www.cathdb.info </small> </small> </p>
        <img src='http://www.cathdb.info/version/v4_2_0/molscript/{{ family }}/L?colour=ss' width='240'/>
    </td>
</tr>
"""

INDEX_TEMPLATE = """
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
    <center>
    <h2> 2DProt visualisation of CATH families generated from AlphaFold database</h2>
    <h3>(October 24, 2022 snapshot)</h3>
    </center>
    <ul>
"""

INDEX_RECORD_TEMPLATE = """
    <li> {{ prefix }} <a href="{{ cluster_path }}">{{ cluster }}</a> {{ text }}</li>
"""

FAMILY_TEMPLATE = """
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
    <center>
    <h2> 2DProt visualisation of CATH families generated from AlphaFold database</h2>
    <h3>(October 24, 2022 snapshot)</h3>
    </center>
    <table>
    <tr>
    <td>
         <small>CATH family:</small>
         <break>
         <a href="{{family}}/index.html">{{family}}</a>
    </td>
    <td>
        <p> <small> <small> AlphaFold dbs - opaque </small> </small> </p>
        <img src='{{ image_path1 }}' width='400'/>
    </td>
    <td>
        <p> <small> <small> AlphaFold dbs - transparent with template {{ family}}</small> </small> </p>
        <img src='{{ image_path2 }}' width='400'/>
    </td>
    <td>
        <p> <small> <small> pdbe dbs - transparent with template </small> </small> </p>
        <img src='https://2dprots.ncbr.muni.cz/files/family/{{ family }}/multiple2.svg' width='240'/>
    </td>
    <td>
        <p> <small> <small> www.cathdb.info 3D diagram</small> </small> </p>
        <img src='http://www.cathdb.info/version/v4_2_0/molscript/{{ family }}/L?colour=ss' width='240'/>
    </td>
</tr>
    <tr><td>-</td><td>-</td><td>-</td><td>-</td><td>-</td></tr>
    <tr>
"""

FAMILY_RECORD_TEMPLATE = """
    <td>
        <p> <small> <small> {{domain_name}} </small> </small> </p>
        <img src='{{ domain_path }}' width='240'/>
    </td>
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


def read_name_data(fd):
    ret = {1: {},
           2: {},
           3: {},
           }
    line = fd.readline()
    name1 = None
    name2 = None
    name3 = None
    while (line):
        c = line.count("\t")
        if (line.count("]") == 0) and (line.count("[") == 0):
            # omit braces
            if c == 2:
                # the first level
                sline = line.split('"')
                if name1 is None:
                    name1 = sline[1]
                else:
                    ret[1][name1] = sline[1]
                    name1 = None
            if c == 4:
                # the second level
                sline = line.split('"')
                if name2 is None:
                    name2 = sline[1]
                else:
                    ret[2][name2] = sline[1]
                    name2 = None
            if c == 6:
                # the second level
                sline = line.split('"')
                if name3 is None:
                    name3 = sline[1]
                else:
                    ret[2][name3] = sline[1]
                    name3 = None

        line = fd.readline()
    return ret


def add_protein_family(web_dir, directory, family):
    """Add info from directory into status web"""
    multi = None
    tree = set()

    for e in os.listdir(f"{directory}"):
        for f in os.listdir(f"{directory}/{e}"):
            if f.endswith(".svg") and f.startswith("image"):
                sp = e.split("-")
                family = sp[1]
                fsplit = family.split(".")
                # create all directory hierarchy
                if not (os.path.isdir(f"{web_dir}/{fsplit[0]}")):
                    os.mkdir(f"{web_dir}/{fsplit[0]}")
                if not (os.path.isdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}")):
                    os.mkdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}")
                if not (os.path.isdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}")):
                    os.mkdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}")
                if not os.path.isdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{family}"):
                    os.mkdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{family}")
                old_image_path = f"{directory}/{e}/{f}"
                aux_f = f.replace(":", "_")
                new_image_path = f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{family}/{aux_f}"

                shutil.copyfile(old_image_path, new_image_path)

            if f.endswith(".png") and not f.endswith("_2.png") and not f.endswith("_3.png") and not f.endswith("_4.png"):
                # this is per family directory add them to a tree, add the record to the list
                fsplit = f.split(".")
                tree.add(f)
                # create all directory hierarchy
                if not (os.path.isdir(f"{web_dir}/{fsplit[0]}")):
                    os.mkdir(f"{web_dir}/{fsplit[0]}")
                if not (os.path.isdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}")):
                    os.mkdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}")
                if not (os.path.isdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}")):
                    os.mkdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}")
                if not os.path.isdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{f[:-4]}"):
                    os.mkdir(f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{f[:-4]}")
                status_web_index = f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/index.html"
                if not (os.path.isfile(status_web_index)):
                    with open(status_web_index, 'w') as fp:
                        fp.write(Template(HEADER_TEMPLATE).render())
                        logging.debug("Created file %s" % status_web_index)

                old_image_path = f"{directory}/{e}/{f}"
                new_image_path = f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{f}"
                shutil.copyfile(old_image_path, new_image_path)
#                old_image_path = f"{directory}/{e}/{f[:-4]}2.png"
#                new_image_path = f"{web_dir}/{fsplit[0]}/{fsplit[1]}/{fsplit[2]}/{f[:-4]}2.png"
#                try:
#                    shutil.copyfile(old_image_path, new_image_path)
                multi = f
                # adds its single page

        if not multi:
            return
        with open(status_web_index, 'a') as fp:
                data = {
                    'image_path1': f"{multi}",
                    'image_path2': f"{multi[:-5]}.png",
                    'family': multi[:-10],
                }
                fp.write(Template(RECORD_MULTIPLE_TEMPLATE).render(data))
                counter = 0

    # create main page for all dbs
    index_page0 = f"{web_dir}/index.html"

    # read clusters name list file cl_names.json
    # is from overprot tool
    with open("cl_names.py", 'r') as fp:
        name_data = read_name_data(fp)

    with open(index_page0, 'w') as fp:
        fp.write(Template(INDEX_TEMPLATE).render())
        logging.debug("Created file %s" % index_page0)

    for e in os.listdir(f"{web_dir}"):
        if (e == "index.html"):
            continue
        with open(index_page0, 'a') as fp:
            data = {
                'prefix': "",
                'cluster_path': f"{e}/index.html",
                'cluster': e,
                'text': name_data[1][e+"."]
            }
            fp.write(Template(INDEX_RECORD_TEMPLATE).render(data))

        index_page1 = f"{web_dir}/{e}/index.html"
        with open(index_page1, 'w') as fp:
            fp.write(Template(INDEX_TEMPLATE).render())
            logging.debug("Created file %s" % index_page1)

        # create second level main page for all dbs
        for f in os.listdir(f"{web_dir}/{e}"):
            if (f == "index.html"):
                continue
            with open(index_page1, 'a') as fp:
                aux = name_data[2][f"{e}.{f}."]
                data = {
                    'prefix': f"{e}.",
                    'cluster_path': f"{f}/index.html",
                    'cluster': f,
                    'text': aux,
                }
                fp.write(Template(INDEX_RECORD_TEMPLATE).render(data))

            index_page2 = f"{web_dir}/{e}/{f}/index.html"
            with open(index_page2, 'w') as fp:
                fp.write(Template(INDEX_TEMPLATE).render())
                logging.debug("Created file %s" % index_page1)

            # create third level main page for all dbs
            for g in os.listdir(f"{web_dir}/{e}/{f}"):
                if (g == "index.html"):
                    continue
                with open(index_page2, 'a') as fp:
                    print("g", g)
                    data = {
                        'prefix': f"{e}.{f}.",
                        'cluster_path': f"{g}/index.html",
                        'cluster': g,
                        'text': name_data[2][f"{e}.{f}.{g}."],
                    }
                    fp.write(Template(INDEX_RECORD_TEMPLATE).render(data))

                for h in os.listdir(f"{web_dir}/{e}/{f}/{g}"):
                    if (h.endswith(".png")) or (h.endswith("html")):
                        continue
                    # create per family pages
                    family_page = f"{web_dir}/{e}/{f}/{g}/{h}/index.html"
                    with open(family_page, 'w') as fp:
                        data = {
                            'family': h,
                            'image_path1': f"../{h}.png",
                            'image_path2': f"../{h}_2.png",
                        }
                        fp.write(Template(FAMILY_TEMPLATE).render(data))
                        logging.debug("Created file %s" % index_page1)
                    print("family_page", family_page)
                    counter = 0

                    for i in os.listdir(f"{web_dir}/{e}/{f}/{g}/{h}"):
                        if i.endswith(".svg") and i.startswith("image"):
                            counter += 1
                            with open(family_page, 'a') as fp:
                                data = {
                                    'domain_name': i[6:-4],
                                    'domain_path': i,
                                    'cluster': g,
                                    'text': "",
                                }
                                fp.write(Template(FAMILY_RECORD_TEMPLATE).render(data))
                                if counter % 5 == 0:
                                    fp.write("</tr><tr>")
                    while counter % 5 != 0:
                        with open(family_page, 'a') as fp:
                            fp.write("</td><td>")
                        counter += 1


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
    logging.debug("Created directory %s" % web_dir)

    add_protein_family(web_dir, args.directory, name)


if __name__ == '__main__':
    main()
