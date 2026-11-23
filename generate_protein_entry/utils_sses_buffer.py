#!/usr/bin/env python3

import os
import hashlib
from pathlib import Path
import shutil
import filecmp


def is_file_in_buffer(protein, buffer_directory, orig):
    md5 = hashlib.md5(open(orig,'rb').read()).hexdigest()
    layout_filename = f"{buffer_directory}/{md5[0]}/{md5[1]}/{md5}/layout-{protein}.json"
    pdb_filename = f"{buffer_directory}/{md5[0]}/{md5[1]}/{md5}/{protein}.pdb"

    if os.path.isdir(buffer_directory):
        if os.path.isdir(f"{buffer_directory}/{md5[0]}"):
            if os.path.isdir(f"{buffer_directory}/{md5[0]}/{md5[1]}"):
                if os.path.isdir(f"{buffer_directory}/{md5[0]}/{md5[1]}/{md5}"):
                    if os.path.isfile(pdb_filename):
                        if filecmp.cmp(orig, pdb_filename):
                          if os.path.isfile(layout_filename):
                            return layout_filename

    return False


def save_layout_to_buffer(protein, buffer_directory, layout_original, pdb_original):
    md5 = hashlib.md5(open(pdb_original,'rb').read()).hexdigest()
    path = Path(f"{buffer_directory}/{md5[0]}/{md5[1]}/{md5}")
    path.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(layout_original, f"{buffer_directory}/{md5[0]}/{md5[1]}/{md5}/layout-{protein}.json")
    shutil.copyfile(pdb_original, f"{buffer_directory}/{md5[0]}/{md5[1]}/{md5}/{protein}.pdb")
    return True