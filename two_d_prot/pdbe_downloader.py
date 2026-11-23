#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import logging
import os.path
import requests
import time


class PDBeDownloadError(Exception):
    """Generic PDBe Download error"""
    pass


class PDBeProteinMissing(Exception):
    """Downloading protein *.cif file from PDBe returned 404"""
    pass


class PDBeDownloader(object):

    PDBE_PROTEIN_URL = "http://www.ebi.ac.uk/pdbe/entry-files/download/%s.cif"
    PDBE_PROTEIN_URL_PDB = "http://www.ebi.ac.uk/pdbe/entry-files/download/pdb%s.ent"

    def protein(self, protein, dest_dir, timeout=100, retries=10):
        url = self.PDBE_PROTEIN_URL % protein
        dest_file = "%s/%s.cif" % (dest_dir, protein)
        if os.path.isfile(dest_file):
            logging.debug("File %s already exists, skipping it" % dest_file)
            return dest_file
        logging.debug("Downloading %s into %s ..." % (url, dest_dir))

        attempts = 0
        attempts_max = retries
        while True:
            try:
                response = requests.get(url, stream=True, timeout=timeout)
            except (requests.exceptions.ConnectTimeout,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.ReadTimeout) as e:
                logging.debug("Failed to get %s: %s" % (url, e))
                attempts += 1

                if attempts >= attempts_max:
                    raise PDBeDownloadError("Out of tries to download %s" % url)
                else:
                    time.sleep(10)
                    continue

            if response.status_code == 404:
                raise PDBeProteinMissing("Protein %s missing on %s" % (protein, url))

            with open(dest_file, "wb") as fp:
                for data in response.iter_content(chunk_size=4096):
                    fp.write(data)
                break

        logging.debug("Downloaded %s B into %s" % (os.path.getsize(dest_file), dest_file))
        return dest_file

    def protein_pdb(self, protein, dest_dir, timeout=100, retries=10):
        url = self.PDBE_PROTEIN_URL_PDB % protein
        dest_file = "%s/%s.pdb" % (dest_dir, protein)
        if os.path.isfile(dest_file):
            logging.debug("File %s already exists, skipping it" % dest_file)
            return dest_file
        logging.debug("Downloading %s into %s ..." % (url, dest_dir))

        attempts = 0
        attempts_max = retries
        while True:
            try:
                response = requests.get(url, stream=True, timeout=timeout)
            except (requests.exceptions.ConnectTimeout,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.ReadTimeout) as e:
                logging.debug("Failed to get %s: %s" % (url, e))
                attempts += 1

                if attempts >= attempts_max:
                    raise PDBeDownloadError("Out of tries to download %s" % url)
                else:
                    time.sleep(10)
                    continue

            if response.status_code == 404:
                raise PDBeProteinMissing("Protein %s missing on %s" % (protein, url))

            with open(dest_file, "wb") as fp:
                for data in response.iter_content(chunk_size=4096):
                    fp.write(data)
                break

        logging.debug("Downloaded %s B into %s" % (os.path.getsize(dest_file), dest_file))
        return dest_file
