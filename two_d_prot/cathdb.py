import logging
import os
import errno
import gzip
import shutil
import urllib.error
import urllib.request
import contextlib


def create_dir_if_missing(directory):
    try:
        logging.debug("Creating directory %s" % directory)
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def download_file_if_missing(url, path):
    if not os.path.isfile(path):
        logging.debug("Downloading %s to %s" % (url, path))
        try:
            with contextlib.closing(urllib.request.urlopen(url)) as r:
                with open(path, 'wb') as f:
                    shutil.copyfileobj(r, f)
        except urllib.error.URLError as e:
            if '550 Failed to change directory.' in e.reason:
                logging.error("Failed to download %s because it is not present on CATH servers" % url)
                raise FileNotFoundError(e.reason)
            else:
                raise


class CathDailyReleaseDomain(object):

    def __init__(self, name, family, boundaries, release=None):
        self._name = name.decode()
        self._family = family.decode()
        self._boundaries = boundaries.decode()
        self._release = release.decode()

    def __key(self):
        return (self._name, self._family, self._boundaries)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __str__(self):
        return ' '.join(self.__key())

    @property
    def family(self):
        return self._family

    @property
    def name(self):
        return self._name

    @property
    def protein(self):
        return self._name[:4]

    @property
    def chain(self):
        return self._name[4:][:-2]

    @property
    def domain(self):
        return self._name[-2:]

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def release(self):
        return self._release


class CathDailyRelease(object):

    STORAGE = 'input/cath-daily-release/'

    def __init__(self, url):
        self.file_name  = os.path.basename(url)
        self.domains = set()

        # Check daily release URL is what we expect
        # e.g. ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/archive/cath-b-20190725-all.gz
        assert (self.file_name.endswith('-all.gz') or self.file_name.endswith('-putative.gz')) and self.file_name.startswith('cath-b-')

        # Create storage directory if it does not exist
        create_dir_if_missing(self.STORAGE)

        # Download release file if we do not have it already
        self.file_path = os.path.join(self.STORAGE, self.file_name)
        download_file_if_missing(url, self.file_path)

        # Unpack the file
        logging.debug("Unpacking file %s" % self.file_path)
        with gzip.open(self.file_path, 'rb') as f:
            self._parse_file(f.read())

    def _parse_file(self, content):
        for row in content.split(b'\n'):
            if row == b'':
                continue
            row_split = row.split()
            try:
                domain = CathDailyReleaseDomain(row_split[0], row_split[2], row_split[3], row_split[1])
            except IndexError:
                logging.warning("Row %s parsing failed" % str(row))
            self.domains.add(domain)

    def get_domains_for_family(self, family):
        out = []
        for d in self.domains:
            if d.family == family:
                out.append(d)
        return out

    def get_family_dor_domain(self, protein, chain, domain):
        for d in self.domains:
            if d.protein == protein and d.chain == chain and d.domain == domain:
                return d.family

    @property
    def families_putative(self):
        """Return set of families with putative release"""
        out = set()
        for d in self.domains:
            if d.release is None or d.release == 'putative':
                out.add(d.family)
        return out

    @property
    def families(self):
        """Return set of families"""
        out = set()
        for d in self.domains:
            out.add(d.family)
        return out


class CathVersionedRelease(object):

    STORAGE = 'input/cath-versioned-release/'

    STORAGE_CATH_CLASSIFICATION_DATA = 'cath-classification-data'
    FILE_CATH_DOMAIN_LIST = 'cath-domain-list-v4_3_0.txt'

    def __init__(self, url):
        self.base_url = url   # e.g. http://download.cathdb.info/cath/releases/all-releases/v4_2_0/
        self.version = [i for i in self.base_url.split('/') if i != ''][-1]
        self.cache_dir = os.path.join(self.STORAGE, self.version)

        # Check we have guessed version correctly
        # (well, this is pretty weak check, but better than nothing)
        assert self.version.startswith('v')

        # Create storage directory if it does not exist
        create_dir_if_missing(self.cache_dir)

        self._cath_domain_list = None

    @property
    def cath_domain_list(self):
        if self._cath_domain_list is None:
            self._cath_domain_list = {}

            # Create caching directory to save the file if needed
            directory = os.path.join(self.cache_dir, self.STORAGE_CATH_CLASSIFICATION_DATA)
            create_dir_if_missing(directory)

            # Download file if needed
            url = os.path.join(self.base_url, self.STORAGE_CATH_CLASSIFICATION_DATA, self.FILE_CATH_DOMAIN_LIST)
            path = os.path.join(directory, self.FILE_CATH_DOMAIN_LIST)
            download_file_if_missing(url, path)

            # Parse file
            logging.debug("Parsing file %s" % path)
            with open(path, 'r') as f:
                for row in f.read().split('\n'):
                    if row == '' or row[0] == '#':
                        continue
                    row_split = row.split()
                    try:
                        domain = (row_split[0][0:4], row_split[0][4:][:-2], row_split[0][-2:])
                        family = row_split[1:8]
                    except IndexError:
                        logging.warning("File %s row '%s' parsing failed" % (path, str(row)))
                    else:
                        self._cath_domain_list[domain] = family
            logging.debug("Parsed file %s with %s entries" % (path, len(self._cath_domain_list)))

        return self._cath_domain_list


    def get_familyS35_for_domain(self, protein, chain, domain):
        """Return S35 family for requested domain"""
        domain = (protein, chain, domain)
        return '.'.join(self.cath_domain_list[domain][:5])
