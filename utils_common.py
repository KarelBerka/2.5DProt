#!/usr/bin/env python3

import os
import time
import logging
import subprocess
import requests
import json


def time_log(message, last_time):
    # create a log message containing a duration timestamp (used in generated log files
    # message   - log message describing this task
    # last_time - the timestamp of start of this task

    now_time = time.time()
    logging.error(f"--- Timestamp: {message} {time.strftime('%Y-%m-%dT%H-%M-%S')} {(now_time-last_time)/60} ---")
    return now_time


def error_report(directory, command, out, err, log_fn):
    # unified error report to log_fn file
    # directory - directory in which the command was used
    # command   - used command
    # out       - returned output
    # err       - returned error message
    # log_fn    - file name of the error log file

    with open(log_fn, "w") as fd:
        timestamp = time.strftime("%Y-%m-%dT%H_%M_%S")
        fd.write(f"* date: {timestamp}\n")
        fd.write(f"* commandline: {command}\n")
        fd.write(f"* dir: {directory}\n")
        fd.write("* stdout:\n")
        if out is not None:
            if type(out) == bytes:
                fd.write(out.decode("utf-8"))
            else:
                fd.write(out)
        else:
            fd.write("None\n")
        fd.write("* stderr:\n")
        if err is not None:
            if type(err) == bytes:
                fd.write(err.decode("utf-8"))
            else:
                fd.write(str(err))
        else:
            fd.write("None\n")


def git_version():
    # returns git version tag of the input

    command = "git log | head -n 1"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    r = p.wait()
    if (r == 0):
        out, err = p.communicate()
        commit_number = (out.decode("utf-8"))[:-1]
    else:
        commit_number = "missing"
    command = "git log | head -n 3 | tail -n 1"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    r = p.wait()
    if (r == 0):
        out, err = p.communicate()
        commit_date = (out.decode("utf-8"))[:-1]
    else:
        commit_date = "missing"
    return commit_number, commit_date


# grab url file using requests package
def get_url_res(url):
    # download file defined by input url using requests package

    attempt = 0
    attempt_max = 10
    while attempt < attempt_max:
        try:
            response = requests.get(url, timeout=30)
        except requests.exceptions.Timeout as e:
            logging.warning("Request to %s timed for a %s time: %s" % (url, attempt, e))
            if attempt >= attempt_max:
                logging.error("Failed the request to %s because it is timeouting" % url)
                raise
            attempt += 1
            time.sleep(1)
        else:
            break
    if response.ok:
        return response
    else:
        logging.error('HTTP request %s failed, status code %s' % (url, response.status_code))
        return response


def get_url_to_json(url):
    # download file defined by input url using requests package
    # and read it to json structure

    response = get_url_res(url)
    if response.ok is not None:
        data = json.loads(response.text)
    else:
        data = None
    return data


def get_url_to_file(url, dest_file):
    # download file defined by input url using requests package

    response = get_url_res(url)
    if response.ok:
        with open(dest_file, "wb") as fp:
            for data in response.iter_content(chunk_size=4096):
                fp.write(data)
        logging.debug("Downloaded %s B into %s" % (os.path.getsize(dest_file), dest_file))
        return dest_file
    else:
        dest_file = None
        return None


# AlphaFold functions
# def download_alphafold_list(directory, input_cif_dir, af_list):
#    # copy af_list set of AlphaFold files (filename AF-{af_list_entry}-*.cif)
#    # form input_cif_dir to directory
#    # af_list       - list of af entries to download
#    # input_cif_dir - directory containing af cifs
#    # directory     - the directory to which af entries will be coppied
#    # af_location="https://alphafold.ebi.ac.uk/files"
#
#    for afe in af_list:
#        found = False
#        for filename in os.listdir(input_cif_dir):
#            if filename.startswith(f"AF-{afe}-") and \
#               filename.endswith(".cif"):
#                shutil.copyfile(f"{input_cif_dir}/{filename}", f"{directory}/{filename}")
#                os.symlink(f"{filename}", f"{directory}/{afe}.cif")
#                shutil.copyfile(f"{input_cif_dir}/{filename[:-4]}.pdb", f"{directory}/{filename[:-4]}.pdb")
#                os.symlink(f"{filename[:-4]}.pdb", f"{directory}/{afe}.pdb")
#                found = True
#        if not found:
#            print(f"Alpha fold entry {afe} is not found in {input_cif_dir}")
#            logging.error(f"Alpha fold entry {afe} is not found in {input_cif_dir}")


def symlink_alphafold_list(directory, input_cif_dir, af_list):
    # symlink af_list set of AlphaFold files (filename AF-{af_list_entry}-*.cif)
    # form input_cif_dir to directory
    # af_list       - list of af entries to download
    # input_cif_dir - directory containing af cifs
    # directory     - the directory to which af entries will be symlinked
    # af_location="https://alphafold.ebi.ac.uk/files"

    for afe in af_list:
        found = False
        for filename in os.listdir(input_cif_dir):
            if filename.startswith(f"AF-{afe}-") and \
               filename.endswith(".cif"):
                os.symlink(f"../../{input_cif_dir}/{filename}", f"{directory}/{afe}.cif")
                os.symlink(f"../../{input_cif_dir}/{filename[:-4]}.pdb", f"{directory}/{afe}.pdb")
                found = True
        if not found:
            print(f"Alpha fold entry {afe} is not found in {input_cif_dir}")
            logging.error(f"Alpha fold entry {afe} is not found in {input_cif_dir}")


def wget_file(link, directory, filename, err_report=True, option=""):
    # using curl wget "link" file to given "filename"

    command = f"wget {option} {link} -O {filename} >/dev/null 2>/dev/null"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    r = p.wait()
    if (r != 0):
        if err_report is True:
            out, err = p.communicate()
            log_fn = f"{directory}/wget-fail-{time.time()}.txt"
            error_report(directory, command, out, err, log_fn)
            logging.error(f"curl of {link} failed, report save to {log_fn}")
        return False
    return True


def curl(directory, err_report=True, options=""):
    # use curl with arbitrary options

    command = f"curl {options} >/dev/null 2>/dev/null"
    counter = 0
    while (counter < 10):
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        r = p.wait()
        if (r != 0):
            counter += 1
        else:
            return True

    if err_report is True:
        out, err = p.communicate()
        log_fn = f"{directory}/curl-fail-{time.time()}.txt"
        error_report(directory, command, out, err, log_fn)
        logging.error(f"curl of {link} failed, report save to {log_fn}")
    return False


def curl_file(link, directory, filename, err_report=True, option=""):
    # using curl wget "link" file to given "filename"

    # If possible, use `download_file(link, directory, filename)` instead

    command = f"curl --connect-timeout 5  -L {option} -f {link} -o {filename} >/dev/null 2>/dev/null"
    counter = 0
    while (counter < 20):
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        r = p.wait()
        if (r != 0):
            counter += 1
        else:
            return True

    if err_report is True:
        out, err = p.communicate()
        log_fn = f"{directory}/curl-fail-{time.time()}.txt"
        error_report(directory, command, out, err, log_fn)
        logging.error(f"curl of {link} failed, report save to {log_fn}")
    return False


def download_file(link, directory, filename):
    """Downloads a file from a URL with retry logic and large file handling.

    Args:
        link: The URL of the file to download.
        directory: The directory to save the file to.
        filename: The name to save the file as.

    Returns:
        True on successful download, False otherwise.
    """

    retries = 20
    retry_delay = 1  # seconds

    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as e:
            logging.error(f"Error creating directory: {e}")
            return False

    for attempt in range(retries):
        try:
            with requests.get(link, stream=True, timeout=30) as response:
                response.raise_for_status()
                with open(filename, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192): # 8KB chunks
                        if chunk:
                            f.write(chunk)
            return True
        except requests.exceptions.RequestException as e:
            logging.warning(f"Download attempt {attempt + 1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(retry_delay)
        except OSError as e:
            logging.error(f"Error writing to file: {e}")
            return False

    logging.error(f"Download failed after {retries} retries.")
    return False


def curl_exists_file(link):
    # using curl test the existence of web file "link"

    command = f"curl -I -f {link} 2>/dev/null| head -n1"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    s = str(out).split(" ")
    if len(s) > 2:
        return s[1]
    else:
        return False


def convert_svg_to_png(old_svg, directory):
    # mogrify leave background transparent during the transformation
    # convert convert it to white
#    command = f"convert {old_svg} {options} -trim {new_png}"
    command = f"mogrify -background none -format png {old_svg}"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        out, err = p.communicate()
        log_fn = f"{directory}/convert_svg_to_png-fail-{time.time()}.txt"
        error_report(directory, command, out, err, log_fn)
        logging.error(f"convert on {old_svg} failed, report save to {log_fn}")
    return r


def unzip_archive(archive, directory):
    # TODO0 - use some python tool for this work
    command = f"unzip {archive} -d {directory}"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = f"{directory}/unzip_-fail-{time.time()}.txt"
        error_report(directory, command, out, err, log_fn)
        logging.error(f"unzip {archive} failed, report save to {log_fn}")
    return r


def gunzip_archive(archive, directory):
    command = f"gunzip {archive}"
    p = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    r = p.wait()
    if (r != 0):
        # not succesfull again, generate error log
        log_fn = f"{directory}/gunzip-fail-{time.time()}.txt"
        error_report(directory, command, out, err, log_fn)
        logging.error(f"gunzip {archive} failed, report save to {log_fn}")
    return r
