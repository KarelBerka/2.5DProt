#!/bin/sh

# This is meant to install/update helper tools used by 2DProt to the checkout

set -x
set -e

PREFFERED_PYTHON=python3.10

update_repos='false'
if echo "$@" | grep --quiet -- "--update"; then
    echo "INFO: Updating"
    update_repos='true'
fi

recreate_venvs='false'
if echo "$@" | grep --quiet -- "--recreate-venvs"; then
    echo "INFO: Recreating venvs"
    recreate_venvs='true'
fi

function do_overprot() {
    [ -d overprot/ ] || git clone https://gitlab.com/midlik/overprot.git
    cd overprot/

    # Update git repo if '--update' was passed
    if $update_repos; then
        git reset --hard
        git clean -d -f
        git pull
    fi

    # Delete venvs so they can be recreated when '--recreate-venvs' was passed
    cd OverProtCore/
    if $recreate_venvs; then
        rm -rf venv
    fi

    # Initiate Python 3 virtualenv
    [ -d venv ] || $PREFFERED_PYTHON -m venv venv/
    source venv/bin/activate
    $PREFFERED_PYTHON -m pip install -U pip
    $PREFFERED_PYTHON  -m pip install -U -r requirements.txt

    # Make sure we have latest numpy
    # Workaround for: RuntimeError: module compiled against API version 0xe but this version of numpy is 0xd
    ###pip install -U numpy

    # An experiment to install pymol to venv - prerequisities are:
    #   dnf install gcc gcc-c++ kernel-devel python3-devel glew-devel msgpack-devel freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel mmtf-cpp-devel netcdf-cxx-devel
    #   rm -rf venv/bin/pymol venv/lib64/python3.9/site-packages/{pymol,pymol2,chempy}
    rpm -q gcc gcc-c++ kernel-devel $PREFFERED_PYTHON-devel glew-devel msgpack-devel freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel mmtf-cpp-devel netcdf-cxx-devel
    [ -d pymol-open-source/ ] || git clone https://github.com/schrodinger/pymol-open-source.git
    cd pymol-open-source/
    if $update_repos; then
        git reset --hard
        git clean -d -f
        git pull
    fi
    python setup.py install --install-lib=$( pwd )/../venv/lib64/$PREFFERED_PYTHON/site-packages/ --install-scripts=$( pwd )/../venv/bin/
    pymol --version
    cd ..

    # # Workaround to add PyMol to our venv
    # DIR=$( echo $VIRTUAL_ENV/lib64/python3.*/site-packages )
    # [ -e "$DIR/pymol" ] || ln -s /usr/lib64/python3.*/site-packages/pymol "$DIR/pymol"
    # [ -e "$DIR/chempy" ] || ln -s /usr/lib64/python3.*/site-packages/chempy "$DIR/chempy"
    # [ -e "$DIR/pymol2" ] || ln -s /usr/lib64/python3.*/site-packages/pymol2 "$DIR/pymol2"

    deactivate

    cd ../..
}

function do_SecStrAnnot2() {
    # Check we have dssp installed as it is used in the config
    if ! type mkdssp; then
        echo "ERROR: mkdssp not in the path - please install dssp package" >&2
        exit 1
    fi

    # Clone SecStrAnnot2 git if it is not present
    mkdir -p overprot/OverProt
    cd overprot/OverProt
    if ! [ -d SecStrAnnot2 ]; then
        git clone https://gitlab.com/midlik/SecStrAnnot2.git
        ###cd SecStrAnnot2
        ###git reset --hard 384958767287accaf8d255a148b2ae2d2953c524   # WORKAROUND: version known to build with dotnet-sdk-3.1
        ###cd ..
    fi
    cd SecStrAnnot2/

    # Update git repo if '--update' was passed
    if $update_repos; then
        git reset --hard
        git clean -d -f
        git pull || return 1
    fi

    # Build .NET project to adjust file paths
    dotnet build

    ## Correct paths of the scripts
    #cat ../SecStrAnnotator_config.json > SecStrAnnotator_config.json

    # Recrete links
    rm -f SecStrAnnotator.dll SecStrAnnotator.runtimeconfig.json
    ln -s bin/Release/net6.0/SecStrAnnotator.dll SecStrAnnotator.dll
    ln -s bin/Release/net6.0/SecStrAnnotator.runtimeconfig.json SecStrAnnotator.runtimeconfig.json

    cd ../../..
}

function do_show2d_layout_venv() {
    # Delete venvs so they can be recreated when '--recreate-venvs' was passed
    if $recreate_venvs; then
        rm -rf venv
    fi

    # In Fedora there is numpy 1.16, but we need 1.17 which have some important bug fixed and we were hitting it when running at scale
    if [ ! -d venv/ ]; then
        $PREFFERED_PYTHON -m venv venv --system-site-packages
        source venv/bin/activate
        $PREFFERED_PYTHON -m pip install -U pip
        $PREFFERED_PYTHON -m pip install numpy
        $PREFFERED_PYTHON -m pip freeze | grep numpy
        deactivate
    fi
}

# Workaround for https://github.com/pypa/pip/issues/7883
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring

do_overprot
do_SecStrAnnot2
do_show2d_layout_venv

echo "All good"
