#!/bin/bash

CONTAINER=bilkoviny-20190805.simg
REPOSITORY=bilkoviny.tar.gz
###DOMAINS=50

MYHOME="/storage/brno2/home/ivarekova/"

: <<'DOCUMENTATION'
# Build and upload container
sudo mount -o remount,size=12G,noatime /tmp
sudo singularity build bilkoviny-20190805.simg Singularity
sudo singularity build bilkoviny-20190805.img bilkoviny-20190805.simg
scp bilkoviny-*.simg ivarekova@skirit.ics.muni.cz:/storage/brno2/home/ivarekova

# Build and upload code
./metacentrum-deploy.sh

# Prepare and upload job script files (also shows how to run):
./metacentrum-jobs.sh
# Examples:
#   ./metacentrum-jobs.sh 50   # first 50 families from randomized list
#   ./metacentrum-jobs.sh 50 2   # second 50 families from same randomized list
#   ./metacentrum-jobs.sh 100 1 5   # first 100 families in 20 job scripts (each script processes 5 families)

# Download results
./metacentrum-download.sh
DOCUMENTATION

# Setup SCRATCH cleaning in case of an error
trap 'clean_scratch' TERM EXIT

function img_bash() {
    singularity exec --bind bilkoviny/:/bilkoviny $CONTAINER bash
}

function doit_family() {
    family=$1
    echo "Working on family: $family"
    echo -e "
        cd /bilkoviny
        export HOME=/bilkoviny
        export OPENBLAS_NUM_THREADS=4
        ./metacentrum-versions.sh
        cd generate_family_entry/
        time ./generate_family_dcs_diagrams.py --family $family
        " | img_bash
    rc=$?
    mv bilkoviny/generate_family_entry/generated-$family-* bilkoviny/
    ###mv "$( find bilkoviny/generated-$family-* -name sses-matches.svg )" bilkoviny/generated-$family-* 2>/dev/null
    ###mv "$( find bilkoviny/generated-$family-* -name consensus-template.sses.json | head -n 1 )" bilkoviny/generated-$family-* 2>/dev/null
    ###mv "$( find bilkoviny/generated-$family-* -name consensus.cif | head -n 1 )" bilkoviny/generated-$family-* 2>/dev/null
    ###mv "$( find bilkoviny/generated-$family-* -name sses-stats.json )" bilkoviny/generated-$family-* 2>/dev/null
    ###mv "$( find bilkoviny/generated-$family-* -name stats.json )" bilkoviny/generated-$family-* 2>/dev/null
    ###rm -rf bilkoviny/generated-$family-*/$family/
    ###find bilkoviny/generated-$family-* -name \*.cif ! -name consensus.cif -delete
    cp bilkoviny/version-*.log bilkoviny/generated-$family-*/
    mv bilkoviny/generated-$family-* $MYHOME
    if [ $rc -ne 0 ]; then
        echo "ERROR: Return code for $family: $rc"
    else
        echo "Return code for $family: $rc"
    fi
}

function doit_protein() {
    protein=$1
    data_dir="generated-$protein-$( date --utc -Ins | sed 's/[^a-zA-Z0-9-]/_/g' )"
    echo "Working on protein: $protein"
    echo -e "
        cd /bilkoviny
        export HOME=/bilkoviny
        export OPENBLAS_NUM_THREADS=4
        ./metacentrum-versions.sh
        cd generate_protein_entry
        time ./generate_protein_diagram.py --protein $protein --output $data_dir
        " | img_bash
    rc=$?
    mv bilkoviny/generate_protein_entry/$data_dir/ $MYHOME
    if [ $rc -ne 0 ]; then
        echo "ERROR: Return code for $family: $rc"
    else
        echo "Return code for $family: $rc"
    fi
}

cd $SCRATCHDIR
rm -rf bilkoviny/ bilkoviny*.tar.gz

echo "Pwd: $( pwd )"
echo "Uname: $( uname -a )"
echo "Du here: $( du -hs . )"

# Prepare dir with the project
if [ ! -d bilkoviny/ ]; then
    cp $MYHOME/$REPOSITORY .
    tar xzf $REPOSITORY
    ###sed -i "s/^DOMAINS=.*/DOMAINS=$DOMAINS/" bilkoviny/generate-new.sh
fi
rm -rf bilkoviny/generated-*
rm -rf bilkoviny/input/sses-stats/*
rm -rf $( find bilkoviny/input/protein-groups/* -type f -name '*.json' | grep -v '^bilkoviny/input/protein-groups/default.json$' )
rm -rf bilkoviny/input/annotations/*
rm -rf bilkoviny/input/cytos-cache-detected/*
rm -rf bilkoviny/input/pdbs/*
:>     bilkoviny/input/tagging.yaml

rm -f bilkoviny/SecStrAnnot2/bin/Debug/netcoreapp2.0/SecStrAnnotator_config.json
pushd bilkoviny/SecStrAnnot2/bin/Debug/netcoreapp2.0/
ln -s ../../../../SecStrAnnotator_config.json SecStrAnnotator_config.json
popd

# Prepare and test container image
if [ ! -f $CONTAINER ]; then
    cp $MYHOME/$CONTAINER .
fi
echo 'date' | img_bash
rc=$?
if [ $rc -ne 0 ]; then
    echo "ERROR: Unable to use singularity"
    exit 1
fi


