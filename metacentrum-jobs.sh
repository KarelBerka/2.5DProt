#!/bin/bash

set -e

FRONTEND=perian.ncbr.muni.cz
FRONTEND=skirit.ics.muni.cz

# How many families to process in each job?
load=${1:-1}

# Cleanup previously existing job scripts
rm -f bilkoviny-*.sh

function gen() {
    local file=$1
    local function=$2

    # Check file exist
    if [ ! -r "$file" ]; then
        echo "ERROR: File $file not readable" >&2
        exit 1
    fi

    # Ensure right line ends
    dos2unix $file 2>/dev/null

    local page=0
    while true; do
        from=$(( $page * $load + 1))
        to=$(( ( $page + 1 ) * $load ))
        end=$(( $to + 1 ))
        tmp=$( mktemp )
        page_long="$( printf %06.0f $page )"
        sed -n "${from},${to}p;${end}q" $file > $tmp
        if [ -s $tmp ]; then
            cp Singularity-script-template.sh bilkoviny-$page_long.sh
            for i in $( cat $tmp ); do
                echo "$function $i" >>bilkoviny-$page_long.sh
            done
            echo bilkoviny-$page_long.sh
            let page+=1
        else
            break
        fi
    done
}


# Prepare job script files
if [ -n "$FAMILIES_LIST" ]; then
    echo "# Generating job script files for families:"
    gen $FAMILIES_LIST doit_family
elif [ -n "$PROTEINS_LIST" ]; then
    echo "# Generating job script files for proteins:"
    gen $PROTEINS_LIST doit_protein
else
    echo "ERROR: Please provide FAMILIES_LIST or PROTEINS_LIST" >&2
    exit 1
fi


# Remove job script files on Metacentrum
echo "# Removing remote job script files"
ssh ivarekova@$FRONTEND "rm -f bilkoviny-*.sh"

# Upload
echo "# Uploading job script files:"
scp bilkoviny-*.sh ivarekova@$FRONTEND:/storage/brno2/home/ivarekova/

# Suggest how to run them
echo '
# To start the jobs:
for f in bilkoviny-*.sh; do
    qsub -l select=1:ncpus=4:mem=8gb:scratch_local=15gb:cgroups=cpuacct -l walltime=24:00:00 $f
    sleep .5
done

# Then to query state of the jobs:
qstat -u ivarekova
'
