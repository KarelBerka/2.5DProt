#!/usr/bin/bash

SOURCE="proteins_new/"
TARGET="web_proteins_new/"
TARGET_IMAGES="images/"
PROTEIN_REGEXP='[0-9a-z][0-9a-z][0-9a-z][0-9a-z]'
ERRORS_COUNTER=0
UPDATED_COUNTER=0

for d in $( find $SOURCE -maxdepth 3 -type d -name 'generated-*' ); do
    echo "DEBUG: Processing '$d'"
    protein=$( echo "$d" | sed "s/^.*\/generated-\($PROTEIN_REGEXP\)-\?.*$/\1/" )
    version=$( echo "$d" | sed "s/^.*\/generated-$PROTEIN_REGEXP-\?\([0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T[0-9_]\+\).*$/\1/" )

    if ! echo "$protein" | grep --quiet "^$PROTEIN_REGEXP$"; then
        echo "ERROR: Failed to parse protein from '$d', skipping" >&2
	let ERRORS_COUNTER+=1
        continue
    fi

    if [ -z "$version" ]; then
        echo "ERROR: Failed to parse version from '$d', skipping" >&2
	let ERRORS_COUNTER+=1
        continue
    fi

    protein_dir="$TARGET/$TARGET_IMAGES/generated-$protein"
    layout_dir="$TARGET/$TARGET_LAYOUT/generated-$protein"
    version_dir="$protein_dir/$version"

    if [ -d "$version_dir" ]; then
        # just for once copy prefamilies.log here
        echo "ERROR: Target directory '$version_dir' already exists, skipping" >&2
	let ERRORS_COUNTER+=1
        continue
    else
        echo "$TARGET_IMAGES"
        echo "$TARGET/$TARGET_IMAGES"
        rm -rf  $TARGET/$TARGET_IMAGES/generated-$protein/*
        mkdir -p "$version_dir"
        echo "copy to " $version_dir/${protein}[_c,_l,].svg

        cp -la $d/${protein}.svg $version_dir/${protein}.svg
#        cp -la $d/${protein}_c.svg $version_dir/${protein}_c.svg
#        cp -la $d/${protein}_l2.svg $version_dir/${protein}_l.svg
    fi
done

echo "ERRORS_COUNTER = $ERRORS_COUNTER"
echo "UPDATED_COUNTER = $UPDATED_COUNTER"
