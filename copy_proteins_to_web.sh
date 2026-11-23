#!/usr/bin/bash

SOURCE="proteins/"
TARGET="web_proteins/"
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
        rm -rf  $TARGET/$TARGET_IMAGES/generated-$protein/*
        mkdir -p "$version_dir"
        for i  in `ls  $d/*.svg 2>/dev/null`; do
            svg_file=$( basename $i )
            c=`echo $svg_file | head -c 5 | tail -c 1`
            echo $svg_file $c
            if [ "$c" = "l" ]; then
                p=`echo $svg_file | head -c 4`
                cp $d/$svg_file $version_dir/${p}_l.svg
            else
                cp $d/$svg_file $version_dir/$svg_file
            fi
        done

        find $d -maxdepth 1 -type f -name 'image-*.svg' -exec cp -t $version_dir '{}' +
	let UPDATED_COUNTER+=1
    fi
done

echo "ERRORS_COUNTER = $ERRORS_COUNTER"
echo "UPDATED_COUNTER = $UPDATED_COUNTER"
