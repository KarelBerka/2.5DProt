#!/usr/bin/bash

SOURCE="processed_new/"
TARGET="web_new/"
TARGET_IMAGES="images/"
TARGET_LAYOUT="layouts/"
FAMILY_REGEXP='[0-9]\+\.[0-9]\+\.[0-9]\+\.[0-9]\+'
ERRORS_COUNTER=0
COUNTER=0

for d in $( find $SOURCE -maxdepth 1 -type d -name 'generated-*' ); do
    family=`echo $d | sed "s|[^-]*-||" | sed "s|-.*||"`
    echo "family:" $family
    # set family and version variable
    family=$( echo "$d" | sed "s/^.*\/generated-\($FAMILY_REGEXP\)-\?.*$/\1/" )
    version=$( echo "$d" | sed "s/^.*\/generated-$FAMILY_REGEXP-\?\([0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T[0-9_]\+\).*$/\1/" )

    if ! echo "$family" | grep --quiet "^$FAMILY_REGEXP$"; then
        echo "ERROR1: Failed to parse family from '$d', skipping" >&2
        let ERRORS_COUNTER+=1
        continue
    fi

    if [ -z "$version" ]; then
        echo "ERROR2: Failed to parse version from '$d', skipping" >&2
        let ERRORS_COUNTER+=1
        continue
    fi

    # create main directories
    family_dir="$TARGET/$TARGET_IMAGES/generated-$family"
    version_dir="$family_dir/$version"

    if [ -d "$version_dir" ]; then
        # just for once copy prefamilies.log here
        echo "ERROR3: Target directory '$version_dir' already exists, skipping" >&2
        # cp $d/prefamilies.log $version_dir
        let ERRORS_COUNTER+=1
        continue
    else
        rm -rf  $family_dir/*
        mkdir -p "$version_dir"

        # copy all multiple family diagrams
        cp $d/${family}_*.svg $version_dir
        cp $d/${family}_*.png $version_dir

        # copy all multiple clusters diagrams
        for dd in $( find $d -maxdepth 1 -type d -name $family* ); do
            cluster=$( basename $dd )
            mkdir -p "$version_dir/$cluster"
            cp $dd/${family}*_*.svg $version_dir/$cluster
            cp $dd/${family}*_*.png $version_dir/$cluster
            cp $dd/domain_list.txt $version_dir/$cluster
            cp $dd/ligands.json $version_dir/$cluster
        done

        # copy single/single layout svgs and json data
        for i in `seq 0 9`; do
#           echo "$d/$i[a-z]*_[A-Z]*[0-9]*.svg"
           for j in `ls  $d/$i[a-z]*_[A-Z]*[0-9]*.svg 2>/dev/null`; do
               svg_file=$( basename $j )
               json_file=`echo $svg_file | sed "s|svg|json|"`
               cp $d/$svg_file $version_dir/$svg_file
               cp $d/$json_file $version_dir/$json_file
           done
        done

        # copy layouts
        find $d -maxdepth 1 -type f -name 'layout-*.json' -exec cp -t $version_dir '{}' +
        # cp $d/prefamilies.log $version_dir TODO
        cp $d/start-layout.json $version_dir
        cp $d/rotation.json $version_dir
        cp $d/${family}.json $version_dir
        cp $d/domain_list.txt $version_dir
        cp $d/ligands.json $version_dir
        let OK+=1
    fi
done

echo "ERRORS_COUNTER = $ERRORS_COUNTER"
echo "OK" = $OK