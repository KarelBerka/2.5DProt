#!/usr/bin/bash

SOURCE="processed/"
TARGET="web/"
TARGET_IMAGES="images/"
TARGET_LAYOUT="layouts/"
FAMILY_REGEXP='[0-9]\+\.[0-9]\+\.[0-9]\+\.[0-9]\+'
ERRORS_COUNTER=0
COUNTER=0

for d in $( find $SOURCE -maxdepth 1 -type d -name 'generated-*' ); do

    # set family and version variable
    family=$( echo "$d" | sed "s/^.*\/generated-\($FAMILY_REGEXP\)-\?.*$/\1/" )
    version=$( echo "$d" | sed "s/^.*\/generated-$FAMILY_REGEXP-\?\([0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T[0-9_]\+\).*$/\1/" )

    if ! echo "$family" | grep --quiet "^$FAMILY_REGEXP$"; then
        echo "ERROR: Failed to parse family from '$d', skipping" >&2
	let ERRORS_COUNTER+=1
        continue
    fi

    if [ -z "$version" ]; then
        echo "ERROR: Failed to parse version from '$d', skipping" >&2
	let ERRORS_COUNTER+=1
        continue
    fi

    # create main directories
    family_dir="$TARGET/$TARGET_IMAGES/generated-$family"
    layout_dir="$TARGET/$TARGET_LAYOUT/generated-$family"
    version_dir="$family_dir/$version"

    if [ -d "$version_dir" ]; then
        # just for once copy prefamilies.log here
        echo "ERROR: Target directory '$version_dir' already exists, skipping" >&2
        # cp $d/prefamilies.log $version_dir
	let ERRORS_COUNTER+=1
        continue
    else
        rm -rf  $family_dir/*
        mkdir -p "$version_dir"

        # copy all multiple family diagrams
        cp $d/$family.svg $version_dir/imageCircle-multiple.svg
        cp $d/${family}_2.svg $version_dir/imageCircle-multiple2.svg
        cp $d/${family}_3.svg $version_dir/imageCircle-multiple3.svg
        cp $d/${family}_4.svg $version_dir/imageCircle-multiple4.svg
        cp $d/$family.png $version_dir/imageCircle-multiple.png
        cp $d/${family}_2.png $version_dir/imageCircle-multiple2.png
        cp $d/${family}_3.png $version_dir/imageCircle-multiple3.png
        cp $d/${family}_4.png $version_dir/imageCircle-multiple4.png

        # copy all multiple clusters diagrams
        for dd in $( find $d -maxdepth 1 -type d -name $family* ); do
            cluster=$( basename $dd )
            mkdir -p "$version_dir/$cluster"
            cp $dd/$family.svg $version_dir/$cluster/imageCircle-multiple.svg
            cp $dd/${family}_2.svg $version_dir/$cluster/imageCircle-multiple2.svg
            cp $dd/${family}_3.svg $version_dir/$cluster/imageCircle-multiple3.svg
            cp $dd/${family}_4.svg $version_dir/$cluster/imageCircle-multiple4.svg
            cp $dd/$family.png $version_dir/$cluster/imageCircle-multiple.png
            cp $dd/${family}_2.png $version_dir/$cluster/imageCircle-multiple2.png
            cp $dd/${family}_3.png $version_dir/$cluster/imageCircle-multiple3.png
            cp $dd/${family}_4.png $version_dir/$cluster/imageCircle-multiple4.png
            cp $dd/domain_list.txt $version_dir/$cluster
            cp $dd/ligands.json $version_dir/$cluster
        done

        # copy single/single layout svgs and json data
        for i in `seq 0 9`; do
           for j in `ls  $d/$i*[A-Z]*l.svg 2>/dev/null`; do
               lsvg_file=$( basename $j )
               svg_file=`echo $lsvg_file | sed "s|l.svg|.svg|"`
               json_file=`echo $lsvg_file | sed "s|l.svg|.json|"`
               cp $d/$lsvg_file $version_dir/image-$lsvg_file
               cp $d/$svg_file $version_dir/image-$svg_file
               cp $d/$json_file $version_dir/image-$json_file
           done
        done

        # copy layouts
        find $d -maxdepth 1 -type f -name 'layout-*.json' -exec cp -t $version_dir '{}' +
        # cp $d/prefamilies.log $version_dir TODO
        cp $d/start-layout.json $version_dir
        cp $d/rotation.json $version_dir
        cp $d/${family}.json $version_dir
        cp $d/domain_list.txt $version_dir   # necessary to check wwhether there were no changes in family domain list
        cp $d/ligands.json $version_dir
        let OK+=1

#        for dd in $( find $d -maxdepth 1 -type d -name $family* -or -name None ); do
#            echo "DEBUG: Processing $dd"
#            subfamily=$( basename $dd )
#	    subfamily_dir="$version_dir/$subfamily"
#
#            if [ -d "$subfamily_dir" ]; then
#    	        echo "ERROR: Target subfamily directory '$subfamily_dir' already exists, skipping" >&2
#	        let ERRORS_COUNTER+=1
#                continue
#            fi
#
#	    mkdir -p "$subfamily_dir"
##	    find $dd -maxdepth 1 -type f -name 'imageCircle-multiple*.svg' -exec cp -t $subfamily_dir '{}' +
##	    find $dd -maxdepth 1 -type f -name 'imageCircle-multiple*.png' -exec cp -t $subfamily_dir '{}' +
#	    for svg in $subfamily_dir/imageCircle-multiple*.svg; do
#                png="$( echo $svg | sed 's/\.svg$/.png/' )"
#
#                if [ ! -f "$png" ]; then
##                    echo "DEBUG: Converting '$svg' to '$png'" >&2
#                    convert $svg -trim $png
#                fi
#	    done
#        done
#        echo "INFO: Copied '$d' to '$version_dir'" >&2
    fi

    # copy layouts as well 
    rm -rf "$layout_dir"
    mkdir -p "$layout_dir"
##    cp $d/layout*.json  $layout_dir
##    cp $d/start-layout.json  $layout_dir
##    cp $d/rotation.json  $layout_dir
##    cp $d/$family.json  $layout_dir
##    cp $d/domain_list.txt $layout_dir

    echo "generated-$family-$version" >$layout_dir/version.txt
#    echo "INFO: Copied '$d' to '$layout_dir'" >&2

done

echo "ERRORS_COUNTER = $ERRORS_COUNTER"
echo "OK" = $OK