#!/bin/sh

set -ex

DOWNSTREAM_USER="varekova"
DOWNSTREAM_PATH="/home/$DOWNSTREAM_USER/radka/novej/tycinky/"
DOWNSTREAM_HOST="192.168.1.51"

skip_this="$1"
tarball="$( date -Iseconds | sed 's/[:+]/_/g' )-generated.tar"

if [ -n "$skip_this" ]; then
    tar cf $tarball --remove-files $( ls | grep '^generated-' | grep -v "$skip_this" )
else
    tar cf $tarball --remove-files $( ls | grep '^generated-' )
fi

scp -C $tarball $DOWNSTREAM_USER@$DOWNSTREAM_HOST:$DOWNSTREAM_PATH && rm -f $tarball
