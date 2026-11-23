#!/bin/sh

set -xe

archive="$( date +%Y-%m-%d )-generated-$( date +%s ).tar.gz"

ssh ivarekova@skirit.ics.muni.cz "tar -czf $archive --remove-files bilkoviny-*.sh.o* bilkoviny-*.sh.e* generated-*"

scp ivarekova@skirit.ics.muni.cz:$archive .
