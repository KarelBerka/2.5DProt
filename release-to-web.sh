#!/bin/sh

set -x -e -o pipefail

DOWNSTREAM_PATH="/home/varekova/radka/novej/tycinky/web/"
DOWNSTREAM_HOST="192.168.1.51"
UPSTREAM_PATH="/var/www/html/2DProt/"
UPSTREAM_HOST="78.128.250.85"

if [ "$USER" == "varekova" ]; then
    rsync -C --info=progress2 -rah --delete $DOWNSTREAM_PATH centos@$UPSTREAM_HOST:$UPSTREAM_PATH
elif [ "$USER" == "pok" ]; then
    echo -e "Run this:\n\ntime rsync -C --info=progress2 -rah --delete $DOWNSTREAM_PATH centos@$UPSTREAM_HOST:$UPSTREAM_PATH && exit\n\nNow logging to $DOWNSTREAM_HOST"
    ssh varekova@$DOWNSTREAM_HOST
else
    echo "ERROR: Running from unknown user/host?" >&2
    exit 1
fi

ssh centos@$UPSTREAM_HOST "cd /home/centos/2dprot-ui && export FLASK_ENV=development && export FLASK_APP=two_d_prot_ui && source venv/bin/activate && time flask init-db"

echo "Looks like it went well, please check http://$UPSTREAM_HOST:5000/"
