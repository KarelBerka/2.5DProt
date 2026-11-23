#!/bin/sh

set -xe

( date; uname --nodename; env ) >version-run.log

rpm -qa | sort > version-rpm-qa.log

git log | head -n 100 > version-git-log-bilkoviny.log

source venv/bin/activate; pip freeze > ../version-pip-freeze-show2d_layout_venv.log; deactivate

cd ubertemplate/; git log | head -n 100 > ../version-git-log-ubertemplate.log; source venv/bin/activate; pip freeze > ../version-pip-freeze-ubertemplate.log; deactivate; cd ..

cd SecStrAnnot2/; git log | head -n 100 > ../version-git-log-SecStrAnnot2.log; cd ..
