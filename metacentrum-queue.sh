#!/bin/sh

set -ex

JOBS_TIMESTAMP=2022-08-22

FRONTEND=skirit.ics.muni.cz
SSH_PARAMS="-i ~/.ssh/libra_id_rsa"

ls -al $JOBS_TIMESTAMP-meta*.txt

ssh $SSH_PARAMS ivarekova@$FRONTEND bash -c "kinit"
ssh $SSH_PARAMS ivarekova@$FRONTEND bash -c ":>jobs.list"

FAMILIES_LIST=$JOBS_TIMESTAMP-meta2.txt ./metacentrum-jobs.sh 100
echo "for f in bilkoviny-*.sh; do qsub -l select=1:ncpus=4:mem=8gb:scratch_local=15gb:cgroups=cpuacct -l walltime=24:00:00 \$f; sleep .1; done >>jobs.list" \
    | ssh $SSH_PARAMS ivarekova@$FRONTEND bash -

FAMILIES_LIST=$JOBS_TIMESTAMP-meta3.txt ./metacentrum-jobs.sh 10
echo "for f in bilkoviny-*.sh; do qsub -l select=1:ncpus=4:mem=8gb:scratch_local=15gb:cgroups=cpuacct -l walltime=24:00:00 \$f; sleep .1; done >>jobs.list" \
    | ssh $SSH_PARAMS ivarekova@$FRONTEND bash -

FAMILIES_LIST=$JOBS_TIMESTAMP-meta4.txt ./metacentrum-jobs.sh 2
echo "for f in bilkoviny-*.sh; do qsub -l select=1:ncpus=6:mem=12gb:scratch_local=15gb:cgroups=cpuacct -l walltime=24:00:00 \$f; sleep .1; done >>jobs.list" \
    | ssh $SSH_PARAMS ivarekova@$FRONTEND bash -

FAMILIES_LIST=$JOBS_TIMESTAMP-meta5.txt ./metacentrum-jobs.sh 1
echo "for f in bilkoviny-*.sh; do qsub -l select=1:ncpus=8:mem=16gb:scratch_local=15gb:cgroups=cpuacct -l walltime=48:00:00 \$f; sleep .1; done >>jobs.list" \
    | ssh $SSH_PARAMS ivarekova@$FRONTEND bash -
