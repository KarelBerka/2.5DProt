#!/bin/sh

set -ex

FRONTEND=perian.ncbr.muni.cz
FRONTEND=skirit.metacentrum.cz

cwd=$( pwd )

rm -rf ../bilkoviny-staging
mkdir ../bilkoviny-staging
cp -r . ../bilkoviny-staging/bilkoviny

cd ../bilkoviny-staging/bilkoviny

# Make sure ubertemplate and SecStrAnnot2 are ready
./bootstrap_tools.sh $@

# Cleanup big unneeded bilkoviny/ files
rm -f *.img *.simg
rm -f *.tar.gz
rm -rf venv/
rm -rf proteins/ web/ processed/
rm -rf .git/

# Cleanup big unneeded opverprot/ files
cd overprot/
rm -rf .git/
rm -rf OverProtServer/ OverProtViewer/
cd ..

## Cleanup big unneeded SecStrAnnot2/ files
#cd SecStrAnnot2/
#rm -rf scripts/weblogo/
#rm -rf venv/
#rm -rf .git/
#cd ..

# Create tar.gz
cd ..
du -hs bilkoviny/
tar czf bilkoviny.tar.gz bilkoviny/

# Deploy to Metacentrum
rsync --compress --info=progress2 --human-readable bilkoviny.tar.gz ivarekova@$FRONTEND:bilkoviny.tar.gz
[ $? -eq 0 ] && rm -rf bilkoviny.tar.gz bilkoviny/
