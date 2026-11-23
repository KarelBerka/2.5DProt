#!/bin/bash

set -e
set -x
set -u


PREFFERED_DOTNET=net10.0


# Check we have dssp installed as it is used in the config
if ! type mkdssp; then
    echo "ERROR: mkdssp not in the path - please install dssp package" >&2
    exit 1
fi

# Clone SecStrAnnot2 git if it is not present
mkdir -p overprot/OverProt
cd overprot/OverProt
if ! [ -d SecStrAnnot2 ]; then
    git clone --depth=1 https://gitlab.com/midlik/SecStrAnnot2.git
fi
cd SecStrAnnot2/

# Patch target platform
sed -i "s/net6.0/$PREFFERED_DOTNET/" SecStrAnnotator.csproj
cat SecStrAnnotator.csproj

# Build .NET project to adjust file paths
rm -rf bin/
dotnet build

## Correct paths of the scripts
#cat ../SecStrAnnotator_config.json > SecStrAnnotator_config.json

# Recrete links
rm -f SecStrAnnotator.dll SecStrAnnotator.runtimeconfig.json
find bin/Debug/
mv bin/Debug/$PREFFERED_DOTNET/SecStrAnnotator.dll .
mv bin/Debug/$PREFFERED_DOTNET/SecStrAnnotator.runtimeconfig.json .

cd ../../..
