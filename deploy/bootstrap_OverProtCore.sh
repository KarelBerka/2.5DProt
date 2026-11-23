#!/bin/bash

set -e
set -x
set -u


PREFFERED_PYTHON=python3.10
PREFFERED_DOTNET=net10.0
PREFFERED_DOTNET_VERSION="10.0.0"   # Minor version floor for roll-forward (dotnet resolves 10.0.0 to installed 10.0.x)
PREFFERED_DOTNET_STRING="Microsoft.NETCore.App 10."


[ -d overprot/.git ] || git clone --depth 1 https://gitlab.com/midlik/overprot.git
cd overprot/OverProtCore/

# Initiate Python 3 virtualenv
[ -d venv ] || $PREFFERED_PYTHON -m venv venv/
source venv/bin/activate
$PREFFERED_PYTHON -m pip install -U pip
$PREFFERED_PYTHON -m pip install -U -r requirements.txt

# Because pymol bellow bumps numpy, we need additional dependencies updated as well
$PREFFERED_PYTHON -m pip install numpy==2.2.6 matplotlib==3.10.9 contourpy==1.3.2 llvmlite==0.49.0 numba==0.67.0 scipy==1.15.3

# Install pymol to venv
#   dnf install gcc gcc-c++ kernel-devel python3-devel glew-devel msgpack-devel freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel mmtf-cpp-devel netcdf-cxx-devel
#   rm -rf venv/bin/pymol venv/lib64/python3.9/site-packages/{pymol,pymol2,chempy}
rpm -q gcc gcc-c++ kernel-devel $PREFFERED_PYTHON-devel glew-devel msgpack-devel freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel mmtf-cpp-devel netcdf-cxx-devel
[ -d pymol-open-source/.git ] || git clone --depth 1 https://github.com/schrodinger/pymol-open-source.git
cd pymol-open-source/
###python setup.py install --install-lib=$( pwd )/../venv/lib64/$PREFFERED_PYTHON/site-packages/ --install-scripts=$( pwd )/../venv/bin/
pip install .
pymol --version
cd ..

deactivate

# Enforce newer DotNET
sed -i "s/net6\.0/$PREFFERED_DOTNET/" overprot/libs/lib_dependencies.py
sed -i "s/Microsoft.NETCore.App 6\.0/$PREFFERED_DOTNET_STRING/" overprot/libs/lib_dependencies.py

# Enforce newer DotNET for prebuilt SecStrAnnotator (shipped in overprot repo,
# not rebuildable here). Its runtimeconfig pins Microsoft.NETCore.App 6.0.0,
# but without 6.x runtime it exits with 150 (framework not found). Roll forward
# to installed 10.x; net6.0 IL is compatible with the 10.x runtime.
sed -i "s/net6\.0/$PREFFERED_DOTNET/" dependencies/SecStrAnnotator/SecStrAnnotator.runtimeconfig.json
sed -i "s/\"version\": \"6\.0\.0\"/\"version\": \"$PREFFERED_DOTNET_VERSION\"/" dependencies/SecStrAnnotator/SecStrAnnotator.runtimeconfig.json
cat dependencies/SecStrAnnotator/SecStrAnnotator.runtimeconfig.json

# Fail fast if dotnet host still cannot load the dll (150 = runtime not found)
dotnet dependencies/SecStrAnnotator/SecStrAnnotator.dll --help

cd dependencies/

# Prepare included StructureCutter
# Patch target platform
cd StructureCutter/
sed -i "s/net6.0/$PREFFERED_DOTNET/" StructureCutter.csproj
cat StructureCutter.csproj

# Build .NET project to adjust file paths
rm -rf bin/
dotnet build -c Release
cd ..

cd ../../..
