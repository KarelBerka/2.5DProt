#!/bin/bash

set -xeuo pipefail

# Check prerequisite packages are installed
rpm -q mono-core

# Download Mole2
mkdir -p Mole2_cmd
cd Mole2_cmd
curl -o Mole2_cmd.zip https://webchem.ncbr.muni.cz/Platform/AppsBin/Mole/2.5.24.6.8/Mole2_cmd.zip
unzip -o -d . Mole2_cmd.zip
rm Mole2_cmd.zip
rm -rf tests/
mv Mole2_cmd/* .
rmdir Mole2_cmd

# Test it works
mono mole2.exe --help
