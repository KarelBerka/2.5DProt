#!/bin/bash

set -e
set -x
set -u

# Check dependencies
rpm -q gcc cjson-devel

# Build
cd gcc
gcc -lm -o sses_main_loop -lcjson utils_matrix_sses_iteration.c
gcc -lm -o ligands_main_loop -lcjson utils_matrix_ligands_iteration.c