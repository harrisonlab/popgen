#!/bin/bash

input=/home/sobczm/popgen/input/mappings

### Rename input files in each folder by prefixing with the strain ID
cd $input/125/
for filename in *; do mv "$filename" "125_$filename"; done;
cd $input/55/
for filename in *; do mv "$filename" "55_$filename"; done;
cd $input/A1-2/
for filename in *; do mv "$filename" "A1-2_$filename"; done;
cd $input/A13/
for filename in *; do mv "$filename" "A13_$filename"; done;
cd $input/A23/
for filename in *; do mv "$filename" "A23_$filename"; done;
cd $input/A28/
for filename in *; do mv "$filename" "A28_$filename"; done;
cd $input/CB3/
for filename in *; do mv "$filename" "CB3_$filename"; done;
cd $input/D2/
for filename in *; do mv "$filename" "D2_$filename"; done;
cd $input/Fus2/
for filename in *; do mv "$filename" "Fus2_$filename"; done;
cd $input/HB6/
for filename in *; do mv "$filename" "HB6_$filename"; done;
cd $input/PG/
for filename in *; do mv "$filename" "PG_$filename"; done;

### Add read group and sample name to each mapped read
