#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH

input=./input_data/evidence_profiler
output=./returned/evidence_profiler
mkdir -p $output

evidence_profiler.py -p $input/raw_data/profile  -E ~pedro/proyectos/pets/source_data/project/ext_data/processed/monarch -o $output/evidence_reports -V $input/raw_data/variants
