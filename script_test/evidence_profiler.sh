#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH

mkdir test
input=./input_data/evidence_profiler
output=./returned/evidence_profiler

mkdir -p $output


evidence_profiler.py -p $input/raw_data/profile  -E ~pedro/proyectos/pets/source_data/project/ext_data/processed/monarch -o $output/evidence_reports -V $input/raw_data/variants
#bundle exec ruby bin/evidence_profiler.rb -p $current/raw_data/profile  -E $current/sample_evid -o $current/evidence_reports
#bundle exec ruby bin/evidence_profiler.rb -p $current/raw_data/profile  -e $current/processed_data/gene_phenotype -o $current/evidence_reports -g $current/processed_data/gene_coordinates_clean
cd $current
