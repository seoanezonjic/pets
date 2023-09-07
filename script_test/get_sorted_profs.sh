#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH

mkdir -p ./returned/get_sorted_profs ./expected/get_sorted_profs ./input_data/get_sorted_profs
#cohorts_path="/mnt/home/users/bio_267_uma/elenarojano/projects/pets/CohortAnalyzer/cohorts/paper_cohorts"
#profile_path="/mnt/home/users/pab_001_uma/pedro/proyectos/pets/pets_reloaded_testing/get_sorted_profs"
cohorts_path="./input_data/get_sorted_profs"
profile_path="./input_data/get_sorted_profs"
output_path="./returned/get_sorted_profs"

get_sorted_profs.py -c 'chr' -s 'start' -e 'stop' -d 'patient_id' -p 'phenotypes' -S ',' \
                    -P $cohorts_path'/pmm2_paco_format.txt' \
                    -r $profile_path/profile \
                    -o $output_path/report.html
