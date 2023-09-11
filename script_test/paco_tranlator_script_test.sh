#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH
data_to_test=./input_data/paco_translator

mkdir -p ./returned/paco_translator

paco_translator.py -P $data_to_test/cohort_toy_dataset.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m paco -t -n -C --n_phens 2 \
    -o ./returned/paco_translator/translated_paco_file.txt

paco_translator.py -P $data_to_test/cohort_toy_dataset.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m paco -n -C --n_phens 2 \
    -o ./returned/paco_translator/untranslated_paco_file.txt

paco_translator.py -P $data_to_test/cohort_toy_dataset.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m default -t -n -C --n_phens 2 \
    -o ./returned/paco_translator/translated_default_file.txt

paco_translator.py -P $data_to_test/cohort_toy_dataset.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m default -n -C --n_phens 2 \
    -o ./returned/paco_translator/untranslated_default_file.txt


for file_to_test in `ls ./expected/paco_translator`; do
 	echo $file_to_test
 	diff ./returned/paco_translator/$file_to_test ./expected/paco_translator/$file_to_test
    rm ./returned/paco_translator/$file_to_test
done