#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH
data_to_test=../test/data


paco_translator.py -P $data_to_test/100_test_dataset_with_sex.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m paco -t -n \
    -o ./returned/paco_translator/translated_paco_file.txt

paco_translator.py -P $data_to_test/100_test_dataset_with_sex.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m paco -n \
    -o ./returned/paco_translator/untranslated_paco_file.txt

paco_translator.py -P $data_to_test/100_test_dataset_with_sex.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m default -t -n \
    -o ./returned/paco_translator/translated_default_file.txt

paco_translator.py -P $data_to_test/100_test_dataset_with_sex.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -S "|" -m default -n \
    -o ./returned/paco_translator/untranslated_default_file.txt


for file_to_test in `ls ./expected/paco_translator`; do
 	echo $file_to_test
 	diff `sort ./returned/paco_translator/$file_to_test` `sort ./expected/paco_translator/$file_to_test`
done