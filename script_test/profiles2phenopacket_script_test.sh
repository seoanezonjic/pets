#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH
data_to_test=./input_data/paco_translator


profiles2phenopacket.py -i $data_to_test/cohort_toy_dataset.txt \
    -c chr -d patient_id -s start -e end -p phenotypes -x "sex" -S "|" -n \
    -o ./returned/profiles2phenopacket


ls ./expected/profiles2phenopacket | sort > ./expected_phenopacket_files.txt
ls ./returned/profiles2phenopacket | sort > ./returned_phenopacket_files.txt
diff ./expected_phenopacket_files.txt ./returned_phenopacket_files.txt
rm ./expected_phenopacket_files.txt ./returned_phenopacket_files.txt

mkdir -p ./returned/profiles2phenopacket

for file_to_test in `ls ./expected/profiles2phenopacket`; do
 	echo $file_to_test
    wc ./returned/profiles2phenopacket/$file_to_test | cut -f 2,4,5 -d " " > ./returned_numbers.txt
    wc ./expected/profiles2phenopacket/$file_to_test | cut -f 2,4,5 -d " " > ./expected_numbers.txt
    diff ./returned_numbers.txt ./expected_numbers.txt
    rm ./returned_numbers.txt ./expected_numbers.txt

    sort ./returned/profiles2phenopacket/$file_to_test > ./tmp_returned_sorted.txt
    sort ./expected/profiles2phenopacket/$file_to_test > ./tmp_expected_sorted.txt
    diff ./tmp_returned_sorted.txt ./tmp_expected_sorted.txt
    rm ./tmp_returned_sorted.txt ./tmp_expected_sorted.txt
 	
    rm ./returned/profiles2phenopacket/$file_to_test
    #diff ./returned/profiles2phenopacket/$file_to_test ./expected/profiles2phenopacket/$file_to_test
done