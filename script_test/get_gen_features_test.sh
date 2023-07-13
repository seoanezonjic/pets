#!/usr/bin/env bash
source ~soft_bio_267/initializes/init_python
export PATH=../bin/:$PATH
single_case_test=./input_data/get_gene_features/one_patient_one_feature.txt
data_to_test=../test/data/cohort_toy_dataset.txt
reference_file=../external_data/gencode.v43.basic.annotation.gtf.gz
mkdir -p ./returned/get_gene_features

######## Testing just one patient with one genomic feature

get_gen_features.py -i $single_case_test \
    -c chr -d patient_id -s start -e end  \
    -o ./returned/get_gene_features/single_case_ensemble.txt \
    -r $reference_file -t gene

echo Single patient with one feature
echo single_case_ensemble.txt
diff ./returned/get_gene_features/single_case_ensemble.txt \
    ./expected/get_gen_features/single_case_ensemble.txt
rm ./returned/get_gene_features/single_case_ensemble.txt

get_gen_features.py -i $single_case_test \
    -c chr -d patient_id -s start -e end  \
    -o ./returned/get_gene_features/single_case_geneName.txt \
    -r $reference_file -t gene -n gene_name

echo single_case_geneName.txt
diff ./returned/get_gene_features/single_case_geneName.txt \
    ./expected/get_gen_features/single_case_geneName.txt
rm ./returned/get_gene_features/single_case_geneName.txt

####### Testing 4 patients (3 of them with several genomic features)

get_gen_features.py -i $data_to_test \
    -c chr -d patient_id -s start -e end  \
    -o ./returned/get_gene_features/patients_ensemble.txt \
    -r $reference_file -t gene

echo Multiple patients
echo patients_ensemble.txt
diff ./returned/get_gene_features/patients_ensemble.txt \
    ./expected/get_gen_features/patients_ensemble.txt
rm ./returned/get_gene_features/patients_ensemble.txt

get_gen_features.py -i $data_to_test \
    -c chr -d patient_id -s start -e end  \
    -o ./returned/get_gene_features/patients_geneName.txt \
    -r $reference_file -t gene -n gene_name

echo patients_geneName.txt
diff ./returned/get_gene_features/patients_geneName.txt \
    ./expected/get_gen_features/patients_geneName.txt
rm ./returned/get_gene_features/patients_geneName.txt