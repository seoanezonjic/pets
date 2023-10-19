import sys, os, json, unittest 
import subprocess
from importlib.resources import files

from subprocess import PIPE
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort
from pets import get_gen_features, get_sorted_profs, paco_translator, profiles2phenopacket, cohort_analyzer, evidence_profiler, HPO_FILE, GENCODE
import warnings
import numpy as np
import pytest
from io import StringIO
import json
import shutil


ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
SCRIPT_DATA_TEST_PATH = os.path.join(ROOT_PATH, 'input_data')
EXPECTED_PATH = os.path.join(ROOT_PATH, 'expected')
RETURNED_PATH = os.path.join(ROOT_PATH, 'returned')


def strng2table(strng, fs="\t", rs="\n"):
	table = [row.split(fs) for row in strng.split(rs)][0:-1]
	return table

def test_evidence_profiler():
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'evidence_profiler')}", exist_ok=True)

    list_of_args = ["-p", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'evidence_profiler', 'raw_data', 'profile')}", 
                    "-E", f"{os.path.join('~pedro', 'proyectos', 'pets', 'source_data', 'project', 'ext_data', 'processed', 'monarch')}",
                    "-o", f"{os.path.join(RETURNED_PATH, 'evidence_profiler', 'evidence_reports')}",
                    "-V", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'evidence_profiler', 'raw_data', 'variants')}"]
    evidence_profiler(list_of_args)
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'evidence_profiler', 'evidence_reports')}")

def test_cohort_analyzer():
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer')}", exist_ok=True)
    list_of_args = ["-i", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'paco_translator', 'cohort_toy_dataset.txt')}",
                    "-o", f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer')}", 
                    "-c", "chr", "-d", "patient_id", "-s", "start", "-e", "end", "-p", "phenotypes", "-S", "|",
                    "-t", "freq", "-D", "-m", "lin", "-a", "-n"]
    
    #Testing first time to check if calculations can be executed without errors
    cohort_analyzer(list_of_args)
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer.html')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer_lin_clusters.html')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'detailed_hpo_profile_evaluation.csv')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'rejected_records.txt')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'pat_hpo_matrix.npy')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'pat_hpo_matrix_x.lst')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'pat_hpo_matrix_y.lst')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'profiles_similarity_lin.txt')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'similarity_matrix_lin.npy')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'similarity_matrix_lin_x.lst')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'lin_raw_cls.npy')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'lin_linkage.npy')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'lin_clusters.txt')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'hpo_cohort_frequency.txt')}")
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'temp', 'cluster_asignation')}")            
    
    #Testing second time to check if files can be loaded (instead of doing the calculations again) without errors
    cohort_analyzer(list_of_args)

    shutil.rmtree(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer')}")

def test_profiles2phenopacket():
    os.environ["hpo_file"] = HPO_FILE
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'profiles2phenopacket')}", exist_ok=True)

    list_of_args = ["-i", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'paco_translator', 'cohort_toy_dataset.txt')}", 
                    "-c", "chr", "-d", "patient_id", "-s", "start", "-e", "end", "-p", "phenotypes", 
                    "-x", "sex", "-S", "|", "-n", "-o", f"{os.path.join(RETURNED_PATH, 'profiles2phenopacket')}"]
    profiles2phenopacket(list_of_args)

    for patient in ["132", "599", "647", "648"]:
        f1 = open(f"{os.path.join(EXPECTED_PATH, 'profiles2phenopacket', patient+'.json')}", "r")
        f2 = open(f"{os.path.join(RETURNED_PATH, 'profiles2phenopacket', patient+'.json')}", "r")
        expected = json.load(f1)
        returned = json.load(f2)
        expected["phenotypicFeatures"] = sorted(expected["phenotypicFeatures"], key=lambda phenotypes: phenotypes['type']['id'])
        returned["phenotypicFeatures"] = sorted(returned["phenotypicFeatures"], key=lambda phenotypes: phenotypes['type']['id'])
        assert expected == returned
        
        f1.close()
        f2.close()
        os.remove(f"{os.path.join(RETURNED_PATH, 'profiles2phenopacket', patient+'.json')}")

def test_paco_translator():
    os.environ["hpo_file"] = HPO_FILE
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'paco_translator')}", exist_ok=True)
    list_of_args = ["-P", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'paco_translator', 'cohort_toy_dataset.txt')}", 
                    "-c", "chr", "-d", "patient_id", "-s", "start", "-e", "end", "-p", "phenotypes", 
                    "-S", "|", "-n", "-C", "--n_phens", "2"]
    
    for file in ["untranslated_paco_file.txt", "translated_paco_file.txt", "untranslated_default_file.txt", "translated_default_file.txt"]:
        list_of_args_copy = list_of_args[:]
        if file == "untranslated_paco_file.txt": 
            list_of_args_copy.extend(["-o", f"{os.path.join(RETURNED_PATH, 'paco_translator', file)}", "-m", "paco"])
        elif file == "translated_paco_file.txt": 
            list_of_args_copy.extend(["-o", f"{os.path.join(RETURNED_PATH, 'paco_translator', file)}", "-m", "paco", "-t"])
        elif file == "untranslated_default_file.txt": 
            list_of_args_copy.extend(["-o", f"{os.path.join(RETURNED_PATH, 'paco_translator', file)}", "-m", "default"])
        elif file == "translated_default_file.txt": 
            list_of_args_copy.extend(["-o", f"{os.path.join(RETURNED_PATH, 'paco_translator', file)}", "-m", "default", "-t"])
        paco_translator(list_of_args_copy)
        f1 = open(f"{os.path.join(EXPECTED_PATH, 'paco_translator', file)}", "r")
        f2 = open(f"{os.path.join(RETURNED_PATH, 'paco_translator', file)}", "r")
        expected = strng2table(f1.read())
        returned = strng2table(f2.read())
        assert expected == returned
        
        f1.close()
        f2.close()
        os.remove(f"{os.path.join(RETURNED_PATH, 'paco_translator', file)}")

def test_get_sorted_profs():
    os.environ["hpo_file"] = HPO_FILE
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'get_sorted_profs')}", exist_ok=True)

    list_of_args = ["-P", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'get_sorted_profs', 'pmm2_paco_format.txt')}", 
                    "-r", f"{os.path.join(SCRIPT_DATA_TEST_PATH, 'get_sorted_profs', 'profile')}",
                    "-c", "chr", "-d", "patient_id", "-s", "start", "-e", "stop", "-p", "phenotypes", 
                    "-S", ",", "-o", f"{os.path.join(RETURNED_PATH, 'get_sorted_profs', 'report.html')}"]
    
    get_sorted_profs(list_of_args)
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'get_sorted_profs', 'report.html')}")

    os.remove(f"{os.path.join(RETURNED_PATH, 'get_sorted_profs', 'report.html')}")


def test_get_gen_features():
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'get_gen_features')}", exist_ok=True)    
    list_of_args = ["-i", f"{os.path.join(DATA_TEST_PATH, 'cohort_toy_dataset.txt')}", 
                "-c", "chr", "-d", "patient_id", "-s", "start", "-e", "end", "-t", "gene",
                "-r", f"{GENCODE}",
                "-o", f"{os.path.join(RETURNED_PATH, 'get_gen_features', 'patients_ensemble.txt')}"]

    #Using ensemble id
    get_gen_features(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'get_gen_features', 'patients_ensemble.txt')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'get_gen_features', 'patients_ensemble.txt')}", "r")
    expected = strng2table(f1.read())
    returned = strng2table(f2.read())
    assert expected == returned
    f1.close()
    f2.close()

    #Using gene name
    list_of_args[-1] = f"{os.path.join(RETURNED_PATH, 'get_gen_features', 'patients_geneName.txt')}"
    list_of_args.extend(["-n", "gene_name"])
    get_gen_features(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'get_gen_features', 'patients_geneName.txt')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'get_gen_features', 'patients_geneName.txt')}", "r")
    expected = strng2table(f1.read())
    returned = strng2table(f2.read())
    assert expected == returned
    f1.close()
    f2.close()

    os.remove(f"{os.path.join(RETURNED_PATH, 'get_gen_features', 'patients_ensemble.txt')}")
    os.remove(f"{os.path.join(RETURNED_PATH, 'get_gen_features', 'patients_geneName.txt')}")


