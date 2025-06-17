import sys, os, json, unittest 
import subprocess
from importlib.resources import files

from subprocess import PIPE
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort
from pets import get_gen_features, paco_translator, profiles2phenopacket, cohort_analyzer, evidence_profiler, diseasome_generator, collapse_terms, filter_omim, MONDO_FILE, HPO_FILE, GENCODE
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
    assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'evidence_profiler', 'evidence_reports', 'pat_174.html')}")

    shutil.rmtree(f"{os.path.join(RETURNED_PATH, 'evidence_profiler')}")

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

    #Testing it can generate static, dynamic and canvas ontoplots
    for plot_type in ["static", "dynamic", "canvas"]:
        os.makedirs(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer')}", exist_ok=True)
        cohort_analyzer(list_of_args + ["--ontoplot_mode", plot_type])
        assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer.html')}")
        assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer_lin_clusters.html')}")
        shutil.rmtree(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer')}")

    #Testing it can generate detailed clusters heatmap with cluster, cohort and cohort_sort phenotypes on the y-axis
    for ySortFunc in ["cluster", "cohort", "cohort_sort"]:
        os.makedirs(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer')}", exist_ok=True)
        cohort_analyzer(list_of_args + ["--detailed_cluster_yaxis", ySortFunc])
        assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer.html')}")
        assert os.path.exists(f"{os.path.join(RETURNED_PATH, 'cohort_analyzer', 'cohort_analyzer_lin_clusters.html')}")
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

def test_diseasome_generator():
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'diseasome_generator')}", exist_ok=True)
    disorder_classes = files('pets.external_data').joinpath('disorder_classes')
    diseases_input = os.path.join(SCRIPT_DATA_TEST_PATH, 'diseasome_generator', 'diseases')
    diseasome_input = os.path.join(SCRIPT_DATA_TEST_PATH, 'diseasome_generator', 'diseasome')

    # Generate diseasome
    list_of_args = ["-i", str(diseases_input), "-C", str(disorder_classes), "-O", str(MONDO_FILE), "-o", f"{os.path.join(RETURNED_PATH, 'diseasome_generator','diseasome')}", "-g"]
    diseasome_generator(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'diseasome_generator', 'diseasome')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'diseasome_generator', 'diseasome')}", "r")
    expected = sorted(strng2table(f1.read()),key = lambda x: x[1])
    returned = sorted(strng2table(f2.read()),key = lambda x: x[1])
    assert expected == returned
    f1.close()
    f2.close()

    # Analyze diseasome
    list_of_args = ["-D", str(diseasome_input),"-O", str(MONDO_FILE), "-o", f"{os.path.join(RETURNED_PATH, 'diseasome_generator','diseasome')}", "-A"]
    diseasome_generator(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'diseasome_generator', 'diseasome_analysis')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'diseasome_generator', 'diseasome_analysis')}", "r")
    expected = sorted(sorted(strng2table(f1.read()), key = lambda x : x[0]), key = lambda x : x[1])
    returned = sorted(sorted(strng2table(f2.read()), key = lambda x : x[0]), key = lambda x : x[1])
    assert expected == returned
    f1.close()
    f2.close()


def test_collapse_terms():
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'collapse_terms')}", exist_ok=True)
    terms2collapse = os.path.join(SCRIPT_DATA_TEST_PATH, 'collapse_terms', 'terms2collapse')
    terms2txt = os.path.join(SCRIPT_DATA_TEST_PATH, 'collapse_terms', 'terms2txt')

    # collapse terms, no uniq parent no just leaves.
    list_of_args = ["-i", str(terms2collapse), "-n", str(terms2txt),"-o", f"{os.path.join(RETURNED_PATH, 'collapse_terms','collapsed_terms')}", "-t", "0.8"]
    collapse_terms(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'collapse_terms', 'collapsed_terms')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'collapse_terms', 'collapsed_terms')}", "r")
    expected = sorted(strng2table(f1.read()),key = lambda x: x[0])
    returned = sorted(strng2table(f2.read()),key = lambda x: x[0])
    assert expected == returned
    f1.close()
    f2.close()

    # TODO: A good example is needed 
    # # collapse terms, uniq parent, no just leaves.
    # list_of_args = ["-i", str(terms2collapse), "-n", str(terms2txt),"-o", f"{os.path.join(RETURNED_PATH, 'collapse_terms','collapsed_terms_uniq_parents')}", "-t", "0.8", "-u"]
    # collapse_terms(list_of_args)
    # f1 = open(f"{os.path.join(EXPECTED_PATH, 'collapse_terms', 'collapsed_terms_uniq_parents')}", "r")
    # f2 = open(f"{os.path.join(RETURNED_PATH, 'collapse_terms', 'collapsed_terms_uniq_parents')}", "r")
    # expected = sorted(strng2table(f1.read()),key = lambda x: x[0])
    # returned = sorted(strng2table(f2.read()),key = lambda x: x[0])
    # assert expected == returned
    # f1.close()
    # f2.close()

    # collapse terms, uniq parent, no just leaves.
    list_of_args = ["-i", str(terms2collapse), "-n", str(terms2txt),"-o", f"{os.path.join(RETURNED_PATH, 'collapse_terms','collapsed_terms_just_leaves')}", "-t", "0.8", "-l"]
    collapse_terms(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'collapse_terms', 'collapsed_terms_just_leaves')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'collapse_terms', 'collapsed_terms_just_leaves')}", "r")
    expected = sorted(strng2table(f1.read()),key = lambda x: x[0])
    returned = sorted(strng2table(f2.read()),key = lambda x: x[0])
    assert expected == returned
    f1.close()
    f2.close()

def test_filter_omim():
    os.makedirs(f"{os.path.join(RETURNED_PATH, 'filter_omim')}", exist_ok=True)
    morbid_file = os.path.join(SCRIPT_DATA_TEST_PATH, 'filter_omim', 'morbidmap.txt')
    list_of_args = ["-i", str(morbid_file),"-o", f"{os.path.join(RETURNED_PATH, 'filter_omim','filtered_morbid_file')}"]
    filter_omim(list_of_args)
    f1 = open(f"{os.path.join(EXPECTED_PATH, 'filter_omim', 'filtered_morbid_file')}", "r")
    f2 = open(f"{os.path.join(RETURNED_PATH, 'filter_omim', 'filtered_morbid_file')}", "r")
    expected = strng2table(f1.read())
    returned = strng2table(f2.read())
    assert expected == returned
    f1.close()
    f2.close()