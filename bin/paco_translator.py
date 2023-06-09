#! /usr/bin/env python

import pets
from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser
from pets.constants import HPO_FILE
import os, sys, argparse

ROOT_PATH = os.path.dirname(__file__)

############################################################################################
## OPTPARSE
############################################################################################
parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')

parser.add_argument("-c", "--chromosome_col", dest="chromosome_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the chromosome")

parser.add_argument("-d", "--pat_id_col", dest="id_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the patient id")

parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")

parser.add_argument("-H", "--header", dest="header", default= None,
                    help="File has a line header. Default true")

parser.add_argument("-o", "--output_file", dest="output_file", default= None,
                    help="Output paco file with HPO names")

parser.add_argument("-P", "--input_file", dest="input_file", default= None,
                    help="Input file with PACO extension")

parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                    help="Column name if header true or 0-based position of the column with the HPO terms")

parser.add_argument("-s", "--start_col", dest="start_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate")

parser.add_argument("-S", "--hpo_separator", dest="separator", default= '|',
                    help="Set which character must be used to split the HPO profile. Default '|'")

parser.add_argument("--n_phens", dest="n_phens", default= None, type=int,
                    help="Remove records with N or less phenotypes")

parser.add_argument("-m","--save_mode", dest="save_mode", default= "default",
                    help="Set output data mode")

parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
                    help="Define if the input HPO are human readable names. Default false")

parser.add_argument("-t","--translate", dest="translate", default= False, action="store_true",
                    help="Set to translate from hpo codes to names. By default, ther is not translation")

opts = parser.parse_args() 
options = vars(opts)

###############
#MAIN
###############
hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
Cohort.load_ontology("hpo", hpo_file)
Cohort.act_ont = "hpo"

patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(options)
if options.get("n_phens"): rejected_patients_by_phen = patient_data.filter_by_term_number(options["n_phens"])
patient_data.save(options["output_file"], options["save_mode"], options["translate"])
