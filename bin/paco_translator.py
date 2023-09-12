#! /usr/bin/env python

import os, sys, argparse
ROOT_PATH=os.path.dirname(__file__)
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser
from pets.common_optparse import Common_optparse

with open(CONSTANTS_PATH) as infile:
    exec(infile.read())

############################################################################################
## OPTPARSE
############################################################################################
parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')
Common_optparse.add_options(parser)

parser.add_argument("-o", "--output_file", dest="output_file", default= None,
                    help="Output paco file with HPO names")

parser.add_argument("-P", "--input_file", dest="input_file", default= None,
                    help="Input file with PACO extension")

parser.add_argument("--n_phens", dest="n_phens", default= None, type=int,
                    help="Remove records with N or less phenotypes")

parser.add_argument("-m","--save_mode", dest="save_mode", default= "default",
                    help="Set output data mode")

parser.add_argument("-t","--translate", dest="translate", default= False, action="store_true",
                    help="Set to translate from hpo codes to names. By default, ther is not translation")

parser.add_argument("-C", "--clean_PACO", dest="clean_PACO", default= False, action="store_true",
                    help="Clean PACO files")

parser.add_argument("-r", "--removed_path", dest="removed_path", default= None,
                    help="Desired path to write removed profiles file")

opts = parser.parse_args() 
options = vars(opts)

###############
#MAIN
###############
hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
Cohort.load_ontology("hpo", hpo_file, options.get("excluded_hpo"))
Cohort.act_ont = "hpo"

patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(options)

if options.get("clean_PACO"):
    removed_terms, removed_profiles = patient_data.check(hard=True)
    if options.get("removed_path") and removed_profiles != None and len(removed_profiles) > 0:
        rejected_file = os.path.basename(options["input_file"]).split(".")[0] +'_excluded_patients'
        file = os.path.join(options["removed_path"], rejected_file)
        with open(file, 'w') as f:
            for profile in removed_profiles:
                f.write(profile+'\n')

if options.get("n_phens"): rejected_patients_by_phen = patient_data.filter_by_term_number(options["n_phens"])
patient_data.save(options["output_file"], options["save_mode"], options["translate"])
