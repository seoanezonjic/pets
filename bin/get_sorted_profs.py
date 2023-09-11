#! /usr/bin/env python

import os, sys
import argparse
from py_report_html import Py_report_html

ROOT_PATH = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort_analyser_methods import get_similarity_matrix
from pets.common_optparse import Common_optparse

CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
with open(CONSTANTS_PATH) as infile:
    exec(infile.read())

#############################################################################################
## OPTPARSE
############################################################################################

parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')
Common_optparse.add_options(parser)

parser.add_argument("-o", "--output_file", dest="output_file", default= 'report.html',
                    help="Output paco file with HPO names")

parser.add_argument("-P", "--input_file", dest="input_file", default= None,
                    help="Input file with PACO extension")

parser.add_argument("-f", "--general_prof_freq", dest="term_freq", default= 0, type= int,
                    help="When reference profile is not given, a general ine is computed with all profiles. If a freq is defined (0-1), all terms with freq minor than limit are removed")

parser.add_argument("-L", "--matrix_limits", dest="matrix_limits", default= [20, 40], type= lambda data: [int(i) for i in data.split(",")],
                    help="Number of rows and columns to show in heatmap defined as 'Nrows,Ncols'. Default 20,40")

parser.add_argument("-r", "--ref_profile", dest="ref_prof", default= None, 
                    type = lambda file: [line.strip() for line in open(file).readlines()],
                    help="Path to reference profile. One term code per line")


opts = parser.parse_args()    
options = vars(opts)


#############################################################################################
## MAIN
############################################################################################

hpo_file = os.environ.get('hpo_file') if os.environ.get('hpo_file') else HPO_FILE
Cohort.load_ontology("hpo", hpo_file)
Cohort.act_ont = "hpo"
hpo = Cohort.get_ontology(Cohort.act_ont)
patient_data, _, _ = Cohort_Parser.load(options)
patient_data.check(hard=True)

clean_profiles = patient_data.profiles

if options.get("ref_prof"):
  ref_profile = hpo.clean_profile_hard(options["ref_prof"])
else:
  ref_profile = patient_data.get_general_profile(options["term_freq"])

hpo.load_profiles({"ref": ref_profile}, reset_stored= True)

similarities = hpo.compare_profiles(external_profiles= clean_profiles, sim_type= "lin", bidirectional= False)

candidate_sim_matrix, candidates, candidates_ids = get_similarity_matrix(ref_profile, similarities["ref"], clean_profiles, hpo, options["matrix_limits"][0], options["matrix_limits"][-1])
candidate_sim_matrix.insert(0, ['HP'] + candidates_ids)

template = open(os.path.join(REPORT_FOLDER, 'similarity_matrix.txt')).read()
container = { "similarity_matrix": candidate_sim_matrix }
report = Py_report_html(container, 'Similarity matrix')
report.build(template)
report.write(options["output_file"])

with open(options["output_file"].replace('.html','') +'.txt', 'w') as f:
  for candidate, value in similarities["ref"].items():
    f.write("\t".join([str(candidate), str(value)])+"\n")
