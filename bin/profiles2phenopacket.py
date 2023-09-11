#! /usr/bin/env python
import argparse, os, sys

ROOT_PATH=os.path.dirname(__file__)
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.common_optparse import Common_optparse
from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser

with open(CONSTANTS_PATH) as infile:
    exec(infile.read())

#############################
## METHODS
#############################
def load_index(path_index):
    vcf_index = {}
    with open(path_index) as f:
        for line in f:
            id, path = line.strip().split("\t")
            vcf_index[id] = path

    return vcf_index

############################################################################################
## OPTPARSE
############################################################################################
parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')
Common_optparse.add_options(parser)

parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                    help="Input file with patient data")

parser.add_argument("-I", "--vcf_index", dest="vcf_index", default= None,
                    help="VCF file with patient id pointing to vcf path")
                    
parser.add_argument("-o", "--output_file", dest="output_folder", default= None,
                    help="Output folder")

opts = parser.parse_args()
options = vars(opts)

#############################################################
## MAIN
#############################################################

hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
Cohort.load_ontology("hpo", hpo_file, options.get("excluded_hpo"))
Cohort.act_ont = "hpo"

patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
rejected_hpos_C, rejected_patients_C = patient_data.check(hard=True)
patient_data.link2ont(Cohort.act_ont)

vcf_index = None
if options.get("vcf_index"): vcf_index = load_index(options["vcf_index"])
patient_data.export_phenopackets(options["output_folder"], options["genome_assembly"], vcf_index= vcf_index)