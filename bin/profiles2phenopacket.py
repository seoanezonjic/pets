#! /usr/bin/env python
import argparse
import pets
from pets import Cohort_Parser, Cohort
import os, sys

ROOT_PATH=os.path.dirname(__file__)


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

#TODO: check how to add the common options (COMMON_OPTPARSE variable not found anywhere)
eval(open(COMMON_OPTPARSE).read())

parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                    help="Input file with patient data")

parser.add_argument("-I", "--vcf_index", dest="vcf_index", default= None,
                    help="VCF file with patient id pointing to vcf path")

parser.add_argument("-n", "--hpo_names", dest="names", default= False, action="store_true",
                    help="Define if the input HPO are human readable names. Default false")

parser.add_argument("-o", "--output_file", dest="output_folder", default= None,
                    help="Output folder")

parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                    help="Column name if header true or 0-based position of the column with the HPO terms")

parser.add_argument("-S", "--hpo_separator", dest="separator", default= "|",
                    help="Set which character must be used to split the HPO profile. Default '|'")

parser.add_argument("-s", "--start_col", dest="start_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate")

parser.add_argument("-h", "--help", help="Show this message", action="store_true")

opts = parser.parse_args()

if opts.help:
    parser.print_help()
    sys.exit()
    
options = vars(opts)

#############################################################
## MAIN
#############################################################

#TODO: ask Pedro about ENV variable
hpo_file = ENV['hpo_file'] if !ENV['hpo_file'].nil? else HPO_FILE
Cohort.load_ontology("hpo", hpo_file, options["excluded_hpo"])
Cohort.act_ont = "hpo"

patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
rejected_hpos_C, rejected_patients_C = patient_data.check(hard=True)
patient_data.link2ont(Cohort.act_ont)

if options.get("vcf_index"): vcf_index = load_index(options["vcf_index"])
patient_data.export_phenopackets(options["output_folder"], options["genome_assembly"], vcf_index= vcf_index)