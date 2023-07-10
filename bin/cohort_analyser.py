import argparse, os, sys

ROOT_PATH=os.path.dirname(__file__)
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.common_optparse import Common_optparse
from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser

with open(CONSTANTS_PATH) as infile:
    exec(infile.read())


############################################################################################
## OPTPARSE
############################################################################################
def tolist(string): return string.split(',')

parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')
Common_optparse.add_options(parser)

parser.add_argument("-a","--coverage_analysis", dest="coverage_analysis", default= True, action="store_false",
  help="Deactivate genome coverage analysis. Default true")
parser.add_argument("-b", "--bin_size", dest="bin_size", default= 50000, type=int,
  help="Maximum number of bins to plot the coverage")
parser.add_argument("-C", "--clusters2show", dest="clusters2show_detailed_phen_data", default= 3, type=int,
  help="How many patient clusters are show in detailed phenotype cluster data section. Default 3")
parser.add_argument("-D","--detailed_clusters", dest="detailed_clusters", default= False, action="store_true",
  help="Show detiled cluster comparation using heatmaps. Default false")
parser.add_argument("-E", "--excluded_hpo", dest="--excluded_hpo", default= None,
  help="Path to file with a list of HPO phenotypes to exclude (low informative)")
parser.add_argument("-M", "--minClusterProportion", dest="minClusterProportion", default= 0.01, type=float,
  help="Minimum percentage of patients per cluster")
parser.add_argument("-f", "--patients_filter", dest="patients_filter", default= 2, type=int,
  help="Minimum number of patients sharing SORs. Default 2")
parser.add_argument("-g", "--clusters2graph", dest="clusters2graph", default= 30, type=int,
  help="How may patient clusters are plotted in cluster plots. Default 30")
parser.add_argument("-i", "--input_file", dest="input_file", default= None,
  help="Input file with patient data")
parser.add_argument("-m", "--clustering_methods", dest="clustering_methods", default=['resnik', 'jiang_conrath', 'lin'], type=tolist,
  help="Clustering methods")
parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
  help="Define if the input HPO are human readable names. Default false")
parser.add_argument("-o", "--output_file", dest="output_file", default= None,
  help="Output file with patient data")
parser.add_argument("-P", "--hpo_file", dest="hpo_file", default= None,
  help="Input HPO file for extracting HPO codes")
parser.add_argument("-p", "--hpo_term_col", dest="hpo_term_col", default= None,
  help="Column name if header true or 0-based position of the column with the HPO terms")
parser.add_argument("-S", "--hpo_separator", dest="separator", default='|',
  help="Set which character must be used to split the HPO profile. Default '|'")
parser.add_argument("-s", "--start_col", dest="start_col", default= None,
  help="Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate")
parser.add_argument("-r", "--root_node", dest="root_node", default= "HP:0000118",
  help="Root node from which the ontology will be represented")
parser.add_argument("-t", "--ic_stats", dest="ic_stats", default= "freq",
  help="'freq' to compute IC based en hpo frequency in the input cohort. 'freq_internal' to use precomputed internal IC stats. 'onto' to compute ic based on hpo ontology structure.. Default freq")
parser.add_argument("-T", "--threads", dest="threads", default= 1, type=int,
  help="Number of threads to be used in calculations. Default 1")
parser.add_argument("--reference_profiles", dest="reference_profiles", default= None,
  help="Path to file tabulated file with first column as id profile and second column with ontology terms separated by separator.")
parser.add_argument("--sim_thr", dest="sim_thr", default=None, type=float,
  help="Keep pairs with similarity value >= FLOAT.")
