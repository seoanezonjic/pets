#! /usr/bin/env python
import argparse, os, sys

ROOT_PATH=os.path.dirname(__file__)
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.common_optparse import Common_optparse
from pets.cohort_analyser_methods import *

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
  help="Show detailed cluster comparation using heatmaps. Default false")
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
parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
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
options = parser.parse_args()

##########################
#MAIN
##########################
opts = vars(options)
if opts['genome_assembly'] == 'hg19' or opts['genome_assembly'] == 'hg37':
  CHR_SIZE = os.path.join(EXTERNAL_DATA, 'chromosome_sizes_hg19.txt')
elif opts['genome_assembly'] == 'hg38':
  CHR_SIZE = os.path.join(EXTERNAL_DATA, 'chromosome_sizes_hg38.txt')
elif opts['genome_assembly'] == 'hg18':
  CHR_SIZE = os.path.join(EXTERNAL_DATA, 'chromosome_sizes_hg18.txt')
else:
  raise Exception('Wrong human genome assembly. Please choose between hg19, hg18 or hg38.')

output_folder = os.path.dirname(opts['output_file'])
detailed_profile_evaluation_file = os.path.join(output_folder, 'detailed_hpo_profile_evaluation.csv')
rejected_file = os.path.join(output_folder, 'rejected_records.txt')
temp_folder = os.path.join(output_folder, 'temp')
matrix_file = os.path.join(temp_folder, 'pat_hpo_matrix.npy')
hpo_ic_file = os.path.join(temp_folder, 'hpo_ic.txt')
hpo_profile_ic_file = os.path.join(temp_folder, 'hpo_ic_profile.txt')
hpo_frequency_file = os.path.join(temp_folder, 'hpo_cohort_frequency.txt')
parents_per_term_file = os.path.join(temp_folder, 'parents_per_term.txt')
clustered_patients_file = os.path.join(temp_folder, 'cluster_asignation')
cluster_ic_data_file = os.path.join(temp_folder, 'cluster_ic_data.txt')
cluster_chromosome_data_file = os.path.join(temp_folder, 'cluster_chromosome_data.txt')
coverage_to_plot_file = os.path.join(temp_folder, 'coverage_data.txt')
sor_coverage_to_plot_file = os.path.join(temp_folder, 'sor_coverage_data.txt')
ronto_file = os.path.join(temp_folder, 'hpo_freq_colour')

if not os.path.exists(temp_folder): os.mkdir(temp_folder)

hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
Cohort.load_ontology("hpo", hpo_file, opts.get("excluded_hpo"))
Cohort.act_ont = "hpo"

patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(opts)
rejected_hpos_C, rejected_patients_C = patient_data.check()
rejected_hpos = list(set(rejected_hpos_L).union(rejected_hpos_C))
rejected_patients = list(set(rejected_patients_L).union(rejected_patients_C))
with open(rejected_file, 'w') as f: f.write("\n".join(rejected_patients))

patient_data.link2ont(Cohort.act_ont) # TODO: check if method load should call to this and use the semtools checking methods (take care to only remove invalid terms)

profile_sizes, parental_hpos_per_profile = patient_data.get_profile_redundancy()
patient_data.check(hard=True)
hpo_stats = patient_data.get_profiles_terms_frequency() # hpo NAME, freq
for stat in hpo_stats: stat[1] = stat[1]*100
with open(hpo_frequency_file, 'w') as f:
  for hpo_code, freq in patient_data.get_profiles_terms_frequency(translate= False): # hpo CODE, freq
    f.write(f"{hpo_code}\t{freq}\n")

suggested_childs, fraction_terms_specific_childs = patient_data.compute_term_list_and_childs()
ontology_levels, distribution_percentage = patient_data.get_profile_ontology_distribution_tables()
onto_ic, freq_ic, onto_ic_profile, freq_ic_profile = patient_data.get_ic_analysis()

if opts['ic_stats'] == 'freq_internal':
  ic_file = os.environ['ic_file'] if os.environ.get('ic_file') else IC_FILE
  freq_ic = load_hpo_ci_values(ic_file)
  phenotype_ic = freq_ic
  freq_ic_profile = {}
  for pat_id, phenotypes in patient_data.each_profile():
    freq_ic_profile[pat_id] = get_profile_ic(phenotypes, phenotype_ic)
elif opts['ic_stats'] == 'freq':
  phenotype_ic = freq_ic
elif opts['ic_stats'] == 'onto':
  phenotype_ic = onto_ic

clustered_patients = dummy_cluster_patients(patient_data.profiles, matrix_file, clustered_patients_file)
all_ics, prof_lengths, clust_by_chr, top_clust_phen, multi_chr_pats = process_dummy_clustered_patients(opts, clustered_patients, patient_data, phenotype_ic)

summary_stats = get_summary_stats(patient_data, rejected_patients, hpo_stats, fraction_terms_specific_childs, rejected_hpos)

all_cnvs_length = []
if not opts.get('chromosome_col') == None:
  summary_stats.append(['Number of clusters with mutations accross > 1 chromosomes', multi_chr_pats])
  
  #----------------------------------
  # Prepare data to plot coverage
  #----------------------------------
  if opts.get('coverage_analysis'):
    patient_data.index_vars()
    all_cnvs_length = patient_data.get_vars_sizes(True)
    cnv_size_average = get_mean_size(all_cnvs_length)
    patients_by_cluster, sors = patient_data.generate_cluster_regions('reg_overlap', 'A', 0)

    ###1. Process CNVs
    raw_coverage, n_cnv, nt, pats_per_region = calculate_coverage(sors)
    summary_stats.extend(
      [['Average variant size', round(cnv_size_average, 4)],
      ['Nucleotides affected by mutations', nt],
      ['Number of genome windows', n_cnv],
      ['Mean patients per genome window', round(pats_per_region, 4)]])
    coverage_to_plot = get_final_coverage(raw_coverage, opts['bin_size'])

    ###2. Process SORs
    raw_sor_coverage, n_sor, nt, pats_per_region = calculate_coverage(sors, opts['patients_filter'] - 1)
    summary_stats.extend(
      [[f"Number of genome window shared by >= {opts['patients_filter']} patients", n_sor],
      ["Number of patients with at least 1 SOR", len(patients_by_cluster)],
      ['Nucleotides affected by mutations', nt]])
    # ['Patient average per region', pats_per_region]
    sor_coverage_to_plot = get_final_coverage(raw_sor_coverage, opts['bin_size'])

    all_sor_length = get_sor_length_distribution(raw_sor_coverage)  

#----------------------------------
# Write files for report
#----------------------------------
write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats)
write_arrays4scatterplot(list(onto_ic.values()), list(freq_ic.values()), hpo_ic_file, 'OntoIC', 'FreqIC') # hP terms
write_arrays4scatterplot(list(onto_ic_profile.values()), list(freq_ic_profile.values()), hpo_profile_ic_file, 'OntoIC', 'FreqIC') #HP profiles
write_arrays4scatterplot(profile_sizes, parental_hpos_per_profile, parents_per_term_file, 'ProfileSize', 'ParentTerms')
write_cluster_ic_data(all_ics, prof_lengths, cluster_ic_data_file, opts['clusters2graph'])

if not os.path.exists(os.path.join(temp_folder, 'hpo_ics.pdf')): system_call(EXTERNAL_CODE, 'plot_scatterplot_simple.R', f"-i {hpo_ic_file} -o {os.path.join(temp_folder, 'hpo_ics.pdf')} -x 'OntoIC' -y 'FreqIC' --x_tag 'HP Ontology IC' --y_tag 'HP Frequency based IC' --x_lim '0,4.5' --y_lim '0,4.5'") 
if not os.path.exists(os.path.join(temp_folder, 'hpo_profile_ics.pdf')): system_call(EXTERNAL_CODE, 'plot_scatterplot_simple.R', f"-i {hpo_profile_ic_file} -o {os.path.join(temp_folder, 'hpo_profile_ics.pdf')} -x 'OntoIC' -y 'FreqIC' --x_tag 'HP Ontology Profile IC' --y_tag 'HP Frequency based Profile IC' --x_lim '0,4.5' --y_lim '0,4.5'")
system_call(EXTERNAL_CODE, 'plot_scatterplot_simple.R', f"-i {parents_per_term_file} -o {os.path.join(temp_folder, 'parents_per_term.pdf')} -x 'ProfileSize' -y 'ParentTerms' --x_tag 'Patient HPO profile size' --y_tag 'Parent HPO terms within the profile'")
if not os.path.exists(ronto_file + '.png'): system_call(EXTERNAL_CODE, 'ronto_plotter.R', f"-i {hpo_frequency_file} -o {ronto_file} --root_node {opts['root_node']} -O {re.sub('.json','.obo', hpo_file)}") ###Cohort frequency calculation
system_call(EXTERNAL_CODE, 'plot_boxplot.R', f"{cluster_ic_data_file} {temp_folder} cluster_id ic 'Cluster size/id' 'Information coefficient' 'Plen' 'Profile size'")


if not opts.get('chromosome_col') == None:
  write_cluster_chromosome_data(clust_by_chr, cluster_chromosome_data_file, opts['clusters2graph'])
  system_call(EXTERNAL_CODE, 'plot_scatterplot.R', f"{cluster_chromosome_data_file} {temp_folder} cluster_id chr count 'Cluster size/id' 'Chromosome' 'Patients'")
  if opts['coverage_analysis']:
    ###1. Process CNVs
    write_coverage_data(coverage_to_plot, coverage_to_plot_file)
    system_call(EXTERNAL_CODE, 'plot_area.R', f"-d {coverage_to_plot_file} -o {temp_folder}/coverage_plot -x V2 -y V3 -f V1 -H -m {CHR_SIZE} -t CNV")
    ###2. Process SORs
    write_coverage_data(sor_coverage_to_plot, sor_coverage_to_plot_file)
    system_call(EXTERNAL_CODE, 'plot_area.R', f"-d {sor_coverage_to_plot_file} -o {temp_folder}/sor_coverage_plot -x V2 -y V3 -f V1 -H -m {CHR_SIZE} -t SOR")

#----------------------------------
# CLUSTER COHORT ANALYZER REPORT
#----------------------------------
get_semantic_similarity_clustering(opts, patient_data, temp_folder, os.path.join(REPORT_FOLDER, 'cluster_report.txt'), EXTERNAL_CODE)

#----------------------------------
# GENERAL COHORT ANALYZER REPORT
#----------------------------------
new_cluster_phenotypes = get_top_dummy_clusters_stats(top_clust_phen)

container = {
  'temp_folder' : temp_folder,
  # 'top_clust_phen' : len(top_clust_phen),
  'summary_stats' : summary_stats,
  'clustering_methods' : opts['clustering_methods'],
  'hpo_stats' : hpo_stats,
  'all_cnvs_length' : all_cnvs_length,
  'all_sor_length' : all_sor_length,
  'new_cluster_phenotypes' : len(new_cluster_phenotypes),
  'ontology_levels' : ontology_levels,
  'distribution_percentage' : distribution_percentage
}

clust_info = []
for clusterID, info in new_cluster_phenotypes.items():
    phens = ', '.join(info[1])
    freqs = ', '.join([ str(round(a,4)) for a in info[2]])
    container[f"clust_{clusterID}"] = [[info[0], phens, freqs]]

report = Py_report_html(container)
report.build(open(os.path.join(REPORT_FOLDER, 'cohort_report.txt')).read())
report.write(opts["output_file"] + '.html')