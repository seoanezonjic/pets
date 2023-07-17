#! /usr/bin/env python
import argparse, os, sys
import re
import numpy as np
import time
from collections import defaultdict
import csv

ROOT_PATH=os.path.dirname(__file__)
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.common_optparse import Common_optparse
from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser

with open(CONSTANTS_PATH) as infile:
    exec(infile.read())

############################################################################################
## METHODS
############################################################################################
def load_hpo_ci_values(information_coefficient_file):
  hpos_ci_values = {}
  with open(information_coefficient_file) as f:
    for line in f:
      hpo_code, ci = line.rstrip().split("\t")
      hpos_ci_values[hpo_code] = float(ci)
  return hpos_ci_values

def get_profile_ic(hpo_names, phenotype_ic):
  ic = 0
  profile_length = 0
  for hpo_id in hpo_names:
    hpo_ic = phenotype_ic.get(hpo_id)
    if hpo_ic == None: raise Exception(f"The term {hpo_id} not exists in the given ic table")
    ic += hpo_ic 
    profile_length += 1
  if profile_length == 0: profile_length = 1 
  return ic/profile_length

def system_call(code_folder, script, args_string):
  cmd = f"{os.path.join(code_folder, script)} {args_string}"
  print(f"==> {cmd}")
  start = time.time()
  os.system(cmd)
  print(f"Execution time: {time.time() - start}")

def dummy_cluster_patients(patient_data, matrix_file, clust_pat_file):
  if not os.path.exists(matrix_file):
    pat_hpo_matrix, pat_id, hp_id  = to_bmatrix(patient_data)
    x_axis_file = re.sub('.npy','_x.lst', matrix_file)
    y_axis_file = re.sub('.npy','_y.lst', matrix_file)
    save(pat_hpo_matrix, matrix_file, hp_id, x_axis_file, pat_id, y_axis_file)
  if not os.path.exists(clust_pat_file):
    system_call(EXTERNAL_CODE, 'get_clusters.R', f"-d {matrix_file} -o {clust_pat_file} -y {re.sub('.npy','', matrix_file)}")
  clustered_patients = load_clustered_patients(clust_pat_file)
  return(clustered_patients)

def load_clustered_patients(file):
  clusters = {}
  with open(file) as f:
    for line in f:
      pat_id, cluster_id = line.rstrip().split("\t")
      query = clusters.get(cluster_id)
      if query == None:
        clusters[cluster_id] = [pat_id]
      else:
        query.append(pat_id)
  return clusters

def process_dummy_clustered_patients(options, clustered_patients, patient_data, phenotype_ic): # get ic and chromosomes
  ont = Cohort.get_ontology(Cohort.act_ont)
  all_ics = []
  all_lengths = []
  top_cluster_phenotypes = []
  cluster_data_by_chromosomes = []
  multi_chromosome_patients = 0
  processed_clusters = 0
  for cluster_id, patient_ids in sorted(list(clustered_patients.items()), key=lambda x: len(x[1]), reverse=True):
    num_of_patients = len(patient_ids)
    if num_of_patients == 1: continue 
    chrs, all_phens, profile_ics, profile_lengths = process_cluster(patient_ids, patient_data, phenotype_ic, options, ont, processed_clusters)
    if processed_clusters < options['clusters2show_detailed_phen_data']: top_cluster_phenotypes.append(all_phens)
    all_ics.append(profile_ics)
    all_lengths.append(profile_lengths)
    if not options.get('chromosome_col') == None:
      if len(chrs) > 1: multi_chromosome_patients += num_of_patients
      for chrm, count in chrs.items(): cluster_data_by_chromosomes.append([cluster_id, num_of_patients, chrm, count])
    processed_clusters += 1
  return all_ics, all_lengths, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients

def process_cluster(patient_ids, patient_data, phenotype_ic, options, ont, processed_clusters):
  chrs = defaultdict(lambda: 0)
  all_phens = []
  profile_ics = []
  profile_lengths = []
  for pat_id in patient_ids:
    phenotypes = patient_data.get_profile(pat_id) 
    profile_ics.append(get_profile_ic(phenotypes, phenotype_ic))
    profile_lengths.append(len(phenotypes))
    if processed_clusters < options['clusters2show_detailed_phen_data']:
      phen_names, rejected_codes = ont.translate_ids(phenotypes) #optional
      all_phens.append(phen_names)
    if not options.get('chromosome_col') == None: 
      for chrm in patient_data.get_vars(pat_id).get_chr(): chrs[chrm] += 1
  return chrs, all_phens, profile_ics, profile_lengths 

def get_summary_stats(patient_data, rejected_patients, hpo_stats, fraction_terms_specific_childs, rejected_hpos):
  stats = [
    ['Unique HPO terms', len(hpo_stats)],
    ['Cohort size', len(patient_data.profiles)],
    ['Rejected patients by empty profile', len(rejected_patients)],
    ['HPOs per patient (average)', patient_data.get_profiles_mean_size()],
    ['HPO terms per patient: percentile 90', patient_data.get_profile_length_at_percentile(perc=90)],
    ['Percentage of HPO with more specific children', round((fraction_terms_specific_childs * 100), 4)],
    ['DsI for uniq HP terms', patient_data.get_dataset_specifity_index('uniq')],
    ['DsI for frequency weigthed HP terms', patient_data.get_dataset_specifity_index('weigthed')],
    ['Number of unknown phenotypes', len(rejected_hpos)]
  ]
  return stats

def get_mean_size(all_sizes):
  accumulated_size = 0
  number = 0
  for size, occurrences in all_sizes: 
    accumulated_size += size *occurrences
    number += occurrences
  return accumulated_size /number

def calculate_coverage(regions_data, delete_thresold = 0):
  raw_coverage = {}
  n_regions = 0
  patients = 0
  nt = 0
  for start, stop, chrm, reg_id in regions_data:
    number_of_patients = int(reg_id.split('.')[-1])
    if number_of_patients <= delete_thresold:
      number_of_patients = 0
    else:
      n_regions += 1
      nt += stop - start      
    add_record(raw_coverage, chrm, [start, stop, number_of_patients])
    patients += number_of_patients
  return raw_coverage, n_regions, nt, patients /n_regions 

def get_final_coverage(raw_coverage, bin_size):
  coverage_to_plot = []
  for chrm, coverages in raw_coverage.items():
    for start, stop, coverage in coverages:
      bin_start = start - start % bin_size
      bin_stop = stop - stop%bin_size
      while bin_start < bin_stop:
        coverage_to_plot.append([chrm, bin_start, coverage])
        bin_start += bin_size
  return coverage_to_plot

def get_sor_length_distribution(raw_coverage):
  all_cnvs_length = []
  cnvs_count = []
  for chrm, coords_info in raw_coverage.items():
    for start, stop, pat_records in coords_info:
      region_length = stop - start + 1
      all_cnvs_length.append([region_length, pat_records])
  all_cnvs_length.sort(key=lambda x: x[1])
  return all_cnvs_length

### AUX

def add_record(dictio, key, record, uniq=False):
  query = dictio.get(key)
  if query == None:
    dictio[key] = [record]
  elif not uniq: # We not take care by repeated entries
    query.append(record)
  elif not record in query: # We want uniq entries
    query.append(record)

def get_hash_values_idx(dictio):
  x_names_indx = {}
  i = 0
  for k, values in dictio.items():
    for val_id in values:
      if type(val_id) is list: val_id = val_id[0]
      query = x_names_indx.get(val_id)
      if query == None:
        x_names_indx[val_id] = i
        i += 1
  return x_names_indx

def to_bmatrix(dictio):
  x_names_indx = get_hash_values_idx(dictio)
  y_names = list(dictio.keys())
  x_names = list(x_names_indx.keys())
   # row (y), cols (x)
  matrix = np.zeros((len(dictio), len(x_names)))
  i = 0
  for id, items in dictio.items():
    for item_id in items: matrix[i, x_names_indx[item_id]] = 1
    i += 1
  return matrix, y_names, x_names

def save(matrix, matrix_filename, x_axis_names=None, x_axis_file=None, y_axis_names=None, y_axis_file=None):
  if not x_axis_names == None:
    with open(x_axis_file, 'w') as f:
      f.write("\n".join(x_axis_names))
  if not y_axis_names == None:
    with open(y_axis_file, 'w') as f:
      f.write("\n".join(y_axis_names))
  np.save(matrix_filename, matrix)

def write_detailed_hpo_profile_evaluation(suggested_childs, detailed_profile_evaluation_file, summary_stats):
  with open(detailed_profile_evaluation_file, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for pat_id, suggestions in suggested_childs.items():
      warning = None
      if len(suggestions) < 4: warning = 'WARNING: Very few phenotypes' 
      csvwriter.writerow([f"PATIENT {pat_id}", f"{warning}"])
      csvwriter.writerow(["CURRENT PHENOTYPES", "PUTATIVE MORE SPECIFIC PHENOTYPES"])
      for parent, childs in suggestions:
        parent_code, parent_name = parent
        if len(childs) == 0:
          csvwriter.writerow([f"{parent_name} ({parent_code})", '-'])
        else:
          parent_writed = False
          for child_code, child_name in childs:
            if not parent_writed:
              parent_field = f"{parent_name} ({parent_code})"
              parent_writed = True
            else:
              parent_field = ""
            csvwriter.writerow([parent_field, f"{child_name} ({child_code})"])
      csvwriter.writerow(["", ""])

def write_arrays4scatterplot(x_axis_value, y_axis_value, filename, x_axis_name, y_axis_name):
  with open(filename, 'w') as f:
    f.write(f"{x_axis_name}\t{y_axis_name}\n")
    for i, value in enumerate(x_axis_value):
      y_value = y_axis_value[i]
      if y_value == None: raise Exception(f"The {i} position is not presented in y_axis_value") 
      f.write(f"{value}\t{y_value}\n")

def write_cluster_ic_data(all_ics, profile_lengths, cluster_ic_data_file, limit):
  with open(cluster_ic_data_file, 'w') as f:
    f.write("\t".join(['cluster_id', 'ic', 'Plen']) + "\n")
    for i, cluster_ics in enumerate(all_ics):
      if i == limit: break
      cluster_length = len(cluster_ics)
      for j, clust_ic in enumerate(cluster_ics):
        f.write(f"{cluster_length}_{i}\t{clust_ic}\t{profile_lengths[i][j]}\n")

def write_cluster_chromosome_data(cluster_data, cluster_chromosome_data_file, limit):
  with open(cluster_chromosome_data_file, 'w') as f:
    f.write("\t".join(['cluster_id', 'chr', 'count']) + "\n")
    index = 0
    if len(cluster_data) > 0: last_id = cluster_data[0][0] 
    for cluster_id, patient_number, chrm, count in cluster_data:
      if cluster_id != last_id: index += 1 
      if index == limit: break 
      f.write("\t".join(["{patient_number}_{index}", chrm, str(count)]) + "\n")
      last_id = cluster_id

def write_coverage_data(coverage_to_plot, coverage_to_plot_file):
  with open(coverage_to_plot_file, 'w') as f:
    for chrm, position, freq in coverage_to_plot:
      f.write(f"{chrm}\t{position}\t{freq}")

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
