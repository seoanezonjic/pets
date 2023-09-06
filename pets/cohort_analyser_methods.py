import os, sys
import re
import numpy as np
import time
from collections import defaultdict
import csv

from py_report_html import Py_report_html
from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser

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

def get_top_dummy_clusters_stats(top_clust_phen):
  new_cluster_phenotypes = {}
  for clusterID, cluster in enumerate(top_clust_phen):
    phenotypes_frequency = defaultdict(lambda: 0)
    total_patients = len(cluster)
    for phenotypes in cluster:
      for p in phenotypes: phenotypes_frequency[p] += 1
    values = [ v / total_patients * 100 for v in phenotypes_frequency.values() ]
    new_cluster_phenotypes[clusterID] = [total_patients, list(phenotypes_frequency.keys()), values ]
  return new_cluster_phenotypes

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

def to_wmatrix(dictio, squared = True, symm = True):
  if squared:
    matrix, element_names = to_wmatrix_squared(dictio, symm=symm)
    return matrix, element_names
  else:
    matrix, y_names, x_names = to_wmatrix_rectangular(dictio, symm=symm)
    return matrix, y_names, x_names

def to_wmatrix_squared(dictio, symm = True):
  element_names = list(dictio.keys())
  matrix = np.zeros((len(dictio), len(dictio)))
  i = 0
  for elementA, relations in dictio.items():
    for j, elementB in enumerate(element_names):
      if elementA != elementB:
        query = relations.get(elementB)
        if query != None:
          matrix[i, j] = query
        elif symm: # TODO: PSZ, lo q se hace aqui no me cuadra
          matrix[i, j] = dictio[elementB][elementA]
    i += 1
  return matrix, element_names

def to_wmatrix_rectangular(dictio, symm = True):
  y_names = list(dictio.keys())
  x_names = list(dictio.get_hash_values_idx().keys())
  matrix = np.zeros((len(y_names), len(x_names)))
  i = 0
  for elementA, relations in dictio.items():
    for j, elementB in enumerate(x_names):
        query = relations.get(elementB)
        if query != None:
          matrix[i, j] = query
        elif symm: # TODO: PSZ, lo q se hace aqui no me cuadra
          query = dig(dictio, elementB, elementA)
          if query != None: matrix[i, j] = query
    i += 1
  return matrix, y_names, x_names

def dig(dictio, *keys):
  try:
    for key in keys:
        dictio = dictio[key]
    return dictio
  except KeyError:
    return None 

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

def write_profile_pairs(similarity_pairs, filename):
  with open(filename, 'w') as f:
    for pairsA, pairsB_and_values in similarity_pairs.items():
      for pairsB, values in pairsB_and_values.items():
        f.write(f"{pairsA}\t{pairsB}\t{values}\n")

def write_patient_hpo_stat(average_hp_per_pat_distribution, output_file):
  with open(output_file, 'w') as f:
    f.write("PatientsNumber\tHPOAverage\n")
    for patient_num, ave in average_hp_per_pat_distribution:
      f.write(f"{patient_num}\t{ave}\n")

def load_profiles(file_path, hpo):
  profiles = {}
  with open(file_path) as f:
    for line in f:
      id, profile = line.rstrip().split("\t")
      hpos = profile.split(',')
      hpos, rejected_hpos = hpo.check_ids(hpos)
      if len(hpos) > 0:
        hpos = hpo.clean_profile(hpos)
        if len(hpos) > 0 : profiles[id] = hpos
  return profiles

def invert_nested_hash(h):
  new_h = {}
  for k1, vals1 in h.items():
    for k2, vals2 in vals1: 
      query = new_h.get(k2)
      if query == None:
        new_h[k2] = {k1 : vals2}
      else:
        query[k1] = vals2
  return new_h

def remove_nested_entries(nested_hash, func):
  empty_root_ids = []
  for root_id, entries in nested_hash.items():
    delete_entries = []
    for k, v in entries:
      if not func(k,v): delete_entries.append(k)
    if len(delete_entries) == len(entries):
      empty_root_ids.append(root_id)
    else:
      for k in delete_entries: entries.pop(k)
  for k in empty_root_ids: nested_hash.pop(k)

def parse_clusters_file(clusters_file, patient_data):
  clusters_info = {}
  with open(clusters_file) as f:
    for line in f:
      patientID, clusterID = line.rstrip().split("\t")
      patientHPOProfile = patient_data.get_profile(patientID)
      query = clusters_info.get(clusterID)
      if query == None :
        clusters_info[clusterID] = {patientID: patientHPOProfile}
      else:
        query[patientID] = patientHPOProfile
  clusters_table = []
  for clusterID, patients_info in clusters_info.items():
    patients_per_cluster = len(patients_info)
    clusters_table.append([clusterID, patients_per_cluster, list(patients_info.keys()), list(patients_info.values())])
  return clusters_table, clusters_info

def get_cluster_metadata(clusters_info):
  average_hp_per_pat_distribution = []
  for cl_id, pat_info in clusters_info.items():
      hp_per_pat_in_clust = [ len(a) for a in pat_info.values() ]
      hp_per_pat_ave = sum(hp_per_pat_in_clust) / len(hp_per_pat_in_clust)
      average_hp_per_pat_distribution.append([len(pat_info), hp_per_pat_ave])
  return average_hp_per_pat_distribution

def translate_codes(clusters, hpo):
  translated_clusters = []
  for clusterID, num_of_pats, patientIDs_ary, patient_hpos_ary in clusters:
        translate_codes = [[ hpo.translate_id(code) for code in profile ] for profile in patient_hpos_ary ]
        translated_clusters.append([clusterID, 
          num_of_pats, 
          patientIDs_ary, 
          patient_hpos_ary, 
          translate_codes
        ])
  return translated_clusters

def get_semantic_similarity_clustering(options, patient_data, temp_folder, template_path, code_folder):
  template = open(template_path).read()
  hpo = Cohort.get_ontology(Cohort.act_ont)
  reference_profiles = None
  if options['reference_profiles'] != None: reference_profiles = load_profiles(options['reference_profiles'], hpo)
  for method_name in options['clustering_methods']:
    matrix_filename = os.path.join(temp_folder, f"similarity_matrix_{method_name}.npy")
    if reference_profiles == None:
      axis_file = re.sub('.npy','.lst', matrix_filename)
    else:
      axis_file_x = re.sub('.npy','_x.lst', matrix_filename)
      axis_file_y = re.sub('.npy','_y.lst', matrix_filename)
    profiles_similarity_filename = os.path.join(temp_folder, f'profiles_similarity_{method_name}.txt')
    clusters_distribution_filename = os.path.join(temp_folder, f'clusters_distribution_{method_name}.txt')

    if not os.path.exists(matrix_filename):
      if reference_profiles == None: 
        profiles_similarity = patient_data.compare_profiles(sim_type = method_name, external_profiles = reference_profiles)
      else: # AS reference profiles are constant, the sematic comparation will be A => B (A reference). So, we have to invert the elements to perform the comparation
        ont = Cohort.get_ontology('hpo')
        pat_profiles = ont.profiles
        ont.load_profiles(reference_profiles, reset_stored = True)
        profiles_similarity = ont.compare_profiles(sim_type = method_name, 
          external_profiles = pat_profiles, 
          bidirectional = False)
        ont.load_profiles(pat_profiles, reset_stored = True)
        profiles_similarity = invert_nested_hash(profiles_similarity)
      if options.get('sim_thr') != None: remove_nested_entries(profiles_similarity, lambda id, sim: sim >= options['sim_thr']) 
      write_profile_pairs(profiles_similarity, profiles_similarity_filename)
      if reference_profiles == None:
        similarity_matrix, axis_names = to_wmatrix(profiles_similarity, squared = True, symm = True)
        save(similarity_matrix, matrix_filename, axis_names, axis_file)
      else:
        similarity_matrix, y_names, x_names = to_wmatrix(profiles_similarity, squared = False, symm = True)
        save(similarity_matrix, matrix_filename, y_names, axis_file_y, x_names, axis_file_x)

    ext_var = ''
    if method_name == 'resnik':
      ext_var = '-m max'
    elif method_name == 'lin':
      ext_var = '-m comp1'
    cluster_file = f"{method_name}_clusters.txt"
    if reference_profiles != None:
      ext_var = ext_var + ' -s' 
      axis_file = f"{axis_file_y},{axis_file_x}"
      cluster_file = f"{method_name}_clusters_rows.txt"
    out_file = os.path.join(temp_folder, method_name)
    if not os.path.exists(out_file +  '_heatmap.png'): system_call(code_folder, 'plot_heatmap.R', f"-y {axis_file} -d {matrix_filename} -o {out_file} -M {options['minClusterProportion']} -t dynamic -H {ext_var}") 
    clusters_codes, clusters_info = parse_clusters_file(os.path.join(temp_folder, cluster_file), patient_data)  
    write_patient_hpo_stat(get_cluster_metadata(clusters_info), clusters_distribution_filename)
    out_file = os.path.join(temp_folder, f"{method_name}_clusters_distribution")
    if not os.path.exists(out_file): system_call(code_folder, 'xyplot_graph.R', f"-d {clusters_distribution_filename} -o {out_file} -x PatientsNumber -y HPOAverage") 
    sim_mat4cluster = {}
    if options['detailed_clusters']:
      for clID, patient_number, patient_ids, hpo_codes in clusters_codes:
        cluster_cohort = Cohort()
        for i, patID in enumerate(patient_ids): cluster_cohort.add_record([patID, hpo_codes[i], []])
        cluster_profiles = cluster_cohort.profiles
        ref_profile = cluster_cohort.get_general_profile()
        hpo.load_profiles({'ref': ref_profile}, reset_stored = True)    
        similarities = hpo.compare_profiles(external_profiles = cluster_profiles, sim_type = 'lin', bidirectional = False)
        candidate_sim_matrix, candidates, candidates_ids = get_similarity_matrix(ref_profile, similarities['ref'], cluster_profiles, hpo, 100, 100)
        candidate_sim_matrix.insert(0, ['HP'] + candidates_ids)
        sim_mat4cluster[clID] = candidate_sim_matrix

    clusters = translate_codes(clusters_codes, hpo)
    container = {
      'temp_folder' : temp_folder,
      'cluster_name' : method_name,
      'clusters' : clusters,
      'hpo' : hpo,
      'sim_mat4cluster' : sim_mat4cluster
     }

    report = Py_report_html(container, title='Patient clusters report')
    report.build(template)
    report.write(options['output_file']+ f"_{method_name}_clusters.html")
    if not os.path.exists(os.path.join(temp_folder, method_name + '_sim_boxplot.png')): system_call(code_folder, 'generate_boxpot.R', f"-i {temp_folder} -m {method_name} -o {os.path.join(temp_folder, method_name + '_sim_boxplot')}") 


def get_similarity_matrix(reference_prof, similarities, evidence_profiles, hpo, term_limit, candidate_limit, other_scores = {}, id2label = {}):
  candidates = [ list(pair) for pair in similarities.items()]
  if len(other_scores) == 0:
    candidates.sort(key=lambda s: s[-1], reverse=True)
    candidates = candidates[:candidate_limit]
  else: # Prioritize first by the external list of scores, select the candidates and then rioritize by similarities
    selected_candidates = []
    for cand in candidates:
      cand_id = cand[0]
      cand_lab = id2label.get(str(cand_id))
      if cand_lab == None: continue
      other_score = other_scores.get(cand_lab)
      if other_score == None: continue
      cand.append(other_score)
      selected_candidates.append(cand)
    selected_candidates.sort(key=lambda e: e[2], reverse=True)
    candidates = selected_candidates[:candidate_limit]
    candidates.sort(key=lambda e: e[1], reverse=True)
    for c in candidates:
      c.pop()

  candidates_ids = [c[0] for c in candidates]
  candidate_similarity_matrix = get_detailed_similarity(reference_prof, candidates, evidence_profiles, hpo)
  for i, row in enumerate(candidate_similarity_matrix):
    row.insert(0,hpo.translate_id(reference_prof[i]))

  candidate_similarity_matrix.sort(key=lambda r: sum(r[1:len(r)]), reverse=True)
  candidate_similarity_matrix = candidate_similarity_matrix[:term_limit]
  return candidate_similarity_matrix, candidates, candidates_ids

def get_detailed_similarity(profile, candidates, evidences, hpo):
  profile_length = len(profile)
  matrix = []
  for times in range(profile_length):
    matrix.append([0]*len(candidates))
  cand_number = 0
  for candidate_id, similarity in candidates:
    local_sim = []
    candidate_evidence = evidences[candidate_id]
    for profile_term in profile:
      for candidate_term in candidate_evidence:
        term_sim = hpo.compare([candidate_term], [profile_term], sim_type = "lin", bidirectional= False)
        local_sim.append([profile_term, candidate_term, term_sim])

    local_sim.sort(key = lambda s: s[-1], reverse=True)
    final_pairs = []
    processed_profile_terms = []
    processed_candidate_terms = []
    for pr_term, cd_term, sim in local_sim:
      if pr_term not in processed_profile_terms and cd_term not in processed_candidate_terms:
        final_pairs.append( [pr_term, cd_term, sim])
        processed_profile_terms.append( pr_term)
        processed_candidate_terms.append( cd_term)
      if profile_length == len(processed_profile_terms): break

    for pr_term, cd_term, similarity in final_pairs:
      matrix[profile.index(pr_term)][cand_number] = similarity
    cand_number += 1
  return matrix    