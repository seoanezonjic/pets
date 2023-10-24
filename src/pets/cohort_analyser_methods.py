import numpy as np
from collections import defaultdict

import pylab
from py_exp_calc import exp_calc
from py_report_html import Py_report_html
from pets.cohort import Cohort
from pets.parsers.cohort_parser import Cohort_Parser

def get_arc_degree_and_radius_values(term, ontology, level_linspace, level_current_index):
  hp_level = ontology.get_term_level(term) - 2
  current_level_idx = level_current_index[hp_level]
  current_level_arc_array = level_linspace[hp_level]
  arc_hp_ont = float(current_level_arc_array[current_level_idx])
  
  level_current_index[hp_level] += 1
  return arc_hp_ont, hp_level

def append_values_to_arrays(arrays, values):
  for idx, value in enumerate(values):
    arrays[idx].append(value)

def prepare_rontoplot_data(hpo_stats_dict, ontology, root_node, reference_node):
  level_terms = ontology.get_ontology_levels()
  hps_to_filter_out = set()
  del level_terms[1] 
  for term in level_terms[2]:
    if term != root_node: hps_to_filter_out.update(ontology.get_descendants(term)+[term])

  cleaned_level_terms = {(level - 2): [term for term in terms if term not in hps_to_filter_out] for level, terms in level_terms.items()}
  level_linspace = {level: np.linspace(0, 2*np.pi, len(terms)) for level, terms in cleaned_level_terms.items()}
  level_current_index = {level: 0 for level in level_linspace.keys()}
  visited_terms = set(root_node)
  terms_to_visit = [term for term in ontology.get_direct_descendants(root_node) if term not in hps_to_filter_out]
  
  color_palette = Py_report_html.get_color_palette(len([ term for term in ontology.get_direct_descendants(reference_node) if term not in hps_to_filter_out]))
  top_parental_colors = {term: color_palette.pop() for term in ontology.get_direct_descendants(reference_node) if term not in hps_to_filter_out}
  grey = (128.0/256, 128.0/256 , 128.0/256, 1.0)
  colors, sizes, radius_values, arc_values = [grey], [1], [0], [0]
  while len(terms_to_visit) > 0:
    term = terms_to_visit.pop(0)
    if term in visited_terms: continue
    visited_terms.add(term)    
    childs = ontology.get_direct_descendants(term)
    if childs != None and len(childs) > 0: terms_to_visit = [term for term in childs if term not in hps_to_filter_out] + terms_to_visit

    arc_hp_ont, hp_level = get_arc_degree_and_radius_values(term, ontology, level_linspace, level_current_index)
    if term in top_parental_colors.keys(): current_color = top_parental_colors[term]
    append_values_to_arrays([colors, sizes, radius_values, arc_values], [grey, 1, hp_level, arc_hp_ont])
    if hpo_stats_dict.get(term) != None: append_values_to_arrays([colors, sizes, radius_values, arc_values], [current_color, 1 + hpo_stats_dict[term], hp_level + 0.3, arc_hp_ont])
  return [colors, sizes, radius_values, arc_values]

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
    exp_calc.add_record(raw_coverage, chrm, [start, stop, number_of_patients])
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
  for chrm, coords_info in raw_coverage.items():
    for start, stop, pat_records in coords_info:
      region_length = stop - start + 1
      for i in range(pat_records):
        all_cnvs_length.append(region_length)
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
def get_cluster_chromosome_data(cluster_data, limit):
  index = 0
  cl_ids = []
  data = {}
  if len(cluster_data) > 0: last_id = cluster_data[0][0] 
  for cluster_id, patient_number, chrm, count in cluster_data:
    if cluster_id != last_id: index += 1 
    if index == limit: break 
    cl_id = f"{patient_number}_{index}"
    if not cl_id in cl_ids: cl_ids.append(cl_id)
    query_chrm = data.get(chrm)
    if query_chrm == None:
      data[chrm] = {cl_id: count}
    else:
      query_chrm[cl_id] = count
    last_id = cluster_id

  cl_chr_data = [['cls_id'] + cl_ids]
  for chrm, ids_counts in data.items():
    chr_counts = [chrm]
    for cl_id in cl_ids:
      count = ids_counts.get(cl_id)
      if count == None:
        chr_counts.append(0)
      else:
        chr_counts.append(count)
    cl_chr_data.append(chr_counts)
  return cl_chr_data

def format_cluster_ic_data(all_ics, profile_lengths, limit):
  ic_data = [['cluster_id', 'ic', 'Plen']]
  for i, cluster_ics in enumerate(all_ics):
    if i == limit: break
    cluster_length = len(cluster_ics)
    for j, clust_ic in enumerate(cluster_ics): ic_data.append([f"{cluster_length}_{i}", clust_ic, profile_lengths[i][j]])
  return ic_data

def parse_clusters_data(patient_clusters, patient_data):
  clusters_info = {}
  for clusterID, patientIDs, in patient_clusters.items():
    for patientID in patientIDs:
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

def get_similarities4boxplot(raw_cls, similarity_matrix):
    sim_table = [['Sims', 'group'] ]
    if raw_cls is not None:
      cl= {}
      for i, item in enumerate(raw_cls): 
        cl_id = item[0]
        query = cl.get(cl_id)
        if query is None:
          cl[cl_id] = [i]
        else:
          query.append(i)
      cl_similarities = []
      for c_id, idxs in cl.items():
        np_ids = np.array(idxs)
        submatrix = similarity_matrix[np_ids[:,None], np_ids[None,:]]
        cl_similarities.extend(submatrix.reshape(submatrix.size).tolist())
      all_similarities = similarity_matrix.reshape(similarity_matrix.size).tolist()
      sim_table.extend([s , 'all'] for s in all_similarities)
      sim_table.extend([ [s, 'cls'] for s in cl_similarities] )
    return sim_table

def get_semantic_similarity_clustering(options, patient_data, reference_profiles, temp_folder, template_path_obj):
  template = open(template_path_obj).read()
  hpo = Cohort.get_ontology(Cohort.act_ont)
  clustering_data = {}
  for method_name in options['clustering_methods']:
    clusters, similarity_matrix, linkage, raw_cls = patient_data.get_similarity_clusters(method_name, 'hpo', options, temp_folder = temp_folder, reference_profiles = reference_profiles)
    clusters_codes, clusters_info = parse_clusters_data(clusters, patient_data)
    clustering_data[method_name] = {'cls': clusters, 'sim': similarity_matrix, 'link': linkage, 'raw_cls': raw_cls, 'boxplot_sims': get_similarities4boxplot(raw_cls, similarity_matrix)}

    sim_mat4cluster = {}
    if options['detailed_clusters']:
      for clID, patient_number, patient_ids, hpo_codes in clusters_codes:
        cluster_cohort = Cohort()
        for i, patID in enumerate(patient_ids): cluster_cohort.add_record([patID, hpo_codes[i], []])
        cluster_profiles = cluster_cohort.profiles
        ref_profile = cluster_cohort.get_general_profile()
        hpo.load_profiles({'ref': ref_profile}, reset_stored = True)    
        candidate_sim_matrix, candidates, candidates_ids = cluster_cohort.calc_sim_term2term_similarity_matrix(ref_profile, 'ref', cluster_profiles, hpo, 
          term_limit = 100, candidate_limit = 100, sim_type = 'lin', bidirectional = False)
        candidate_sim_matrix.insert(0, ['HP'] + candidates_ids)
        sim_mat4cluster[clID] = candidate_sim_matrix

    container = {
      'temp_folder' : temp_folder,
      'cluster_name' : method_name,
      'clusters' : translate_codes(clusters_codes, hpo),
      'hpo' : hpo,
      'sim_mat4cluster' : sim_mat4cluster
     }

    report = Py_report_html(container, title='Patient clusters report')
    report.build(template)
    report.write(options['output_file']+ f"_{method_name}_clusters.html")
    
  return clustering_data