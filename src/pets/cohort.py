import json
import os, sys
import re
import numpy as np
from collections import defaultdict
import csv
from py_exp_calc import exp_calc
from py_exp_calc.exp_calc import intersection
from py_semtools.ontology import Ontology
from pets.genomic_features import Genomic_Feature

class Cohort():
    #TODO: Ask Pedro about these first lines of code
    ont = {}
    act_ont = "hpo"
    profiles = {}
    """
    class << self # https://www.ruby-forum.com/t/attr-accessor-for-class-variable/136693
        attr_accessor :act_ont # Which ontology use for ont related operations
    end

    attr_accessor :profiles
    """

    @classmethod
    def get_ontology(cls, ont_id):
        return cls.ont[ont_id]

    @classmethod
    def load_ontology(cls, ont_name, ont_file, excluded_terms_file = None):
        if excluded_terms_file == None:
            excluded_terms = []
        else:
            excluded_terms = Cohort.read_excluded_ont_file(excluded_terms_file)
        ont = Ontology(file= ont_file, removable_terms= excluded_terms) 
        cls.ont[ont_name] = ont

    @classmethod
    def read_excluded_ont_file(cls, file):
        excluded_hpo = []
        with open(file) as f:
            for line in f:
                excluded_hpo.append(line.strip())
        return excluded_hpo

    def __init__(self):
        self.profiles = {}
        self.vars = {}
        self.extra_attr = {}
        self.var_idx = Genomic_Feature([])

    def add_record(self, rec, extra_attr = None): #rec= [id, [profile], [[chr1, start1, stop1],[chr1, start1, stop1]]]
        id, profile, vars = rec
        if profile: self.profiles[id] = profile
        if extra_attr: self.extra_attr[id] = extra_attr 
        if vars: self.add_gen_feat(id, vars)

    def delete(self, id):
        del(self.profiles[id])
        del(self.vars[id])

    def select_by_profile(self, func):
        self.profiles = dict( filter(lambda id_profile_pair: func(id_profile_pair[0], id_profile_pair[1]), self.profiles.items()) )
        current_ids = self.profiles.keys()
        self.vars = dict( filter(lambda id_var_pair: id_var_pair[0] in current_ids, self.vars.items() ) )

    def select_by_var(self, func):
        self.vars = dict( filter(lambda id_var_pair: func(id_var_pair[0], id_var_pair[1]), self.vars.items()) )
        current_ids = self.vars.keys()
        self.profiles = dict( filter(lambda id_profile_pair: id_profile_pair[0] in current_ids, self.profiles.items() ) )

    def filter_by_term_number(self, n_terms):
        self.select_by_profile(lambda id, profile, n_terms=n_terms: len(profile) >= n_terms)
    
    def remove_incomplete_records(self): # remove resc that lacks of vars or phenotypes
        ids_with_terms = []
        for id, profile in self.profiles.items():
            if len(profile) > 0: ids_with_terms.append(id)
        ids_with_vars = []
        for id, regs in self.vars.items():
            if regs.len() > 0: ids_with_vars.append(id)

        full_ids = intersection(ids_with_vars, ids_with_terms)

        self.profiles = dict( filter( 
            lambda id_prof_pair: id_prof_pair[0] in full_ids, 
            self.profiles.items() ) )

        self.vars = dict( filter(
            lambda id_var_pair: id_var_pair[0] in full_ids,
            self.vars.items() ) )
        

    def add_gen_feat(self, id, feat_array): # [[chr1, start1, stop1],[chr1, start1, stop1]]
        self.vars[id] = Genomic_Feature(feat_array)

    def get_profile(self, id):
        return self.profiles[id]

    def get_vars(self, id):
        return self.vars[id]

    def each_profile(self):
        for id, profile in self.profiles.items():
            yield(id, profile)

    def each_var(self):
        for id, var_info in self.vars.items():
            yield(id, var_info)

    def get_general_profile(self,thr=0): # TODO move funcionality to semtools
        term_count = defaultdict(lambda: 0)
        for id, prof in self.each_profile():
            for term in prof:
                term_count[term] += 1

        records = len(self.profiles)
        general_profile = []
        for term, count in term_count.items():
            if count / float(records) >= thr: general_profile.append(term)
        
        #TODO: check if this is the correct way to access the ontology (migrated from ruby)
        ont = Cohort.ont[Cohort.act_ont]
        # ont = @@ont[Cohort.act_ont]
        return ont.clean_profile_hard(general_profile)

    def check(self, hard=False): # OLD format_patient_data
        #TODO: check if this is the correct way to access the ontology (migrated from ruby)
        ont = Cohort.ont[Cohort.act_ont]
        #ont = @@ont[Cohort.act_ont]
        rejected_terms = []
        rejected_recs = []
        for id, terms in self.profiles.items():
            if hard:
                terms = ont.clean_profile_hard(terms)
                rejec_terms = []
            else:
                terms, rejec_terms = ont.check_ids(terms)

            if rejec_terms and len(rejec_terms) > 0:
                sys.stderr.write(f"WARNING: record {id} has the unknown CODES '{','.join(rejec_terms)}'. Codes removed.")
                rejected_terms.extend(rejec_terms)

            if not terms or len(terms) == 0:
                rejected_recs.append(id)
            else:
                self.profiles[id] = terms

        self.profiles = dict(filter(
            lambda id_record_pair: id_record_pair[0] not in rejected_recs, 
            self.profiles.items()))
        
        self.vars = dict(filter(
            lambda id_record_pair: id_record_pair[0] not in rejected_recs,
            self.vars.items()))

        return list(set(rejected_terms)), rejected_recs

    def link2ont(self, ont_id):
        Cohort.ont[ont_id].load_profiles(self.profiles)

    def get_profile_redundancy(self):
        ont = Cohort.ont[Cohort.act_ont]
        profile_sizes, parental_terms_per_profile = ont.get_profile_redundancy()
        return profile_sizes, parental_terms_per_profile

    def get_profiles_terms_frequency(self, **options):
        ont = Cohort.ont[Cohort.act_ont]
        term_stats = ont.get_profiles_terms_frequency(**options) #https://www.ruby-lang.org/en/news/2019/12/12/separation-of-positional-and-keyword-arguments-in-ruby-3-0/
        return term_stats

    def compute_term_list_and_childs(self, file = None):
        ont = Cohort.ont[Cohort.act_ont]
        suggested_childs, term_with_childs_ratio = ont.compute_term_list_and_childs()
        if file != None: self.write_detailed_hpo_profile_evaluation(suggested_childs, file)
        return suggested_childs, term_with_childs_ratio

    def get_profile_ontology_distribution_tables(self):
        ont = Cohort.ont[Cohort.act_ont]
        ontology_levels, distribution_percentage = ont.get_profile_ontology_distribution_tables()
        ontology_levels.insert(0, ["level", "ontology", "cohort"])
        distribution_percentage.insert(0, ["level", "ontology", "weighted cohort", "uniq terms cohort"])
        return ontology_levels, distribution_percentage

    def get_ic_analysis(self):
        ont = Cohort.ont[Cohort.act_ont]
        onto_ic, freq_ic = ont.get_observed_ics_by_onto_and_freq() # IC for TERMS
        onto_ic_profile, freq_ic_profile = ont.get_profiles_resnik_dual_ICs() # IC for PROFILES
        return onto_ic, freq_ic, onto_ic_profile, freq_ic_profile

    def get_profiles_mean_size(self):
        ont = Cohort.ont[Cohort.act_ont]
        profile_mean_size = ont.get_profiles_mean_size()
        return profile_mean_size

    def get_profile_length_at_percentile(self, perc=50, increasing_sort=False):
        ont = Cohort.ont[Cohort.act_ont]
        length_percent = ont.get_profile_length_at_percentile(perc=perc, increasing_sort=increasing_sort)
        return length_percent

    def get_dataset_specifity_index(self, type):
        ont = Cohort.ont[Cohort.act_ont]
        dsi = ont.get_dataset_specifity_index(type)
        return dsi

    def compare_profiles(self, **options):
        ont = Cohort.ont[Cohort.act_ont]
        similarities = ont.compare_profiles(**options)
        return similarities

    def index_vars(self): # equivalent to process_patient_data
        for id, var in self.each_var():
            self.var_idx.merge(var, id)

    def get_vars_sizes(self, summary=False):
        if summary:
            return self.var_idx.get_summary_sizes()
        else:
            return self.var_idx.get_sizes()

    def generate_cluster_regions(self, meth, tag, lim):
        ids_by_cluster, annotated_full_ref = self.var_idx.generate_cluster_regions(meth, tag, lim)
        return ids_by_cluster, annotated_full_ref

    def save(self, output_file, mode = "default", translate = False):
        with open(output_file, "w") as f:
            if mode == 'paco': f.write("id\tchr\tstart\tstop\tterms\n")
            ont = Cohort.ont[Cohort.act_ont]
            for id, terms in self.profiles.items():
                if translate: terms, rejected = ont.translate_ids(terms)
                id_variants = self.vars.get(id)
                variants = []
                if id_variants == None or id_variants.len() == 0:
                    variants.append(['-', '-', '-'])
                else:
                    for chrm, reg in id_variants.each():
                        variants.append([chrm, reg["start"], reg["stop"]])

                for var in variants:
                    vars_joined = "\t".join([str(item) for item in var])
                    if mode == "default":
                        f.write(f"{id}\t{'|'.join(sorted(terms))}\t{vars_joined}\n")
                    elif mode == "paco":
                        f.write(f"{id}\t{vars_joined}\t{'|'.join(sorted(terms))}\n")
                    else:
                        raise Exception('Wrong save mode option, please try default or paco')
    
    #TODO: test the method
    def export_phenopackets(self, output_folder, genome_assembly, vcf_index= None):
        ont = Cohort.ont[Cohort.act_ont]
        metaData = {
            "createdBy": "PETS",
            "resources": [{
                "id": "hp",
                "name": "human phenotype ontology",
                "namespacePrefix": "HP",
                "url": "http://purl.obolibrary.org/obo/hp.owl",
                #"version" => "2018-03-08",
                "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
            }]
        }

        for id, terms in self.profiles.items():
            phenopacket = {"metaData": metaData}
            query_sex = self.extra_attr.get(id)
            query_sex = query_sex.get("sex") if query_sex != None else None
            sex = 'UNKNOWN_SEX' if query_sex == None else query_sex
            phenopacket["subject"] = {
                "id": id,
                "sex": sex
            }
            phenotypicFeatures = []
            for term in terms:
                term_name = ont.translate_id(term)
                phenotypicFeatures.append({
                    "type": { "id": term, "label": term_name},
                    "classOfOnset": {"id": "HP:0003577", "label": "Congenital onset"}
                })

            phenopacket["phenotypicFeatures"] = phenotypicFeatures
            if vcf_index and id in vcf_index:
                htsFiles = []
                htsFiles.append({
                    "uri": "file:/" + vcf_index[id],
                    "description": id,
                    "htsFormat": "VCF",
                    "genomeAssembly": genome_assembly,
                    "individualToSampleIdentifiers": { "patient1": id }
                })
                phenopacket["htsFiles"] = htsFiles

            with open(os.path.join(output_folder, str(id) + ".json"), "w") as f:
                f.write(json.dumps(phenopacket, indent=4))

    def process_dummy_clustered_patients(self, options, phenotype_ic, temp_folder = './'): # get ic and chromosomes
        if len(self.profiles) > 1 : clustered_patients = self.dummy_cluster_patients(temp_folder = temp_folder)
        ont = self.get_ontology(Cohort.act_ont)
        all_ics = []
        all_lengths = []
        top_cluster_phenotypes = []
        cluster_data_by_chromosomes = []
        multi_chr_clusters = 0
        processed_clusters = 0
        if len(self.profiles) > 1 :
            for cluster_id, patient_ids in sorted(list(clustered_patients.items()), key=lambda x: len(x[1]), reverse=True):
                num_of_patients = len(patient_ids)
                if num_of_patients == 1: continue 
                chrs, all_phens, profile_ics, profile_lengths = self.process_cluster(patient_ids, phenotype_ic, options, ont, processed_clusters)
                if processed_clusters < options['clusters2show_detailed_phen_data']: top_cluster_phenotypes.append(all_phens)
                all_ics.append(profile_ics)
                all_lengths.append(profile_lengths)
                if not options.get('chromosome_col') == None:
                    if len(chrs) > 1: multi_chr_clusters += num_of_patients
                    for chrm, count in chrs.items(): cluster_data_by_chromosomes.append([cluster_id, num_of_patients, chrm, count])
                processed_clusters += 1
        return all_ics, all_lengths, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chr_clusters

    def dummy_cluster_patients(self, temp_folder = "./"):
        clust_pat_file = os.path.join(temp_folder, 'cluster_asignation')
        matrix_file = os.path.join(temp_folder, 'pat_hpo_matrix.npy')
        if not os.path.exists(clust_pat_file):
            x_axis_file = re.sub('.npy','_x.lst', matrix_file)
            y_axis_file = re.sub('.npy','_y.lst', matrix_file)
            if not os.path.exists(matrix_file):
                pat_hpo_matrix, pat_id, hp_id  = exp_calc.to_bmatrix(self.profiles)
                exp_calc.save(pat_hpo_matrix, matrix_file, hp_id, x_axis_file, pat_id, y_axis_file)
            else:
                pat_hpo_matrix, hp_id, pat_id = exp_calc.load(matrix_file, x_axis_file, y_axis_file)
            clustered_patients, _ = exp_calc.get_hc_clusters(pat_hpo_matrix, dist = 'euclidean', method = 'ward', height = [1.5], item_list=pat_id)
            with open(clust_pat_file, 'w') as f:
                for clusterID, pat_ids in clustered_patients.items(): 
                    f.write(f"{clusterID}\t{','.join(pat_ids)}\n")
        else:
            clustered_patients = {}
            with open(clust_pat_file) as f:
                for line in f:
                    clusterID, pat_ids = line.rstrip().split('\t')
                    clustered_patients[int(clusterID)] = pat_ids.split(',')
        return(clustered_patients)



    def process_cluster(self, patient_ids, phenotype_ic, options, ont, processed_clusters):
      chrs = defaultdict(lambda: 0)
      all_phens = []
      profile_ics = []
      profile_lengths = []
      for pat_id in patient_ids:
        phenotypes = self.get_profile(pat_id) 
        profile_ics.append(self.get_profile_ic(phenotypes, phenotype_ic))
        profile_lengths.append(len(phenotypes))
        if processed_clusters < options['clusters2show_detailed_phen_data']:
          phen_names, rejected_codes = ont.translate_ids(phenotypes) #optional
          all_phens.append(phen_names)
        if not options.get('chromosome_col') == None: 
          for chrm in self.get_vars(pat_id).get_chr(): chrs[chrm] += 1
      return chrs, all_phens, profile_ics, profile_lengths 

    def get_profile_ic(self, hpo_names, phenotype_ic):
        ic = 0
        profile_length = 0
        for hpo_id in hpo_names:
            hpo_ic = phenotype_ic.get(hpo_id)
            if hpo_ic == None: raise Exception(f"The term {hpo_id} not exists in the given ic table")
            ic += hpo_ic 
            profile_length += 1
        if profile_length == 0: profile_length = 1 
        return ic/profile_length

    def get_matrix_similarity(self, method_name, options, ontology = 'hpo', reference_profiles=None, profiles_similarity_filename=None, matrix_filename = None):
        if reference_profiles == None: 
            profiles_similarity = self.compare_profiles(sim_type = method_name, external_profiles = reference_profiles)
        else: # AS reference profiles are constant, the sematic comparation will be A => B (A reference). So, we have to invert the elements to perform the comparation
            ont = Cohort.get_ontology(ontology)
            pat_profiles = ont.profiles # TEmporal copy to preserve patient profiles and inject reference profiles
            ont.load_profiles(reference_profiles, reset_stored = True)
            profiles_similarity = ont.compare_profiles(sim_type = method_name, 
                external_profiles = pat_profiles, 
                bidirectional = False)
            ont.load_profiles(pat_profiles, reset_stored = True)
            profiles_similarity = exp_calc.invert_nested_hash(profiles_similarity)
        if options.get('sim_thr') != None: exp_calc.remove_nested_entries(profiles_similarity, lambda id, sim: sim >= options['sim_thr']) 
        if profiles_similarity_filename != None: self.write_profile_pairs(profiles_similarity, profiles_similarity_filename)
        axis_file_x = re.sub('.npy','_x.lst', matrix_filename)
        axis_file_y = re.sub('.npy','_y.lst', matrix_filename)
        if reference_profiles == None:
            y_names = None
            similarity_matrix, x_names = exp_calc.to_wmatrix(profiles_similarity, squared = True, symm = True)
        else:
            similarity_matrix, y_names, x_names = exp_calc.to_wmatrix(profiles_similarity, squared = False, symm = False)
        exp_calc.save(similarity_matrix, matrix_filename, 
            x_axis_names=x_names, x_axis_file=axis_file_x, 
            y_axis_names=y_names, y_axis_file=axis_file_y)
        return similarity_matrix, y_names, x_names

    def get_similarity_clusters(self, method_name, ontology, options, temp_folder = None, reference_profiles = None):
        clusters = {}
        similarity_matrix = None
        linkage = None
        raw_cls = None
        if len(self.profiles) > 1:
            if temp_folder != None: # To save and load results from disk
                matrix_filename = os.path.join(temp_folder, f"similarity_matrix_{method_name}.npy")
                axis_file = re.sub('.npy','_x.lst', matrix_filename)
                axis_file_y = None if reference_profiles == None else re.sub('.npy','_y.lst', matrix_filename)
                profiles_similarity_filename = os.path.join(temp_folder, f'profiles_similarity_{method_name}.txt')
                cluster_file = os.path.join(temp_folder, f"{method_name}_clusters.txt")
                linkage_file = os.path.join(temp_folder, f"{method_name}_linkage.npy")
                raw_cls_file = os.path.join(temp_folder, f"{method_name}_raw_cls.npy")
            if temp_folder == None or not os.path.exists(matrix_filename):
                similarity_matrix, y_names, x_names = self.get_matrix_similarity(method_name, options, 
                    ontology = ontology, reference_profiles=reference_profiles,  
                    profiles_similarity_filename=profiles_similarity_filename, 
                    matrix_filename = matrix_filename)
            elif temp_folder != None or os.path.exists(matrix_filename):
                similarity_matrix, x_names, y_names = exp_calc.load(matrix_filename, x_axis_file=axis_file, y_axis_file=axis_file_y)
            
            if temp_folder == None or not os.path.exists(cluster_file):
                if method_name == 'resnik':
                    dist_matrix = np.amax(similarity_matrix) - similarity_matrix
                elif method_name == 'lin':
                    dist_matrix = 1 - similarity_matrix
                clusters, cls_objects = exp_calc.get_hc_clusters(dist_matrix, dist = 'custom', method = 'ward', identify_clusters='max_avg', n_clusters=3, item_list = x_names)
                linkage = cls_objects['link']
                raw_cls = cls_objects['cls']
                
                if temp_folder != None:
                    with open(cluster_file, 'w') as f:
                        for clusterID, patientIDs in clusters.items(): f.write(f"{clusterID}\t{','.join(patientIDs)}\n")
                    np.save(linkage_file, linkage)
                    np.save(raw_cls_file, np.array(raw_cls, dtype=np.int32))
                    #with open(raw_cls_file, 'w') as f: f.write(json.dumps(raw_cls))
            elif temp_folder != None or os.path.exists(cluster_file):
                with open(cluster_file) as f:
                    for l in f: 
                        clusterID, patientIDs = l.rstrip().split("\t")
                        clusters[int(clusterID)] = patientIDs.split(",")
                linkage = np.load(linkage_file)
                raw_cls = np.load(raw_cls_file)
                #with open(raw_cls_file) as f: raw_cls = json.loads(f.read())
        return clusters, similarity_matrix, linkage, raw_cls

    def calc_sim_term2term_similarity_matrix(self, ref_profile, ref_profile_id, external_profiles, ontology, term_limit = 100, candidate_limit = 100, sim_type = 'lin', bidirectional = True):
        similarities = ontology.compare_profiles(external_profiles = external_profiles, sim_type = sim_type, bidirectional = bidirectional)
        candidate_sim_matrix, candidates, candidates_ids = self.get_term2term_similarity_matrix(ref_profile, similarities[ref_profile_id], external_profiles, ontology, term_limit, candidate_limit)
        return candidate_sim_matrix, candidates, candidates_ids

    def get_term2term_similarity_matrix(self, reference_prof, similarities, evidence_profiles, hpo, term_limit, candidate_limit, other_scores = {}, id2label = {}):
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
            for c in candidates: c.pop()

        candidates_ids = [c[0] for c in candidates]
        candidate_similarity_matrix = self.get_detailed_similarity(reference_prof, candidates, evidence_profiles, hpo)
        for i, row in enumerate(candidate_similarity_matrix):
            row.insert(0,hpo.translate_id(reference_prof[i]))

        candidate_similarity_matrix.sort(key=lambda r: sum(r[1:len(r)]), reverse=True)
        candidate_similarity_matrix = candidate_similarity_matrix[:term_limit]
        return candidate_similarity_matrix, candidates, candidates_ids

    def get_detailed_similarity(self, profile, candidates, evidences, hpo):
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

    def write_detailed_hpo_profile_evaluation(self, suggested_childs, detailed_profile_evaluation_file):
        with open(detailed_profile_evaluation_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for pat_id, suggestions in suggested_childs.items():
                warning = ""
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

    def write_profile_pairs(self, similarity_pairs, filename):
        with open(filename, 'w') as f:
            for pairsA, pairsB_and_values in similarity_pairs.items():
                for pairsB, values in pairsB_and_values.items():
                    f.write(f"{pairsA}\t{pairsB}\t{values}\n")
