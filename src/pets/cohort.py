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
from pets.io import load_hpo_ci_values

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
    def load_ontology(cls, ont_name, ont_file, excluded_terms_file = None, inplace = True):
        if excluded_terms_file == None:
            excluded_terms = []
        else:
            excluded_terms = Cohort.read_excluded_ont_file(excluded_terms_file)
        ont = Ontology(file= ont_file, removable_terms= excluded_terms)
        if inplace: 
            cls.ont[ont_name] = ont
        else:
            return ont

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

        for rr in rejected_recs:
            self.profiles.pop(rr)
            if self.vars: self.vars.pop(rr)

        return list(set(rejected_terms)), rejected_recs

    def link2ont(self, ont_id):
        Cohort.ont[ont_id].load_profiles(self.profiles)

    def get_profile_redundancy(self):
        ont = Cohort.ont[Cohort.act_ont]
        profile_sizes, parental_terms_per_profile = ont.get_profile_redundancy()
        return profile_sizes, parental_terms_per_profile

    def compute_term_list_and_childs(self, file = None):
        ont = Cohort.ont[Cohort.act_ont]
        suggested_childs, term_with_childs_ratio = ont.compute_term_list_and_childs()
        if file != None: self.write_detailed_hpo_profile_evaluation(suggested_childs, file)
        return suggested_childs, term_with_childs_ratio

    def get_ic_analysis(self, freq_type = None, ic_file = None):
        ont = Cohort.ont[Cohort.act_ont]
        onto_ic, freq_ic = ont.get_observed_ics_by_onto_and_freq() # IC for TERMS
        if freq_type == 'freq_internal':
            freq_ic = load_hpo_ci_values(ic_file)
            ont.ics['resnik_observed'] = freq_ic
            for pat_id, phenotypes in self.each_profile():
                ont.dicts['prof_IC_observ'][pat_id] = self.get_profile_ic(phenotypes, freq_ic)
        ont.get_profiles_resnik_dual_ICs() # IC for PROFILES
        if freq_type == 'freq_internal' or freq_type == 'freq':
            phenotype_ic = freq_ic
        elif freq_type == 'onto':
            phenotype_ic = onto_ic
        return phenotype_ic

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
    def export_phenopackets(self, output_folder, genome_assembly, vcf_index= None, attr_index=None, attr_name = None, v2= False):
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
                phen = {"type": { "id": term, "label": term_name}}
                if not v2: phen["classOfOnset"] = {"id": "HP:0003577", "label": "Congenital onset"}
                phenotypicFeatures.append(phen)
            extra_attr = self.extra_attr.get(id)
            if extra_attr != None and extra_attr.get('neg_hpo'):
                neg_hpos = extra_attr.get("neg_hpo")
                for term in neg_hpos:
                    term_name = ont.translate_id(term)
                    phen = {"type": { "id": term, "label": term_name}, "excluded": True}
                    if not v2: phen["classOfOnset"] = {"id": "HP:0003577", "label": "Congenital onset"}
                    phenotypicFeatures.append(phen)

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
    
            #attr_name 
            if attr_index and id in attr_index:
                if attr_name == 'hgvsc':
                    interpretation = {
                        'id': id,
                        'diagnosis': {
                            'genomicInterpretations': [{
                                'subjectOrBiosampleId': id,
                                "interpretationStatus": "CAUSATIVE",
                                "variantInterpretation": {
                                    "variationDescriptor": {
                                        "expressions": [
                                            {"syntax": "hgvs.c", "value": attr_index[id]}
                                        ],
                                        "moleculeContext": "genomic",
                                        "allelicState": { # needed for GPSEA
                                          "id": "GENO:0000135",
                                          "label": "heterozygous"
                                        }
                                    }
                                }
                            }]
                        }
                    }

                    phenopacket["interpretations"] = [interpretation]

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
