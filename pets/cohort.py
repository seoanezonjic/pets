import json
import semtools
import os, sys
from collections import defaultdict
from semtools.ontology import Ontology
from pets.genomic_features import Genomic_Feature

class Cohort():
    #TODO: Ask Pedro about these first lines of code
    ont = {}
    
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
        ont = None
        if ".json" not in ont_file:
            if excluded_terms_file:
                ont = Ontology(file= ont_file, load_file= True, removable_terms= Cohort.read_excluded_ont_file(excluded_terms_file))
            else:
                ont = Ontology(file= ont_file, load_file= True)
        else:
            ont = Ontology()
            ont.read(ont_file)
            if excluded_terms_file:
                ont.add_removable_terms(Cohort.read_excluded_ont_file(excluded_terms_file))
                ont.remove_removable()
                ont.build_index()
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

    def add_record(self, rec, extra_attr = None): #[id, [profile], [[chr1, start1, stop1],[chr1, start1, stop1]]]
        id, profile, vars = rec
        if profile: self.profiles[id] = profile
        if extra_attr: self.extra_attr[id] = extra_attr 
        if vars: self.add_gen_feat(id, vars)

    def delete(self, id):
        del(self.profiles[id])
        del(self.vars[id])

    def select_by_profile(self, func):
        self.profiles = dict( filter(lambda id, profile: func(id, profile), self.profiles.items()) )
        current_ids = self.profiles.keys()
        self.vars = dict( filter(lambda id, var: id in current_ids, self.vars.items() ) )

    def select_by_var(self, func):
        self.vars = dict( filter(lambda id, var: func(id, var), self.vars.items()) )
        current_ids = self.vars.keys()
        self.profiles = dict( filter(lambda id, profile: id in current_ids, self.profiles.items() ) )

    def filter_by_term_number(self, n_terms):
        self.select_by_profile(lambda id, profile, n_terms=n_terms: len(profile) >= n_terms)
    
    def remove_incomplete_records(self): # remove resc that lacks of vars or phenotypes
        ids_with_terms = self.profiles.keys()
        ids_with_vars = []
        for id, regs in self.vars.items():
            if len(regs) > 0: ids_with_vars.append(id)

        full_ids = self.intersection(ids_with_vars, ids_with_terms)

        self.profiles = dict( filter( 
            lambda id, prof: id in full_ids, 
            self.profiles.items() ) )

        self.vars = dict( filter(
            lambda id, var: id in full_ids,
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
                sys.stderr.write(f"WARNING: record {id} has the unknown CODES '{rejec_terms.join(',')}'. Codes removed.")
                rejected_terms.extend(rejec_terms)

            if not terms or len(terms) == 0:
                rejected_recs.append(id)
            else:
                self.profiles[id] = terms

        self.profiles = dict(filter(
            lambda id, record: id not in rejected_recs, 
            self.profiles))
        
        self.vars = dict(filter(
            lambda id, record: id not in rejected_recs,
            self.vars))

        return list(set(rejected_terms)), rejected_recs

    def link2ont(self, ont_id):
        Cohort.ont[ont_id].load_profiles(self.profiles)

    def get_profile_redundancy(self):
        ont = Cohort.ont[Cohort.act_ont]
        profile_sizes, parental_terms_per_profile = ont.get_profile_redundancy()
        return profile_sizes, parental_terms_per_profile

    def get_profiles_terms_frequency(self, options={}):
        ont = Cohort.ont[Cohort.act_ont]
        term_stats = ont.get_profiles_terms_frequency(**options) #https://www.ruby-lang.org/en/news/2019/12/12/separation-of-positional-and-keyword-arguments-in-ruby-3-0/
        return term_stats

    def compute_term_list_and_childs(self):
        ont = Cohort.ont[Cohort.act_ont]
        suggested_childs, term_with_childs_ratio = ont.compute_term_list_and_childs()
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

    def compare_profiles(self, options={}):
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
                if id_variants == None or len(id_variants) == 0:
                    variants.append(['-', '-', '-'])
                else:
                    for chrm, reg in id_variants.each():
                        variants.append([chrm, reg["start"], reg["stop"]])

                for var in variants:
                    vars_joined = "\t".join(var)
                    if mode == "default":
                        f.write(f"{id}\t{terms.join('|')}\t{vars_joined}")
                    elif mode == "paco":
                        f.write(f"{id}\t{vars_joined}\t{terms.join('|')}")
                    else:
                        raise Exception('Wrong save mode option, please try default or paco')

    def export_phenopackets(self, output_folder, genome_assembly, vcf_index: None):
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
            phenopacket = {metaData: metaData}
            query_sex = self.extra_attr.get(id).get("sex")
            sex = 'UNKNOWN_SEX' if query_sex == None else query_sex
            phenopacket["subject"] = {
                id: id,
                sex: sex
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
            id_variants = self.vars.get(id)
            variants = []
            if id_variants == None or len(id_variants) == 0:
                variants.append(['-', '-', '-'])
            else:
                for chrm, reg in id_variants.each():
                    variants.append([chrm, reg["start"], reg["stop"]])

    #Supplementary functions
    def intersection(self, arr1, arr2):
         return [item for item in arr1 if item in arr2] 