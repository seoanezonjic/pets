import sys, os, json, unittest 
from importlib.resources import files
import subprocess
from subprocess import PIPE
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort
from pets import HPO_FILE
import warnings
import numpy as np

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
PATIENTS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, '100_test_dataset_with_sex.txt'))


class CohortTestSuite(unittest.TestCase):
    def setUp(self):
        warnings.filterwarnings(action='ignore', category=UserWarning, module="py_semtools")

        patients = [
                    ["132", ["HP:0000717", "HP:0001252", "HP:0001249"],	[["22", 50, 100]], {"sex": "M"}],
                    ["599", ["HP:0001249"],	[["15", 1, 200]], {"sex": "F"}],
                    ["647", ["HP:0001249", "HP:0000262"], [["22", 500, 1000]], {"sex": "M"}],
                    ["648", ["HP:0000365", "HP:0000252", "HP:0001249"], [["22", 52, 102], ["22", 10000, 20000]], {"sex": "F"}]
                    ]
        self.n_patients = len(patients)

        self.expected_y_names = ['132', '599', '647', '648']
        self.expected_x_names = ['HP:0000717', 'HP:0001252', 'HP:0001249', 'HP:0000262', 'HP:0000365', 'HP:0000252']
        self.expected_matrix = [[1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 0],
                                [0, 0, 1, 1, 0, 0],
                                [0, 0, 1, 0, 1, 1]]
        
        self.ic_hpos = {"HP:0000717": 1.0876789846883463, "HP:0001252": 1.3467739354833985,
                        "HP:0001249": 0.4615785999681353, "HP:0000262": 3.2505419780102724,
                        "HP:0000365": 1.8888141419926796, "HP:0000252": 1.1510339842823074}

        hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
        Cohort.load_ontology("hpo", hpo_file)
        self.patient_data = Cohort()

        for patient in patients:
            id, terms, variants, other_attr = patient
            self.patient_data.add_record([id, terms, variants], other_attr)

        self.patient_data.link2ont(Cohort.act_ont)

    def test_load(self):
        self.assertEqual(len(self.patient_data.profiles), 4)
        self.assertEqual(len(self.patient_data.vars), 4)

        self.assertEqual(self.patient_data.profiles.get("132"), ["HP:0000717", "HP:0001252", "HP:0001249"])
        self.assertEqual(self.patient_data.vars.get("132").regions.get("22"), [{'chrm': '22', 'start': 50, 'stop': 100, 'to': 0}])
        self.assertEqual(self.patient_data.extra_attr["132"]["sex"], "M")


    def test_add_record(self):
        new_record = [2000, ["HP:0025500"], [["21", 20, 25], ["X", 1000, 2000]]]
        extra_attr = {"sex": "M"}

        self.patient_data.add_record(new_record, extra_attr)

        self.assertEqual(self.patient_data.profiles.get(2000),  ["HP:0025500"])
        self.assertEqual(self.patient_data.vars.get(2000).regions.get("21"), [{'chrm': '21', 'start': 20, 'stop': 25, 'to': 0}])
        self.assertEqual(self.patient_data.vars.get(2000).regions.get("X"), [{'chrm': 'X', 'start': 1000, 'stop': 2000, 'to': 1}])
        self.assertEqual(self.patient_data.extra_attr.get(2000), {"sex": "M"})


    def test_delete(self):
        new_record = [2000, ["HP:0025500"], [["21", 20, 25], ["X", 1000, 2000]]]
        extra_attr = {"sex": "M"}

        self.patient_data.add_record(new_record, extra_attr)
        self.assertTrue(self.patient_data.profiles.get(2000))
        self.assertTrue(self.patient_data.vars.get(2000))

        self.patient_data.delete(2000)
        self.assertEqual(self.patient_data.profiles.get(2000), None)
        self.assertEqual(self.patient_data.vars.get(2000), None)
        #TODO: Ask PSZ wether extra_attr should be removed as well


    def test_select_by_profile(self):
        def check_hp_in_pat(query):
            def _check_hp_in_pat(pat, hps):
                return query in hps
            return _check_hp_in_pat
            
        #Checking if the patients have the term HP:0000262
        self.patient_data.select_by_profile(check_hp_in_pat("HP:0000262")) #Only patient 647 has this term

        self.assertEqual(list(self.patient_data.profiles.keys()), ["647"])
        self.assertEqual(list(self.patient_data.vars.keys()), ["647"])
        #TODO: Ask PSZ wether extra_attr should be removed as well
        #self.assertEqual(len(self.patient_data.extra_attr.keys()), 4)

        for pat, hpterms in self.patient_data.profiles.items():
            self.assertIn("HP:0000262", hpterms)


    def test_select_by_var(self):
        def check_var_in_pat(query):
            def _check_var_in_pat(pat, vars):
                return query in vars.regions.keys()
            return _check_var_in_pat
        
        #Checking if the patients have a mutation in chromosome 15
        self.patient_data.select_by_var(check_var_in_pat("15")) #Only patient 599 has a mutation in this chromosome

        self.assertEqual(list(self.patient_data.profiles.keys()), ["599"])
        self.assertEqual(list(self.patient_data.vars.keys()), ["599"])

        for pat, vars in self.patient_data.vars.items():
            self.assertIn("15", vars.regions.keys())


    def test_filter_by_term_number(self):
        #We are filtering out patients with less than 3 phenotypes, so only patients 132 and 648 should remain
        self.patient_data.filter_by_term_number(3)

        self.assertEqual(list(self.patient_data.profiles.keys()), ["132", "648"])
        for pat, hpterms in self.patient_data.profiles.items():
            self.assertGreaterEqual(len(hpterms), 3)
        

    def test_remove_incomplete_records(self):
        #We are checking that we have the 4 patients in the toy dataset before the cleaning
        self.assertEqual(list(self.patient_data.profiles.keys()), ["132", "599", "647", "648"])
        self.assertEqual(list(self.patient_data.vars.keys()), ["132", "599", "647", "648"])

        self.patient_data.profiles["132"] = [] #making a patient to have no phenotype
        self.patient_data.vars["648"] = Genomic_Feature([]) #making another patient to have no genomic regions
        
        #We are filtering out patients without phenotype or genomic regions
        self.patient_data.remove_incomplete_records()

        #Now there should be 2 patients, 599 and 647
        self.assertEqual(list(self.patient_data.profiles.keys()), ["599", "647"])
        self.assertEqual(list(self.patient_data.vars.keys()), ["599", "647"])



    def test_add_gen_feat(self):
        self.assertEqual(list(self.patient_data.vars.get("132").regions.keys()), ["22"]) #Patient 132 has only one region in chromosome 22

        patient_new_gen_feat = [["21", 20, 25, 0], ["22", 49932021, 51187844]]
        self.patient_data.add_gen_feat("132", patient_new_gen_feat)
        self.assertEqual(list(self.patient_data.vars.get("132").regions.keys()), ["21", "22"])
        self.assertEqual(self.patient_data.vars.get("132").regions.get("21"), [{'chrm': '21', 'start': 20, 'stop': 25, 'to': 0}])
        self.assertEqual(self.patient_data.vars.get("132").regions.get("22"), [{'chrm': '22', 'start': 49932021, 'stop': 51187844, 'to': 1}])


    def test_get_profile(self):
        self.assertEqual(self.patient_data.get_profile("132"), ["HP:0000717", "HP:0001252", "HP:0001249"])


    def test_get_vars(self):
        self.assertEqual(self.patient_data.get_vars("132"), Genomic_Feature([["22", 50, 100]]))


    def test_each_profile(self):
        for pat, hpterms in self.patient_data.each_profile():
            self.assertIn(pat, ["132", "599", "647", "648"])
            self.assertIn(hpterms, self.patient_data.profiles.values())


    def test_each_var(self):
        for pat, vars in self.patient_data.each_var():
            self.assertIn(pat, ["132", "599", "647", "648"])
            self.assertIn(vars, self.patient_data.vars.values())


    def test_get_general_profile(self):
        #Intellectual disability HP:0001249 is the term with a frequency higher than 0.5 in the toy dataset
        general_profiles = self.patient_data.get_general_profile(thr=0.5)
        self.assertEqual(general_profiles, ["HP:0001249"])


    def test_check(self):
        #Adding a new record with an incorrect HP term
        new_record_bad = [2000, ["HP:BADCODE"], [["21", 20, 25], ["X", 1000, 2000]]]
        
        #Adding a new record with an obsolete HP term
        new_record_obs = [2001, ["HP:0000284"], [["21", 300, 2300], ["X", 3000, 4000]]] 

        #adding new records to the profiles
        self.patient_data.profiles[new_record_obs[0]] = new_record_obs[1]
        self.patient_data.profiles[new_record_bad[0]] = new_record_bad[1]

        #adding new records to genomic regions
        self.patient_data.add_gen_feat(2000, new_record_bad[2])
        self.patient_data.add_gen_feat(2001, new_record_obs[2])
        
        self.assertEqual(len(self.patient_data.profiles.items()), 6)
        self.assertEqual(len(self.patient_data.vars.items()), 6)

        #Checking the dataset. The incorrect term should be rejected (and also the patient, because s/he only has one HP term), 
        # but not the obsolete one, that will be replaced by its current equivalent
        rejected_terms, rejected_recs = self.patient_data.check(hard=False)

        self.assertEqual(rejected_recs, [2000])
        self.assertEqual(rejected_terms, ["HP:BADCODE"])

        self.assertEqual(len(self.patient_data.profiles.items()), 5)
        self.assertEqual(len(self.patient_data.vars.items()), 5)
        self.assertEqual(self.patient_data.profiles[2001], ["HP:0000315"])


    def test_link2ont(self):
        self.patient_data.link2ont(Cohort.act_ont)
        #We check that the profiles of the patients are correctly saved as profiles in the ontology object
        self.assertEqual(Cohort.ont[Cohort.act_ont].profiles.items(),  self.patient_data.profiles.items())


    def test_get_profile_redundancy(self):
        #Adding a new record with a child and parent HP terms
        new_record = [2000, ["HP:0001513", "HP:0025500"], [["21", 20, 25], ["X", 1000, 2000]]]
        extra_attr = {"sex": "M"}
        self.patient_data.add_record(new_record, extra_attr)

        self.patient_data.link2ont(Cohort.act_ont)
        profile_sizes, parental_terms_per_profile = self.patient_data.get_profile_redundancy()
        expected_profile_sizes = (3, 3, 2, 2, 1) #From patients sorted from higher to lower number of HP terms
        expected_parental_terms_per_profile = (0, 0, 1, 0, 0) #Same order as above

        self.assertEqual(profile_sizes, expected_profile_sizes)
        self.assertEqual(parental_terms_per_profile, expected_parental_terms_per_profile)


    def test_get_profiles_terms_frequency(self):
        self.patient_data.link2ont(Cohort.act_ont)
        term_stats = self.patient_data.get_profiles_terms_frequency()
        expected_term_stats = [['Intellectual disability', 1.0], ['Autism', 0.25], ['Hypotonia', 0.25], 
                               ['Turricephaly', 0.25], ['Hearing impairment', 0.25], ['Microcephaly', 0.25]]
        self.assertEqual(term_stats, expected_term_stats)


    def test_compute_term_list_and_childs(self):
        self.patient_data.link2ont(Cohort.act_ont)
        suggested_childs, term_with_childs_ratio = self.patient_data.compute_term_list_and_childs()
        expected_suggested_childs = {'132': [[['HP:0000717', 'Autism'], []], [['HP:0001252', 'Hypotonia'], [['HP:0000297', 'Facial hypotonia'], ['HP:0001290', 'Generalized hypotonia'], ['HP:0001319', 'Neonatal hypotonia'], ['HP:0003397', 'Generalized hypotonia due to defect at the neuromuscular junction'], ['HP:0006829', 'Severe muscular hypotonia'], ['HP:0006852', 'Episodic generalized hypotonia'], ['HP:0008935', 'Generalized neonatal hypotonia'], ['HP:0008936', 'Axial hypotonia'], ['HP:0008947', 'Infantile muscular hypotonia'], ['HP:0009062', 'Infantile axial hypotonia'], ['HP:0012389', 'Appendicular hypotonia'], ['HP:0030190', 'Oral motor hypotonia'], ['HP:0031139', 'Frog-leg posture']]], [['HP:0001249', 'Intellectual disability'], [['HP:0001256', 'Intellectual disability, mild'], ['HP:0002187', 'Intellectual disability, profound'], ['HP:0002342', 'Intellectual disability, moderate'], ['HP:0006887', 'Intellectual disability, progressive'], ['HP:0006889', 'Intellectual disability, borderline'], ['HP:0010864', 'Intellectual disability, severe']]]], 
                            '599': [[['HP:0001249', 'Intellectual disability'], [['HP:0001256', 'Intellectual disability, mild'], ['HP:0002187', 'Intellectual disability, profound'], ['HP:0002342', 'Intellectual disability, moderate'], ['HP:0006887', 'Intellectual disability, progressive'], ['HP:0006889', 'Intellectual disability, borderline'], ['HP:0010864', 'Intellectual disability, severe']]]], 
                            '647': [[['HP:0001249', 'Intellectual disability'], [['HP:0001256', 'Intellectual disability, mild'], ['HP:0002187', 'Intellectual disability, profound'], ['HP:0002342', 'Intellectual disability, moderate'], ['HP:0006887', 'Intellectual disability, progressive'], ['HP:0006889', 'Intellectual disability, borderline'], ['HP:0010864', 'Intellectual disability, severe']]], [['HP:0000262', 'Turricephaly'], [['HP:0000244', 'Brachyturricephaly'], ['HP:0000263', 'Oxycephaly']]]], 
                            '648': [[['HP:0000365', 'Hearing impairment'], [['HP:0000399', 'Prelingual sensorineural hearing impairment'], ['HP:0011474', 'Childhood onset sensorineural hearing impairment'], ['HP:0000407', 'Sensorineural hearing impairment'], ['HP:0000405', 'Conductive hearing impairment'], ['HP:0000408', 'Progressive sensorineural hearing impairment'], ['HP:0001730', 'Progressive hearing impairment'], ['HP:0000410', 'Mixed hearing impairment'], ['HP:0001757', 'High-frequency sensorineural hearing impairment'], ['HP:0005101', 'High-frequency hearing impairment'], ['HP:0008504', 'Moderate sensorineural hearing impairment'], ['HP:0012713', 'Moderate hearing impairment'], ['HP:0008513', 'Bilateral conductive hearing impairment'], ['HP:0008527', 'Congenital sensorineural hearing impairment'], ['HP:0008542', 'Low-frequency hearing loss'], ['HP:0008573', 'Low-frequency sensorineural hearing impairment'], ['HP:0008587', 'Mild neurosensory hearing impairment'], ['HP:0012712', 'Mild hearing impairment'], ['HP:0008591', 'Congenital conductive hearing impairment'], ['HP:0008596', 'Postlingual sensorineural hearing impairment'], ['HP:0008598', 'Mild conductive hearing impairment'], ['HP:0008607', 'Progressive conductive hearing impairment'], ['HP:0008610', 'Infantile sensorineural hearing impairment'], ['HP:0008615', 'Adult onset sensorineural hearing impairment'], ['HP:0008619', 'Bilateral sensorineural hearing impairment'], ['HP:0008625', 'Severe sensorineural hearing impairment'], ['HP:0012714', 'Severe hearing impairment'], ['HP:0009900', 'Unilateral deafness'], ['HP:0011476', 'Profound sensorineural hearing impairment'], ['HP:0012715', 'Profound hearing impairment'], ['HP:0011975', 'Aminoglycoside-induced hearing loss'], ['HP:0012716', 'Moderate conductive hearing impairment'], ['HP:0012717', 'Severe conductive hearing impairment'], ['HP:0012779', 'Transient hearing impairment'], ['HP:0012781', 'Mid-frequency hearing loss'], ['HP:0040113', 'Old-aged sensorineural hearing impairment'], ['HP:0040119', 'Unilateral conductive hearing impairment']]], [['HP:0000252', 'Microcephaly'], [['HP:0000253', 'Progressive microcephaly'], ['HP:0005484', 'Secondary microcephaly'], ['HP:0004485', 'Cessation of head growth'], ['HP:0011451', 'Primary microcephaly']]], [['HP:0001249', 'Intellectual disability'], [['HP:0001256', 'Intellectual disability, mild'], ['HP:0002187', 'Intellectual disability, profound'], ['HP:0002342', 'Intellectual disability, moderate'], ['HP:0006887', 'Intellectual disability, progressive'], ['HP:0006889', 'Intellectual disability, borderline'], ['HP:0010864', 'Intellectual disability, severe']]]]}
        expected_term_with_childs_ratio = 0.8888888888888888
        self.assertEqual(suggested_childs, expected_suggested_childs)
        self.assertEqual(term_with_childs_ratio, expected_term_with_childs_ratio)


    def test_get_profile_ontology_distribution_tables(self):
        self.patient_data.link2ont(Cohort.act_ont)
        ontology_levels, distribution_percentage = self.patient_data.get_profile_ontology_distribution_tables()
        expected_ontology_levels = [['level', 'ontology', 'cohort'], [1, 1, 0], [2, 6, 0], [3, 62, 0], [4, 304, 0], [5, 882, 0], [6, 2234, 5], [7, 3810, 2], [8, 3802, 2], [9, 3008, 0], [10, 1914, 0], [11, 709, 0], [12, 348, 0], [13, 106, 0], [14, 31, 0], [15, 13, 0], [16, 2, 0]]
        expected_distribution_percentage = [['level', 'ontology', 'weighted cohort', 'uniq terms cohort'], [1, 0.006, 0.0, 0.0], [2, 0.035, 0.0, 0.0], [3, 0.36, 0.0, 0.0], [4, 1.764, 0.0, 0.0], [5, 5.118, 0.0, 0.0], [6, 12.964, 55.556, 33.333], [7, 22.11, 22.222, 33.333], [8, 22.064, 22.222, 33.333], [9, 17.456, 0.0, 0.0], [10, 11.107, 0.0, 0.0], [11, 4.114, 0.0, 0.0], [12, 2.019, 0.0, 0.0], [13, 0.615, 0.0, 0.0], [14, 0.18, 0.0, 0.0], [15, 0.075, 0.0, 0.0], [16, 0.012, 0.0, 0.0]]
        self.assertEqual(ontology_levels, expected_ontology_levels)
        self.assertEqual(distribution_percentage, expected_distribution_percentage)

    def test_get_ic_analysis(self):
        self.patient_data.link2ont(Cohort.act_ont)
        onto_ic, freq_ic, onto_ic_profile, freq_ic_profile = self.patient_data.get_ic_analysis()
        expected_onto_ic = {'HP:0001252': 3.0902076502756683, 'HP:0000365': 2.6681339618869115, 'HP:0001249': 3.3912376459396496, 'HP:0000262': 3.759214431234244, 'HP:0000252': 3.537365681617888, 'HP:0000717': 4.2363356859539065}
        expected_freq_ic = {'HP:0001252': 0.9030899869919435, 'HP:0000365': 0.9030899869919435, 'HP:0001249': 0.3010299956639812, 'HP:0000262': 0.9030899869919435, 'HP:0000252': 0.9030899869919435, 'HP:0000717': 0.9030899869919435}
        expected_onto_ic_profile = {'132': 3.5725936607230744, '599': 3.3912376459396496, '647': 3.5752260385869468, '648': 3.198912429814816}
        expected_freq_ic_profile = {'132': 0.7024033232159561, '599': 0.3010299956639812, '647': 0.6020599913279624, '648': 0.7024033232159561}

        self.assertEqual(onto_ic, expected_onto_ic)
        self.assertEqual(freq_ic, expected_freq_ic)
        self.assertEqual(onto_ic_profile, expected_onto_ic_profile)
        self.assertEqual(freq_ic_profile, expected_freq_ic_profile)


    def test_get_profiles_mean_size(self):
        self.patient_data.link2ont(Cohort.act_ont)
        profile_mean_size = self.patient_data.get_profiles_mean_size()
        self.assertEqual(profile_mean_size, 2.25)


    def test_get_profile_length_at_percentile(self):
        self.patient_data.link2ont(Cohort.act_ont)
        length_percent = self.patient_data.get_profile_length_at_percentile()
        self.assertEqual(length_percent, 2.5)


    def test_get_dataset_specifity_index(self):
        self.patient_data.link2ont(Cohort.act_ont)
        
        dsi_uniq = self.patient_data.get_dataset_specifity_index("uniq")
        dsi_weigthed = self.patient_data.get_dataset_specifity_index("weigthed")

        expected_dsi_uniq = 0.16867992874998128
        expected_dsi_weigthed = 0.0014407344880051807
        
        self.assertEqual(dsi_uniq, expected_dsi_uniq)
        self.assertEqual(dsi_weigthed, expected_dsi_weigthed)


    def test_compare_profiles(self):
        self.patient_data.link2ont(Cohort.act_ont)
        similarities = self.patient_data.compare_profiles()
        exptected_similarities = {'132': {'132': 3.5725936607230744, '599': 2.104921447549329, '647': 1.9166626779429385, '648': 1.6416691014914389}, 
                                  '599': {'132': 2.104921447549329, '599': 3.3912376459396496, '647': 2.2639587126068683, '648': 1.9089226834606448}, 
                                  '647': {'132': 1.9166626779429385, '599': 2.2639587126068683, '647': 3.5752260385869468, '648': 2.0721949277359677}, 
                                  '648': {'132': 1.6416691014914389, '599': 1.9089226834606448, '647': 2.0721949277359677, '648': 3.198912429814816}}
        self.assertEqual(similarities, exptected_similarities)


    def test_index_vars(self):
        #expected_number_of_regions = sum([len(regions) for pat_gen_feats in self.patient_data.vars.values() 
        #                                for chrm, regions in pat_gen_feats.regions.items()])
        expected_number_of_regions = 5
        
        #Mixing the genomic regions of the patients in a unique Genomic Feature
        self.patient_data.index_vars()
        returned_number_of_regions = sum([len(regions) for chrm, regions in self.patient_data.var_idx.regions.items()])
        #Checking that the number of regions in this mixed genomic feature is the same as the total number of regions of the patients
        self.assertEqual(expected_number_of_regions, returned_number_of_regions)


    def test_get_vars_sizes(self):
        expected_sizes = [51, 200, 501, 51, 10001] #For every region of the patients toy dataset: region["stop"] - region["start"] + 1

        #Mixing the genomic regions of the patients in a unique Genomic Feature (var_idx attribute)
        self.patient_data.index_vars()
        #Getting the sizes of each of the genomic regions of the patients
        returned_sizes = self.patient_data.get_vars_sizes()

        self.assertEqual(sorted(expected_sizes), sorted(returned_sizes))

        #Getting summary sizes, that is a 2D list with the size and number of patients with that size
        returned_sizes = self.patient_data.get_vars_sizes(summary=True)
        expected_sizes = [[51, 2], [10001, 1], [501, 1], [200, 1]]    
        self.assertEqual(expected_sizes, returned_sizes)


    def test_generate_cluster_regions(self):
        #Mixing the genomic regions of the patients in a unique Genomic Feature
        self.patient_data.index_vars()

        returned_ids_by_cluster_1, returned_annotated_full_red_1 = self.patient_data.var_idx.generate_cluster_regions(meth="reg_overlap", tag="coh", ids_per_reg = 0)

        #We expect all the patients to be returned when the filtering threshold is 0 and at least one region per patient, and each region to have 4 elements: chrm, start, stop and the region_id
        self.assertDictEqual(returned_ids_by_cluster_1, {'132': ['22.0.coh.1', '22.1.coh.2'], '648': ['22.1.coh.2', '22.2.coh.1', '22.4.coh.1'], '647': ['22.3.coh.1'], '599': ['15.0.coh.1']})
        self.assertEqual(sorted(returned_annotated_full_red_1), sorted([[50, 52, '22', '22.0.coh.1'], [52, 100, '22', '22.1.coh.2'], [100, 102, '22', '22.2.coh.1'], [500, 1000, '22', '22.3.coh.1'], [10000, 20000, '22', '22.4.coh.1'], [1, 200, '15', '15.0.coh.1']]))

        ##We expect only one (overlapping) regions with more than 1 patients (region 52-100 in chromosome 22, shared by patients 132 and 648)
        returned_ids_by_cluster_2, returned_annotated_full_red_2 = self.patient_data.var_idx.generate_cluster_regions(meth="reg_overlap", tag="coh", ids_per_reg = 1) #ids_per_reg is exclusive, so we expect more than 1 patient per region
        self.assertDictEqual(returned_ids_by_cluster_2, {'132': ['22.0.coh.2'], '648': ['22.0.coh.2']})
        self.assertEqual(sorted(returned_annotated_full_red_2), sorted([[52, 100, '22', '22.0.coh.2']]))

        for reference_region in returned_annotated_full_red_2:
            #We expect the last number of the region tag (number of patients overlapping the region) to be greater than the threshold we used
            self.assertGreaterEqual(int(reference_region[3].split(".")[-1]), 1)
            #And also that the first number in the region tag corresponds to the chromosome
            self.assertEqual(reference_region[3].split(".")[0], reference_region[2]) 


        returned_ids_by_cluster_1000, returned_annotated_full_red_1000 = self.patient_data.var_idx.generate_cluster_regions(meth="reg_overlap", tag="coh", ids_per_reg = 1000)
        self.assertEqual(returned_ids_by_cluster_1000, {}) #We expect no patients to be returned when the filtering threshold is 1000 (not that number of overlaps in the toy dataset)
        self.assertEqual(returned_annotated_full_red_1000, []) #Neither we expect any reference regions to be returned


    def test_save(self):
        ### Creating the tmp folder and saving the toy dataset in paco format ###
        os.makedirs(os.path.join(ROOT_PATH, "tmp"), exist_ok=True)
        tmp_file = os.path.join(ROOT_PATH, "tmp", "test_paco_format.txt")
        self.patient_data.save(f"{tmp_file}", mode="paco", translate=True)

        ### Checking that the file has been created and has the expected number of lines (equals to input file)###
        check_n_lines = f'wc -l {tmp_file} | cut -d " " -f 1'
        ps2 = subprocess.Popen(check_n_lines,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output2 = ps2.communicate()[0]
        total_lines = str(output2).replace("b'","").replace("\\n'", "").strip()
        self.assertEqual(total_lines, "6")

        ### Checking that the HP terms has been translated to descriptions, so there are no HP terms in the file ###
        check_no_HP = f'grep HP: {tmp_file} | wc -l'
        ps = subprocess.Popen(check_no_HP,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        lines_HP = str(output).replace("b'","").replace("\\n'", "").strip()
        self.assertEqual(lines_HP, "0")

        ### Loading back the saved file and checking that the dataset is the same as the original one ###
        # (and also that descriptions have been translated back to HP terms)
        options = {"input_file": tmp_file,
                   "id_col":"id", "chromosome_col": "chr", 
                   "start_col":"start", "end_col":"stop", "ont_col":"terms",
                   "separator": "|", "header": True, "names": True}

        patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
        
        #Checking that patient_id and phenotypic profile are the same
        self.assertEqual([[pair[0], sorted(pair[1])] for pair in patient_data.profiles.items()], 
                         [[pair[0], sorted(pair[1])] for pair in self.patient_data.profiles.items()])
        
        #Checking that patient_id and genomic regions are the same
        #(doing for loop because the region identifiers ("to") could be different because of
        # the order of saving the patients in the file)
        self.assertEqual(patient_data.vars.keys(), self.patient_data.vars.keys())
        for pat_id, genomic_region in patient_data.vars.items():
            for chrm, regions in genomic_region.regions.items():
                original_region = self.patient_data.vars[pat_id].regions[chrm]
                saved_region = patient_data.vars[pat_id].regions[chrm]
                self.assertEqual([[region["start"],region["stop"],region["chrm"]] for region in original_region],
                                 [[region["start"],region["stop"],region["chrm"]] for region in saved_region])

        os.remove(f"{tmp_file}")

    def test_export_phenopackets(self):
        os.makedirs(os.path.join(ROOT_PATH, "tmp", "phenopackets"), exist_ok=True)
        tmp_folder = os.path.join(ROOT_PATH, "tmp", "phenopackets")
        self.patient_data.export_phenopackets(f"{tmp_folder}", "hg38")

        ### Checking that the file has been created
        ### and loading the file back and checking that is can be correctly parsed as a json ###
        for id in self.patient_data.profiles.keys():
            pheno_id_file = os.path.join(tmp_folder, id + ".json")
            self.assertTrue(os.path.exists(pheno_id_file))
            with open(pheno_id_file, "r") as f:
                data = json.load(f)
                self.assertTrue(isinstance(data, dict))
        
        ###Checking that there are the same number of files as patients in the dataset ###
        check_n_files = f'ls {tmp_folder} | wc -l'
        ps2 = subprocess.Popen(check_n_files,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output2 = ps2.communicate()[0]
        total_files = int(str(output2).replace("b'","").replace("\\n'", "").strip())

        self.assertEqual(total_files, 4)
        
        for id in self.patient_data.profiles.keys():
            pheno_id_file = os.path.join(tmp_folder, id + ".json")
            os.remove(f"{pheno_id_file}")

    def test_process_dummy_clustered_patients(self):
        tmp_folder = os.path.join(ROOT_PATH, "tmp", "dummy_cluster")
        os.makedirs(tmp_folder, exist_ok=True)

        options = {"chromosome_col": True, "clusters2show_detailed_phen_data": 3}
        all_ics, all_lengths, cluster_data_by_chromosomes, top_cluster_phenotypes, multi_chromosome_patients = self.patient_data.process_dummy_clustered_patients(options, phenotype_ic=self.ic_hpos, temp_folder=tmp_folder)
        
        self.assertEqual(all_ics, [[0.4615785999681353, 1.856060288989204]])
        self.assertEqual(all_lengths, [[1, 2]])
        self.assertEqual(cluster_data_by_chromosomes, [[1, 2, '15', 1], [1, 2, '22', 1]])
        self.assertEqual(top_cluster_phenotypes, [[['Intellectual disability'], ['Intellectual disability', 'Turricephaly']]])
        self.assertEqual(multi_chromosome_patients, 2)
        
        for file in os.listdir(tmp_folder):
            os.remove(os.path.join(tmp_folder, file))

    def test_dummy_cluster_patients(self):        
        tmp_folder = os.path.join(ROOT_PATH, "tmp", "dummy_cluster")
        os.makedirs(tmp_folder, exist_ok=True)
        
        clustered_patients = self.patient_data.dummy_cluster_patients(temp_folder=tmp_folder)
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, "cluster_asignation")))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, "pat_hpo_matrix.npy")))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, "pat_hpo_matrix_x.lst")))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, "pat_hpo_matrix_y.lst")))
        self.assertEqual(clustered_patients, {0: ['132'], 1: ['599', '647'], 2: ['648']})

        for file in os.listdir(tmp_folder):
            os.remove(os.path.join(tmp_folder, file))

    def test_process_cluster(self):
        patients_id = ["599", "647"]
        options = {"chromosome_col": True, "clusters2show_detailed_phen_data": 3}
        ont = self.patient_data.get_ontology(Cohort.act_ont)
        
        chrs, all_phens, profile_ics, profile_lengths = self.patient_data.process_cluster(patients_id, self.ic_hpos, options, ont, 0)
        self.assertEqual(chrs, {'15': 1, '22': 1})
        self.assertEqual(all_phens, [['Intellectual disability'], ['Intellectual disability', 'Turricephaly']])
        self.assertEqual(profile_ics, [0.4615785999681353, 1.856060288989204])
        self.assertEqual(profile_lengths, [1, 2])

    def test_get_profile_ic(self):
        hpos = ['HP:0000717', 'HP:0001252']
        profile_ic = (self.ic_hpos[hpos[0]] + self.ic_hpos[hpos[1]]) / 2
        returned = self.patient_data.get_profile_ic(hpos, self.ic_hpos)
        self.assertEqual(returned, profile_ic)

    def test_get_matrix_similarity(self): #Helper function tested in get_similarity_clusters (its higher order function)
        pass
        #self.assertFalse(True, "Helper function tested in get_similarity_clusters")

    def test_get_similarity_clusters(self):
        options = {"sim_thr": 0.3}
        method_name = "resnik"
        tmp_folder = os.path.join(ROOT_PATH, "tmp", "similarity_cluster")
        os.makedirs(tmp_folder, exist_ok=True)

        expected_clusters = {6: ['648', '132', '599', '647']}
        expected_sim_matrix =  np.array([[0.         , 1.6416691  , 1.90892268 , 2.07219493],
                                        [1.6416691  , 0.         , 2.10492145 , 1.91666268],
                                        [1.90892268 , 2.10492145 , 0.         , 2.26395871],
                                        [2.07219493 , 1.91666268 , 2.26395871 , 0.        ]])
        
        expected_linkage = np.array([[2.         , 3.         , 0.         , 2.        ],
                                    [1.         , 4.         , 0.31188394 , 3.        ],
                                    [0.         , 5.         , 0.50071574 , 4.        ]])
        
        expected_raw_cls = np.array([[6], [6], [6], [6]])

        #Testing for the first time (no files in tmp_folder, so values have to be calculated)
        clusters, similarity_matrix, linkage, raw_cls = self.patient_data.get_similarity_clusters(method_name, "hpo", options, temp_folder=tmp_folder, reference_profiles=None)

        self.assertEqual(clusters, expected_clusters)
        self.assertTrue(np.isclose(similarity_matrix, expected_sim_matrix).all())
        self.assertTrue(np.isclose(linkage, expected_linkage).all())
        self.assertTrue((raw_cls == expected_raw_cls).all())

        self.assertTrue(os.path.exists(os.path.join(tmp_folder, f"similarity_matrix_{method_name}.npy")))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, f"similarity_matrix_{method_name}_x.lst")))
        self.assertFalse(os.path.exists(os.path.join(tmp_folder, f"similarity_matrix_{method_name}_y.lst")))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, f"{method_name}_clusters.txt")))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, f'profiles_similarity_{method_name}.txt')))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, f'{method_name}_linkage.npy')))
        self.assertTrue(os.path.exists(os.path.join(tmp_folder, f'{method_name}_raw_cls.npy')))

        #Testing for the second time (files in tmp_folder, so values will be loaded from files instead of being calculated)
        clusters, similarity_matrix, linkage, raw_cls = self.patient_data.get_similarity_clusters(method_name, "hpo", options, temp_folder=tmp_folder, reference_profiles=None)
        self.assertEqual(clusters, expected_clusters)
        self.assertTrue(np.isclose(similarity_matrix, expected_sim_matrix).all())
        self.assertTrue(np.isclose(linkage, expected_linkage).all())
        self.assertTrue((raw_cls == expected_raw_cls).all())

        for file in os.listdir(tmp_folder):
            os.remove(os.path.join(tmp_folder, file))



    def test_calc_sim_term2term_similarity_matrix(self):
        reference_profile = ["HP:0000365", "HP:0001249"]
        ref_profile_id = "657"
        ref_profile_dict = {ref_profile_id: reference_profile}
        ont = self.patient_data.get_ontology(Cohort.act_ont)
        ont.load_profiles(ref_profile_dict)
        
        candidate_sim_matrix, candidates, candidates_ids = self.patient_data.calc_sim_term2term_similarity_matrix(reference_profile, ref_profile_id, self.patient_data.profiles, ont)
        self.assertEqual(candidate_sim_matrix, [['Intellectual disability', 1.0, 1.0, 1.0, 1.0], ['Hearing impairment', 1.0, 0, 0.002925264157569972, 0.0032651226948571363]])
        self.assertEqual(candidates, [['648', 0.8487148451790173], ['599', 0.6677009703920632], ['647', 0.5015070438334399], ['132', 0.48667053298635815]])
        self.assertEqual(candidates_ids, ['648', '599', '647', '132'])


    def test_get_term2term_similarity_matrix(self): #Helper function of calc_sim_term2term_similarity_matrix
        pass
        #self.assertFalse(True, "Helper function of calc_sim_term2term_similarity_matrix")

    def test_get_detailed_similarity(self): #Helper function of get_term2term_similarity_matrix
        pass
        #self.assertFalse(True, "Helper function of get_term2term_similarity_matrix")

    def test_write_detailed_hpo_profile_evaluation(self): #TODO: test this function
        tmp_folder = os.path.join(ROOT_PATH, "tmp", "profile_evaluation")
        os.makedirs(tmp_folder, exist_ok=True)
        suggested_childs = {
         "132":[[["HP:0000717","Autism"],[["HP:0001290", "Generalized hypotonia"]]], 
                [["HP:0001252","Hypotonia"],[]], 
                [["HP:0001249","Intellectual disability"],[["HP:0001263", "Global developmental delay"], ["HP:0005280", "Depressed nasal bridge"]]], 
                [["HP:0000262","Turricephaly"],[["HP:0011800", "Midface retrusion"], ["HP:0005280", "Depressed nasal bridge"], ["HP:0000358", " Posteriorly rotated ears "]]]],
        "599": [ [["HP:0000365","Hearing impairment"],[["HP:0000175", "Cleft palate"]] ] ]}
        self.patient_data.write_detailed_hpo_profile_evaluation(suggested_childs, os.path.join(tmp_folder, "file.csv"))

        self.assertTrue(os.path.exists(os.path.join(tmp_folder, "file.csv")))
        lines = subprocess.run("wc -l ./tests/tmp/profile_evaluation/file.csv", shell=True, capture_output=True, text=True).stdout.split(" ")[0]
        patients = subprocess.run("grep PATIENT ./tests/tmp/profile_evaluation/file.csv | wc -l", shell=True, capture_output=True, text=True).stdout.split(" ")[0]
        few_HPs_warning = subprocess.run("grep WARNING ./tests/tmp/profile_evaluation/file.csv | wc -l", shell=True, capture_output=True, text=True).stdout.split(" ")[0]
                
        self.assertEqual(int(lines), 14)
        self.assertEqual(int(patients), 2)
        self.assertEqual(int(few_HPs_warning), 1)
        
        os.remove(os.path.join(tmp_folder, "file.csv"))
        

    def test_write_profile_pairs(self):
        pairs = {"A": {"B": 3, "C": 4}, "B": {"C": 5, "D": 9}, "C": {"D": 6}}
        tmp_folder = os.path.join(ROOT_PATH, "tmp", "dummy_profile")
        filename = os.path.join(tmp_folder, "profile_pairs.txt")
        os.makedirs(tmp_folder, exist_ok=True)

        self.patient_data.write_profile_pairs(pairs, filename)
        self.assertTrue(os.path.exists(filename))

        filee = open(filename)
        self.assertEqual("A\tB\t3\nA\tC\t4\nB\tC\t5\nB\tD\t9\nC\tD\t6\n", filee.read())
        filee.close()

        for file in os.listdir(tmp_folder):
            os.remove(os.path.join(tmp_folder, file))