import sys, os, json, unittest 
import subprocess
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
PATIENTS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, '100_test_dataset_with_sex.txt'))
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))

with open(CONSTANTS_PATH) as infile:
    exec(infile.read())



class CohortTestSuite(unittest.TestCase):
    def setUp(self):
        patients = [
                    ["132", ["HP:0000717", "HP:0001252", "HP:0001249"],	[["22", 50, 100]], {"sex": "M"}],
                    ["599", ["HP:0001249"],	[["15", 1, 200]], {"sex": "F"}],
                    ["647", ["HP:0001249", "HP:0000262"], [["22", 500, 1000]], {"sex": "M"}],
                    ["648", ["HP:0000365", "HP:0000252", "HP:0001249"], [["22", 52, 102], ["22", 10000, 20000]], {"sex": "F"}]
                    ]
        self.n_patients = len(patients)


        hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
        Cohort.load_ontology("hpo", hpo_file)
        self.patient_data = Cohort()

        for patient in patients:
            id, terms, variants, other_attr = patient
            self.patient_data.add_record([id, terms, variants], other_attr)

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


#    def test_get_profile_redundancy(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        profile_sizes, parental_terms_per_profile = self.patient_data.get_profile_redundancy()
#        print(profile_sizes)
#        print(parental_terms_per_profile)
#
#
##    def test_get_profiles_terms_frequency(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        term_stats = self.patient_data.get_profiles_terms_frequency()
##        print(term_stats)
##
##
##    def test_compute_term_list_and_childs(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        #suggested_childs, term_with_childs_ratio = self.patient_data.compute_term_list_and_childs()
##        #print(suggested_childs)
##        #print(term_with_childs_ratio)
##
##
##    def test_get_profile_ontology_distribution_tables(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        #ontology_levels, distribution_percentage = self.patient_data.get_profile_ontology_distribution_tables()
##        #print(ontology_levels)
##        #print(distribution_percentage)
##
##
##    def test_get_ic_analysis(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        onto_ic, freq_ic, onto_ic_profile, freq_ic_profile = self.patient_data.get_ic_analysis()
##        print(onto_ic)
##        print(freq_ic)
##        print(onto_ic_profile)
##        print(freq_ic_profile)
##
##
##    def test_get_profiles_mean_size(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        profile_mean_size = self.patient_data.get_profiles_mean_size()
##        print(profile_mean_size)
##
##
##    def test_get_profile_length_at_percentile(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        length_percent = self.patient_data.get_profile_length_at_percentile()
##        print(length_percent)
##
##
##    def test_get_dataset_specifity_index(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        #dsi_uniq = self.patient_data.get_dataset_specifity_index("uniq")
##        #dsi_weigthed = self.patient_data.get_dataset_specifity_index("weigthed")
##        #print(dsi_uniq, dsi_weigthed)
##
##    def test_compare_profiles(self):
##        self.patient_data.link2ont(Cohort.act_ont)
##        similarities = self.patient_data.compare_profiles()
#
#
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
        expected_sizes = [(51, 2), (10001, 1), (501, 1), (200, 1)]    
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
        self.assertEqual(patient_data.profiles.items(), self.patient_data.profiles.items())
        
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
            pheno_id_file = os.path.join(tmp_file, id + ".json")
            os.remove(f"{pheno_id_file}")
