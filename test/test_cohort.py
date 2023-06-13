import sys 
import unittest
import os
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
        hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
        #Cohort.load_ontology("hpo", hpo_file, "./test/data/excluded.txt") #TODO: Ask PSZ about this
        Cohort.load_ontology("hpo", hpo_file )
        Cohort.act_ont = "hpo"
        options = {"input_file": PATIENTS_FILE,
                   "id_col":"patient_id", "chromosome_col": "chr", 
                   "header": True, "separator": "|", "names": True,
                   "start_col":"start", "end_col":"end", "ont_col":"phenotypes",
                   "sex_col": "sex"}
        self.patient_data, self.rejected_hpos_L, self.rejected_patients_L = Cohort_Parser.load(options)

        self.n_patients = 84


    def test_load(self):
        chroms = [str(chrm) for chrm in range(1,23)] + ['X', 'Y']
        
        self.assertEqual(sorted(self.rejected_hpos_L),
                         sorted(['Generalized tonic-clonic seizures', 'Stereotypy', 'Sparse and thin eyebrow', 'obsolete Prenatal short stature', 'Abnormality of the pinna', 'Capillary hemangiomas']))
        self.assertEqual(self.rejected_patients_L, [])
        
        self.assertEqual(len(self.patient_data.profiles.keys()), self.n_patients) #Checking that all patients are loaded with their phenotypes
        self.assertEqual(len(self.patient_data.vars.keys()), self.n_patients) #Checking that all patients are loaded with their genomic regions


        for patient, regions in self.patient_data.vars.items():
            self.assertIsInstance(regions, Genomic_Feature) #Checking that all patients are loaded with their genomic regions (as Genomic_Feature objects)
            for chromosome in regions.regions.keys():
                self.assertIn(chromosome, chroms) #Checking that the chromosomes are correct
                for region in regions.regions[chromosome]:
                    self.assertEqual(set(region.keys()), 
                                     {'start', 'stop', 'to', 'chrm'}) #Checking that each region has the correct attributes
                    self.assertIs(type(region['start']), int) #Checking that the start position is an integer
                    self.assertIs(type(region['stop']), int) #Checking that the stop position is an integer


        for patient, hpterms in self.patient_data.profiles.items(): 
            self.assertGreater(len(hpterms), 0) #Checking that all patients are loaded with their phenotypes. Otherwise, they should have been discarded in "rejected_patients_L"
            self.assertEqual(len(set(hpterms)), len(hpterms)) #Checking that there are no repeated phenotypes for each patient
            for hpterm in hpterms:
                self.assertIn("HP:", hpterm) #Checking HPs has been translated to HP:XXXXXXX if name/description was given


    
        for patient, attrs in self.patient_data.extra_attr.items():
            for attr, value in attrs.items():
                self.assertIn(attr, ["sex"]) #Checking that additional attributes are correctly saved (only sex for now)
                self.assertIn(value, ["M", "F"])


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
        #HP:0000486: Strabismus, just 8 patients in the toy dataset
        def check_hp_in_pat(query):
            def _check_hp_in_pat(pat, hps):
                return query in hps
            return _check_hp_in_pat
            
        self.patient_data.select_by_profile(check_hp_in_pat("HP:0000486"))

        self.assertEqual(len(self.patient_data.profiles.keys()), 8)
        self.assertEqual(len(self.patient_data.vars.keys()), 8)
        #TODO: Ask PSZ wether extra_attr should be removed as well
        #self.assertEqual(len(self.patient_data.extra_attr.keys()), 8)

        for pat, hpterms in self.patient_data.profiles.items():
            self.assertIn("HP:0000486", hpterms)


    def test_select_by_var(self):
        #There are only 5 patients with mutations in the X chromosome in the toy dataset
        def check_var_in_pat(query):
            def _check_var_in_pat(pat, vars):
                return query in vars.regions.keys()
            return _check_var_in_pat
        
        self.patient_data.select_by_var(check_var_in_pat("X"))

        self.assertEqual(len(self.patient_data.profiles.keys()), 5)
        self.assertEqual(len(self.patient_data.vars.keys()), 5)

        for pat, vars in self.patient_data.vars.items():
            self.assertIn("X", vars.regions.keys())


    def test_filter_by_term_number(self):
        #We are filtering out patients with less than 4 phenotypes
        self.patient_data.filter_by_term_number(4)

        for pat, hpterms in self.patient_data.profiles.items():
            self.assertGreaterEqual(len(hpterms), 4)

        self.assertGreaterEqual(84, len(self.patient_data.profiles.keys())) #There are 84 patients in the toy dataset before the filtering, so now it should be less


    def test_remove_incomplete_records(self):
        #We are checking that there are 84 patients in the toy dataset before the cleaning
        self.assertEqual(len(self.patient_data.profiles.keys()), self.n_patients)
        self.assertEqual(len(self.patient_data.vars.keys()), self.n_patients)

        self.patient_data.profiles["130"] = [] #making a patient to have no phenotype
        self.patient_data.vars["132"] = Genomic_Feature([]) #making another patient to have no genomic regions
        
        #We are filtering out patients without phenotype or genomic regions
        self.patient_data.remove_incomplete_records()

        #Now there should be 82 patients
        self.assertEqual(len(self.patient_data.profiles.keys()), self.n_patients - 2)
        self.assertEqual(len(self.patient_data.vars.keys()), self.n_patients - 2)


    def test_add_gen_feat(self):
        new_gen_feat = [["21", 20, 25, 0]]
        self.patient_data.add_gen_feat(2000, new_gen_feat)
        self.assertEqual(self.patient_data.vars.get(2000).regions.get("21"), [{'chrm': '21', 'start': 20, 'stop': 25, 'to': 0}])


    def test_get_profile(self):
        self.assertEqual(self.patient_data.get_profile("130"), self.patient_data.profiles["130"])


    def test_get_vars(self):
        self.assertEqual(self.patient_data.get_vars("130"), self.patient_data.vars["130"])


    def test_each_profile(self):
        for pat, hpterms in self.patient_data.each_profile():
            self.assertIn(pat, self.patient_data.profiles.keys())
            self.assertIn(hpterms, self.patient_data.profiles.values())


    def test_each_var(self):
        for pat, vars in self.patient_data.each_var():
            self.assertIn(pat, self.patient_data.vars.keys())
            self.assertIn(vars, self.patient_data.vars.values())


    def test_get_general_profile(self):
        #Intellectual disability HP:0001249 is the term with a frequency higher than 0.5 in the toy dataset
        general_profiles = self.patient_data.get_general_profile(thr=0.5)
        self.assertEqual(general_profiles, ["HP:0001249"])


    def test_check(self):
        #Adding a new record with an incorrect HP term
        new_record_bad = [2000,#id
                        ["HP:BADCODE"], #phenotype
                        [["21", 20, 25], ["X", 1000, 2000]]]
        
        #Adding a new record with an obsolete HP term
        new_record_obs = [2001,#id 
                      ["HP:0000284"], #phenotype
                      [["21", 20, 25], ["X", 1000, 2000]]] #genomic regions
        self.patient_data.profiles[new_record_obs[0]] = new_record_obs[1] #adding the new record to the profiles
        self.patient_data.profiles[new_record_bad[0]] = new_record_bad[1]
        self.patient_data.add_gen_feat(2000, new_record_bad[2]) #adding the new record to the genomic regions
        self.patient_data.add_gen_feat(2001, new_record_obs[2])
        
        self.assertEqual(len(self.patient_data.profiles.items()), self.n_patients + 2)
        self.assertEqual(len(self.patient_data.vars.items()), self.n_patients + 2)

        #Checking the dataset. The incorrect term should be rejected (and also the patient, because s/he only has one HP term), 
        # but not the obsolete one, that will be replaced by its current equivalent
        rejected_terms, rejected_recs = self.patient_data.check(hard=False)

        self.assertEqual(rejected_recs, [2000])
        self.assertEqual(rejected_terms, ["HP:BADCODE"])

        self.assertEqual(len(self.patient_data.profiles.items()), self.n_patients + 1)
        self.assertEqual(len(self.patient_data.vars.items()), self.n_patients + 1)
        self.assertEqual(self.patient_data.profiles[2001], ["HP:0000315"])


    def test_link2ont(self):
        self.patient_data.link2ont(Cohort.act_ont)
        #We check that the profiles of the patients are correctly saved as profiles in the
        #ontology object
        self.assertEqual(Cohort.ont[Cohort.act_ont].profiles.items(),
                         self.patient_data.profiles.items())


    def test_get_profile_redundancy(self):
        pass


    def test_get_profiles_terms_frequency(self):
        pass


    def test_compute_term_list_and_childs(self):
        pass


    def test_get_profile_ontology_distribution_tables(self):
        pass


    def test_get_ic_analysis(self):
        pass


    def test_get_profiles_mean_size(self):
        pass


    def test_get_profile_length_at_percentile(self):
        pass


    def test_get_dataset_specifity_index(self):
        pass


    def test_compare_profiles(self):
        pass


    def test_index_vars(self):
        expected_number_of_regions = sum([len(regions) for pat_gen_feats in self.patient_data.vars.values() 
                                        for chrm, regions in pat_gen_feats.regions.items()])
        
        #Mixing the genomic regions of the patients in a unique Genomic Feature
        self.patient_data.index_vars()
        returned_number_of_regions = sum([len(regions) for chrm, regions in self.patient_data.var_idx.regions.items()])
        #Checking that the number of regions in this mixed genomic feature is the same as the total number of regions of the patients
        self.assertEqual(expected_number_of_regions, returned_number_of_regions)


    def test_get_vars_sizes(self):
        self.maxDiff = None
        expected_sizes = []
        for regions in self.patient_data.vars.values():
            for chrm, region in regions.each():
                    expected_sizes.append(region["stop"] - region["start"] + 1)
        #Mixing the genomic regions of the patients in a unique Genomic Feature (var_idx attribute)
        self.patient_data.index_vars()
        #Getting the sizes of each of the genomic regions of the patients
        returned_sizes = self.patient_data.get_vars_sizes()

        self.assertEqual(sorted(expected_sizes), sorted(returned_sizes))

        #Getting summary sizes, which is a 2D list with the size and number of patients with that size
        returned_sizes = self.patient_data.get_vars_sizes(summary=True)
        expected_sizes = sorted([(size, expected_sizes.count(size)) for size in set(expected_sizes)], 
                                key=lambda x: [x[1],x[0]], reverse=True)
        self.assertEqual(expected_sizes, returned_sizes)


    def test_generate_cluster_regions(self):
        #Mixing the genomic regions of the patients in a unique Genomic Feature
        self.patient_data.index_vars()

        returned_ids_by_cluster_0, returned_annotated_full_red_0 = self.patient_data.var_idx.generate_cluster_regions(meth="reg_overlap", tag="coh", ids_per_reg = 0)
        self.assertEqual(returned_ids_by_cluster_0.keys(), self.patient_data.vars.keys()) #We expect all the patients to be returned when the filtering threshold is 0
        for regions_id in returned_ids_by_cluster_0.values():
            self.assertGreaterEqual(len(regions_id), 1) #We expect at least one region per patient

        for reference_region in returned_annotated_full_red_0:
            self.assertEqual(len(reference_region), 4) #We expect the reference region to have 4 elements: chrm, start, stop and the region_id

        returned_ids_by_cluster_3, returned_annotated_full_red_3 = self.patient_data.var_idx.generate_cluster_regions(meth="reg_overlap", tag="coh", ids_per_reg = 3)
        self.assertLessEqual(len(returned_annotated_full_red_3), len(returned_annotated_full_red_0)) #We expect fewer overlapping regions with more than 3 patients than with 1 patients
        self.assertLessEqual(len(returned_ids_by_cluster_3), len(returned_ids_by_cluster_0)) #The same is true for the patients
        
        for regions_id in returned_ids_by_cluster_3.values():
            self.assertGreaterEqual(len(regions_id), 1) #We expect at least one region per patient, so after the thresholding, patients in the dictionary should have at leat 1 region with more than 3 patients

        for reference_region in returned_annotated_full_red_3:
            self.assertGreaterEqual(int(reference_region[3].split(".")[-1]), 3) #We expect the last number of the region tag (number of patients overlapping the region) to be equal or greater to the threshold we used
            self.assertEqual(reference_region[3].split(".")[0], reference_region[2]) #And also that the first number in the region tag corresponds to the chromosome 

        returned_ids_by_cluster_1000, returned_annotated_full_red_1000 = self.patient_data.var_idx.generate_cluster_regions(meth="reg_overlap", tag="coh", ids_per_reg = 1000)
        self.assertEqual(returned_ids_by_cluster_1000, {}) #We expect no patients to be returned when the filtering threshold is 1000 (not that many overlaps in the toy dataset)
        self.assertEqual(returned_annotated_full_red_1000, []) #Neither we expect any reference regions to be returned


    def test_save(self):
        ### Creating the tmp folder and saving the toy dataset in paco format ###
        tmp_file = os.path.join(ROOT_PATH, "tmp", "test_save")
        self.patient_data.save(f"{tmp_file}", mode="paco", translate=True)

        ### Checking that the file has been created and has the expected number of lines (equals to input file)###
        check_n_lines = f'wc -l {tmp_file} | cut -d " " -f 1'
        ps2 = subprocess.Popen(check_n_lines,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output2 = ps2.communicate()[0]
        total_lines = str(output2).replace("b'","").replace("\\n'", "").strip()
        self.assertEqual(total_lines, "101")

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

        #Because filtering was applied with cohort_parsers the first time we loaded the dataset, 
        # when loading again the dataset the rejected patients and HP terms should have been cleared
        self.assertEqual(len(rejected_hpos_L), 0)
        self.assertEqual(len(rejected_patients_L), 0)

        os.remove(f"{tmp_file}")

    def test_export_phenopackets(self):
        pass