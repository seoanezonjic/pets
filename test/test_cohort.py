import sys 
import unittest
import os
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

    def test_load(self):
        n_patients = 84
        chroms = [str(chrm) for chrm in range(1,23)] + ['X', 'Y']
        
        self.assertEqual(sorted(self.rejected_hpos_L),
                         sorted(['Generalized tonic-clonic seizures', 'Stereotypy', 'Sparse and thin eyebrow', 'obsolete Prenatal short stature', 'Abnormality of the pinna', 'Capillary hemangiomas']))
        self.assertEqual(self.rejected_patients_L, [])
        
        self.assertEqual(len(self.patient_data.profiles.keys()), 84) #Checking that all patients are loaded with their phenotypes
        self.assertEqual(len(self.patient_data.vars.keys()), 84) #Checking that all patients are loaded with their genomic regions


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

        self.assertGreaterEqual(84, len(self.patient_data.profiles.keys())) #There are 84 patients in the toy dataset before the filtering
