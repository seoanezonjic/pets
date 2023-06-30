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
        hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
        Cohort.load_ontology("hpo", hpo_file, "./test/data/excluded.txt")
        #Cohort.load_ontology("hpo", hpo_file )
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