import sys, os, json, unittest 
from importlib.resources import files

import subprocess
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
PATIENTS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, 'cohort_toy_dataset.txt'))

HPO_FILE = files('pets.external_data').joinpath('hp.json')



class CohortParserTestSuite(unittest.TestCase):
    def setUp(self):
        self.hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
        #Cohort.load_ontology("hpo", self.hpo_file, "./test/data/excluded.txt") #Excluded phenotype is Obesity HP:0001513, that is present in patient 625 and 627. 
        ##Line above commented: it seems to only affect semtools related methods, not the parser itself 
        Cohort.load_ontology("hpo", self.hpo_file )
        Cohort.act_ont = "hpo"
       
        self.n_patients = 4
        self.patient_624 = {"regions": {0:{"start":100963222, "stop":101153990, "chrm":"5", "to":0}}, "phenotypes": ["HP:0001249"], "sex":"M"}
        self.patient_625 = {"regions": {0:{"start":102358320, "stop":105487655, "chrm":"7", "to":0}}, "phenotypes": ["HP:0001249", "HP:0001513"], "sex":"M"}
        self.patient_626 = {"regions": {0:{"start":31923988, "stop":32092796, "chrm":"6", "to":0},
                                        1:{"start":44083882, "stop":44210195, "chrm":"17", "to":1}}, "phenotypes": ["HP:0001249"], "sex":"M"}
        self.patient_627 = {"regions": {0:{"start":154208417, "stop":154309447, "chrm":"X", "to":0}}, "phenotypes": ["HP:0001249", "HP:0000929", "HP:0010461", "HP:0001513"], "sex":"F"}

        self.patients = {"624": self.patient_624, "625": self.patient_625, "626": self.patient_626, "627": self.patient_627}

        self.options = {"input_file": PATIENTS_FILE, "header": True, "separator": "|", "names": True,
                   "id_col":"patient_id", "chromosome_col": "chr", "start_col":"start", "end_col":"end", "ont_col":"phenotypes", "sex_col": "sex"}    


    def test_load(self):
        self.patient_data, self.rejected_hpos_L, self.rejected_patients_L = Cohort_Parser.load(self.options)
        self.assertEqual(self.rejected_hpos_L, [])
        self.assertEqual(self.rejected_patients_L, [])

        for patientID in ["624", "625", "626", "627"]:
            self.assertEqual(sorted(self.patients[patientID]["phenotypes"]), sorted(self.patient_data.profiles[patientID]))
            self.assertDictEqual(self.patients[patientID]["regions"], self.patient_data.vars[patientID].reg_by_to)
            self.assertEqual(self.patients[patientID]["sex"], self.patient_data.extra_attr[patientID]["sex"])


    def test_load_wrong_patient_data_and_patient_with_alternative_hpo_name(self):
        self.options.update({"input_file": os.path.join(DATA_TEST_PATH, "cohort_toy_dataset_with_wrong_patient_data.txt")})
        self.patient_data, self.rejected_hpos_L, self.rejected_patients_L = Cohort_Parser.load(self.options)
        self.assertEqual(self.rejected_hpos_L, ["Wrong Phenotype"])
        self.assertEqual(self.rejected_patients_L, ["625"])
    
        #This time "HP:0001249" was described ad Mental Retardation in the input file instead of "Intellectual Disability"
        self.assertEqual(sorted(self.patients["624"]["phenotypes"]), sorted(self.patient_data.profiles["624"]))
        self.assertDictEqual(self.patients["624"]["regions"], self.patient_data.vars["624"].reg_by_to)
        self.assertEqual(self.patients["624"]["sex"], self.patient_data.extra_attr["624"]["sex"])
          
