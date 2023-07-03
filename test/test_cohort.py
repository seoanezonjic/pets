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
        self.patient_data.link2ont(Cohort.act_ont)
        profile_sizes, parental_terms_per_profile = self.patient_data.get_profile_redundancy()
        print(profile_sizes)
        print(parental_terms_per_profile)


#    def test_get_profiles_terms_frequency(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        term_stats = self.patient_data.get_profiles_terms_frequency()
#        print(term_stats)
#        self.assertEqual(term_stats, open(os.path.join(DATA_TEST_PATH,"profiles_terms_frequency.txt")).read().strip())
#
#
#    def test_compute_term_list_and_childs(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        #suggested_childs, term_with_childs_ratio = self.patient_data.compute_term_list_and_childs()
#        #print(suggested_childs)
#        #print(term_with_childs_ratio)
#
#
#    def test_get_profile_ontology_distribution_tables(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        #ontology_levels, distribution_percentage = self.patient_data.get_profile_ontology_distribution_tables()
#        #print(ontology_levels)
#        #print(distribution_percentage)
#        #TODO: Fix
#        #Traceback (most recent call last):
#        #  File "/mnt/home/users/bio_267_uma/jperezg/dev_py/pets/test/test_cohort.py", line 249, in test_get_profile_ontology_distribution_tables
#        #    ontology_levels, distribution_percentage = self.patient_data.get_profile_ontology_distribution_tables()
#        #  File "/mnt/home/users/bio_267_uma/jperezg/dev_py/pets/pets/cohort.py", line 170, in get_profile_ontology_distribution_tables
#        #    ontology_levels, distribution_percentage = ont.get_profile_ontology_distribution_tables()
#        #  File "/mnt/home/soft/soft_bio_267/programs/x86_64/pyenv/.pyenv/versions/3.10.8/lib/python3.10/site-packages/py_semtools/ontology.py", line 1104, in get_profile_ontology_distribution_tables
#        #    cohort_ontology_levels = self.get_ontology_levels_from_profiles(uniq=False)
#        #  File "/mnt/home/soft/soft_bio_267/programs/x86_64/pyenv/.pyenv/versions/3.10.8/lib/python3.10/site-packages/py_semtools/ontology.py", line 1092, in get_ontology_levels_from_profiles
#        #    terms_levels = self.dicts['level']['byValue']
#        #KeyError: 'level'
#
#
#    def test_get_ic_analysis(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        onto_ic, freq_ic, onto_ic_profile, freq_ic_profile = self.patient_data.get_ic_analysis()
#        self.assertEqual(onto_ic, {'HP:0002236': 4.2363356859539065, 'HP:0000189': 3.9353056902899253, 'HP:0000426': 4.2363356859539065, 'HP:0001320': 3.759214431234244, 'HP:0010492': 4.2363356859539065, 'HP:0008070': 3.1571544399062814, 'HP:0001385': 4.2363356859539065, 'HP:0001089': 4.2363356859539065, 'HP:0002078': 3.9353056902899253, 'HP:0000365': 2.6681339618869115, 'HP:0010465': 4.2363356859539065, 'HP:0002373': 3.634275694625944, 'HP:0004467': 4.2363356859539065, 'HP:0001572': 3.759214431234244, 'HP:0001252': 3.0902076502756683, 'HP:0006695': 3.3912376459396496, 'HP:0002121': 3.4581844355702627, 'HP:0000154': 4.2363356859539065, 'HP:0001344': 4.2363356859539065, 'HP:0000286': 3.537365681617888, 'HP:0011567': 4.2363356859539065, 'HP:0010772': 3.3332456989619628, 'HP:0001642': 3.634275694625944, 'HP:0006934': 3.759214431234244, 'HP:0000194': 3.9353056902899253, 'HP:0004279': 4.2363356859539065, 'HP:0000463': 4.2363356859539065, 'HP:0000864': 2.430155711970019, 'HP:0001540': 4.2363356859539065, 'HP:0000176': 4.2363356859539065, 'HP:0010818': 4.2363356859539065, 'HP:0007328': 3.9353056902899253, 'HP:0004691': 4.2363356859539065, 'HP:0000219': 4.2363356859539065, 'HP:0000960': 4.2363356859539065, 'HP:0011304': 3.634275694625944, 'HP:0000414': 4.2363356859539065, 'HP:0000262': 3.759214431234244, 'HP:0001156': 3.1223923336470696, 'HP:0002553': 4.2363356859539065, 'HP:0000444': 4.2363356859539065, 'HP:0000687': 3.9353056902899253, 'HP:0001339': 2.9353056902899253, 'HP:0000615': 3.2820931765145813, 'HP:0000028': 3.759214431234244, 'HP:0001608': 2.8561244442423, 'HP:0001249': 3.3912376459396496, 'HP:0001260': 3.634275694625944, 'HP:0001305': 4.2363356859539065, 'HP:0000151': 4.2363356859539065, 'HP:0000276': 4.2363356859539065, 'HP:0001629': 3.2363356859539065, 'HP:0000280': 3.9353056902899253, 'HP:0009909': 4.2363356859539065, 'HP:0002186': 3.3332456989619628, 'HP:0000343': 4.2363356859539065, 'HP:0000465': 4.2363356859539065, 'HP:0000717': 4.2363356859539065, 'HP:0000494': 4.2363356859539065, 'HP:0000821': 3.3912376459396496, 'HP:0001250': 1.6947564420073253, 'HP:0000411': 4.2363356859539065, 'HP:0000437': 4.2363356859539065, 'HP:0001388': 3.759214431234244, 'HP:0000581': 3.759214431234244, 'HP:0000708': 1.7463772065290717, 'HP:0000486': 2.436995136500325, 'HP:0010773': 3.9353056902899253, 'HP:0000752': 3.9353056902899253, 'HP:0000402': 3.634275694625944, 'HP:0000664': 4.2363356859539065, 'HP:0000023': 4.2363356859539065, 'HP:0004322': 2.8383956772818686, 'HP:0001636': 3.3912376459396496, 'HP:0000252': 3.537365681617888, 'HP:0001611': 4.2363356859539065, 'HP:0002079': 3.9353056902899253, 'HP:0000545': 3.4581844355702627, 'HP:0001363': 3.060244426898225, 'HP:0000218': 3.9353056902899253, 'HP:0002571': 4.2363356859539065, 'HP:0001822': 4.2363356859539065, 'HP:0002750': 3.759214431234244, 'HP:0000272': 4.2363356859539065, 'HP:0100335': 3.194943000795681, 'HP:0000776': 3.4581844355702627, 'HP:0000479': 1.8765002036140184, 'HP:0004209': 3.9353056902899253, 'HP:0000058': 3.060244426898225, 'HP:0000490': 4.2363356859539065, 'HP:0000089': 3.759214431234244, 'HP:0000582': 4.2363356859539065, 'HP:0000322': 4.2363356859539065, 'HP:0011667': 4.2363356859539065, 'HP:0000256': 3.537365681617888, 'HP:0000243': 4.2363356859539065, 'HP:0009803': 2.8746078499363135, 'HP:0003196': 4.2363356859539065, 'HP:0001631': 3.3912376459396496, 'HP:0000294': 4.2363356859539065, 'HP:0000076': 3.194943000795681, 'HP:0010461': 2.1124840449868207, 'HP:0010538': 4.2363356859539065, 'HP:0001513': 3.3912376459396496, 'HP:0000193': 4.2363356859539065, 'HP:0000369': 3.9353056902899253, 'HP:0003508': 3.634275694625944, 'HP:0000736': 3.9353056902899253, 'HP:0000054': 4.2363356859539065, 'HP:0009894': 4.2363356859539065, 'HP:0009908': 4.2363356859539065, 'HP:0000215': 4.2363356859539065, 'HP:0002561': 4.2363356859539065, 'HP:0000238': 3.537365681617888, 'HP:0010721': 3.634275694625944, 'HP:0000786': 4.2363356859539065, 'HP:0000518': 2.592883009467719, 'HP:0000047': 3.3332456989619628, 'HP:0000601': 4.2363356859539065, 'HP:0001251': 2.9353056902899253, 'HP:0000232': 4.2363356859539065, 'HP:0008872': 3.4581844355702627, 'HP:0000135': 3.537365681617888, 'HP:0006335': 4.2363356859539065, 'HP:0000316': 4.2363356859539065, 'HP:0000472': 4.2363356859539065, 'HP:0008734': 4.2363356859539065, 'HP:0000358': 3.9353056902899253, 'HP:0000750': 3.2363356859539065, 'HP:0000160': 4.2363356859539065, 'HP:0000534': 2.8561244442423, 'HP:0008551': 3.634275694625944, 'HP:0000311': 4.2363356859539065, 'HP:0000929': 1.784549250429616, 'HP:0005326': 4.2363356859539065, 'HP:0009889': 3.537365681617888, 'HP:0000954': 3.9353056902899253, 'HP:0010055': 3.9353056902899253, 'HP:0001371': 2.2632078323542077, 'HP:0002360': 2.804971921794919, 'HP:0000470': 4.2363356859539065, 'HP:0008034': 3.1571544399062814, 'HP:0000448': 4.2363356859539065, 'HP:0002650': 3.194943000795681, 'HP:0000974': 3.634275694625944, 'HP:0000077': 1.5031384208473368, 'HP:0000303': 4.2363356859539065, 'HP:0007766': 4.2363356859539065, 'HP:0001176': 4.2363356859539065, 'HP:0000337': 4.2363356859539065, 'HP:0001763': 3.9353056902899253, 'HP:0001609': 3.9353056902899253, 'HP:0001199': 3.634275694625944, 'HP:0009748': 4.2363356859539065, 'HP:0001943': 3.2363356859539065, 'HP:0008069': 2.3730128258334506, 'HP:0004305': 2.4959729964596624, 'HP:0000574': 4.2363356859539065, 'HP:0000508': 3.3332456989619628, 'HP:0000527': 3.634275694625944, 'HP:0000455': 4.2363356859539065, 'HP:0000098': 3.4581844355702627, 'HP:0001833': 4.2363356859539065, 'HP:0001518': 4.2363356859539065, 'HP:0000248': 3.759214431234244, 'HP:0010469': 4.2363356859539065, 'HP:0000341': 4.2363356859539065, 'HP:0005280': 4.2363356859539065, 'HP:0008065': 2.8939130051317, 'HP:0010847': 4.2363356859539065, 'HP:0000431': 4.2363356859539065, 'HP:0001643': 3.634275694625944, 'HP:0000957': 3.759214431234244, 'HP:0000767': 3.9353056902899253})
#        self.assertEqual(freq_ic, {'HP:0002236': 2.649334858712142, 'HP:0000189': 2.3483048630481607, 'HP:0000426': 2.3483048630481607, 'HP:0001320': 2.3483048630481607, 'HP:0010492': 2.649334858712142, 'HP:0008070': 2.649334858712142, 'HP:0001385': 2.649334858712142, 'HP:0001089': 2.649334858712142, 'HP:0002078': 2.649334858712142, 'HP:0000365': 2.0472748673841794, 'HP:0010465': 2.649334858712142, 'HP:0002373': 2.649334858712142, 'HP:0004467': 2.649334858712142, 'HP:0001572': 2.649334858712142, 'HP:0001252': 1.695092349272817, 'HP:0006695': 2.1722136039924793, 'HP:0002121': 2.3483048630481607, 'HP:0000154': 2.3483048630481607, 'HP:0001344': 2.649334858712142, 'HP:0000286': 2.0472748673841794, 'HP:0011567': 2.649334858712142, 'HP:0010772': 2.3483048630481607, 'HP:0001642': 2.1722136039924793, 'HP:0006934': 2.649334858712142, 'HP:0000194': 2.649334858712142, 'HP:0004279': 2.649334858712142, 'HP:0000463': 2.1722136039924793, 'HP:0000864': 2.649334858712142, 'HP:0001540': 2.649334858712142, 'HP:0000176': 2.649334858712142, 'HP:0010818': 2.649334858712142, 'HP:0007328': 2.649334858712142, 'HP:0004691': 2.1722136039924793, 'HP:0000219': 2.0472748673841794, 'HP:0000960': 2.3483048630481607, 'HP:0011304': 2.3483048630481607, 'HP:0000414': 2.649334858712142, 'HP:0000262': 2.3483048630481607, 'HP:0001156': 2.649334858712142, 'HP:0002553': 2.649334858712142, 'HP:0000444': 2.649334858712142, 'HP:0000687': 2.3483048630481607, 'HP:0001339': 2.649334858712142, 'HP:0000615': 2.649334858712142, 'HP:0000028': 1.7462448717201982, 'HP:0001608': 2.0472748673841794, 'HP:0001249': 0.76284413353966, 'HP:0001260': 2.649334858712142, 'HP:0001305': 2.649334858712142, 'HP:0000151': 2.649334858712142, 'HP:0000276': 2.3483048630481607, 'HP:0001629': 2.3483048630481607, 'HP:0000280': 2.649334858712142, 'HP:0009909': 2.649334858712142, 'HP:0002186': 2.649334858712142, 'HP:0000343': 2.1722136039924793, 'HP:0000465': 2.649334858712142, 'HP:0000717': 1.4452148760562171, 'HP:0000494': 2.1722136039924793, 'HP:0000821': 2.649334858712142, 'HP:0001250': 1.695092349272817, 'HP:0000411': 2.3483048630481607, 'HP:0000437': 2.649334858712142, 'HP:0001388': 1.8711836083284983, 'HP:0000581': 2.3483048630481607, 'HP:0000708': 1.0365510019924065, 'HP:0000486': 1.7462448717201982, 'HP:0010773': 2.649334858712142, 'HP:0000752': 2.0472748673841794, 'HP:0000402': 2.649334858712142, 'HP:0000664': 2.649334858712142, 'HP:0000023': 2.649334858712142, 'HP:0004322': 1.3069121778899357, 'HP:0001636': 2.649334858712142, 'HP:0000252': 1.649334858712142, 'HP:0001611': 2.3483048630481607, 'HP:0002079': 2.3483048630481607, 'HP:0000545': 2.649334858712142, 'HP:0001363': 2.649334858712142, 'HP:0000218': 2.1722136039924793, 'HP:0002571': 2.649334858712142, 'HP:0001822': 2.649334858712142, 'HP:0002750': 2.649334858712142, 'HP:0000272': 2.649334858712142, 'HP:0100335': 2.649334858712142, 'HP:0000776': 2.649334858712142, 'HP:0000479': 2.649334858712142, 'HP:0004209': 2.1722136039924793, 'HP:0000058': 2.649334858712142, 'HP:0000490': 2.3483048630481607, 'HP:0000089': 2.649334858712142, 'HP:0000582': 2.0472748673841794, 'HP:0000322': 2.3483048630481607, 'HP:0011667': 2.649334858712142, 'HP:0000256': 2.649334858712142, 'HP:0000243': 2.3483048630481607, 'HP:0009803': 2.649334858712142, 'HP:0003196': 2.3483048630481607, 'HP:0001631': 2.1722136039924793, 'HP:0000294': 2.649334858712142, 'HP:0000076': 2.649334858712142, 'HP:0010461': 1.3271155639782226, 'HP:0010538': 2.649334858712142, 'HP:0001513': 2.1722136039924793, 'HP:0000193': 2.1722136039924793, 'HP:0000369': 2.0472748673841794, 'HP:0003508': 2.0472748673841794, 'HP:0000736': 2.649334858712142, 'HP:0000054': 2.0472748673841794, 'HP:0009894': 2.649334858712142, 'HP:0009908': 2.649334858712142, 'HP:0000215': 2.649334858712142, 'HP:0002561': 2.649334858712142, 'HP:0000238': 2.3483048630481607, 'HP:0010721': 2.1722136039924793, 'HP:0000786': 2.649334858712142, 'HP:0000518': 2.1722136039924793, 'HP:0000047': 2.0472748673841794, 'HP:0000601': 2.649334858712142, 'HP:0001251': 2.3483048630481607, 'HP:0000232': 2.3483048630481607, 'HP:0008872': 1.570153612664517, 'HP:0000135': 2.649334858712142, 'HP:0006335': 2.649334858712142, 'HP:0000316': 1.8711836083284983, 'HP:0000472': 2.649334858712142, 'HP:0008734': 2.649334858712142, 'HP:0000358': 2.649334858712142, 'HP:0000750': 1.5353915064053052, 'HP:0000160': 2.649334858712142, 'HP:0000534': 1.9503648543761232, 'HP:0008551': 2.3483048630481607, 'HP:0000311': 2.649334858712142, 'HP:0000929': 1.2876070226945489, 'HP:0005326': 2.3483048630481607, 'HP:0009889': 2.649334858712142, 'HP:0000954': 2.649334858712142, 'HP:0010055': 2.649334858712142, 'HP:0001371': 2.3483048630481607, 'HP:0002360': 2.1722136039924793, 'HP:0000470': 2.3483048630481607, 'HP:0008034': 2.649334858712142, 'HP:0000448': 2.649334858712142, 'HP:0002650': 1.8711836083284983, 'HP:0000974': 2.649334858712142, 'HP:0000077': 2.3483048630481607, 'HP:0000303': 2.649334858712142, 'HP:0007766': 2.649334858712142, 'HP:0001176': 2.649334858712142, 'HP:0000337': 2.649334858712142, 'HP:0001763': 2.649334858712142, 'HP:0001609': 2.649334858712142, 'HP:0001199': 2.649334858712142, 'HP:0009748': 2.649334858712142, 'HP:0001943': 2.3483048630481607, 'HP:0008069': 2.649334858712142, 'HP:0004305': 2.649334858712142, 'HP:0000574': 2.3483048630481607, 'HP:0000508': 2.649334858712142, 'HP:0000527': 2.649334858712142, 'HP:0000455': 2.3483048630481607, 'HP:0000098': 2.1722136039924793, 'HP:0001833': 2.649334858712142, 'HP:0001518': 2.3483048630481607, 'HP:0000248': 2.649334858712142, 'HP:0010469': 2.3483048630481607, 'HP:0000341': 2.649334858712142, 'HP:0005280': 2.649334858712142, 'HP:0008065': 2.649334858712142, 'HP:0010847': 2.649334858712142, 'HP:0000431': 2.3483048630481607, 'HP:0001643': 2.3483048630481607, 'HP:0000957': 2.649334858712142, 'HP:0000767': 2.3483048630481607})
#        self.assertEqual(onto_ic_profile,{'105': 3.780058831652254, '106': 3.622256642870539, '126': 3.382798109980759, '127': 3.926541472047526, '128': 3.5450787238346773, '129': 3.9362827509371123, '130': 3.4370223497298937, '131': 3.3690143350049926, '132': 3.5725936607230744, '133': 3.1652280024927766, '135': 3.3958520114953212, '138': 2.8945563702149544, '139': 3.822302231318444, '140': 3.8865046305749464, '141': 3.6432345441307317, '147': 3.1399865912663416, '184': 3.6003186748663376, '221': 3.789632191528686, '306': 3.4892668153410775, '308': 4.115923687688314, '309': 3.928417561349144, '310': 3.698784138999596, '311': 3.6213030059491538, '313': 3.813786665946778, '314': 3.6280715081381936, '318': 3.736950749369396, '320': 3.2257410364189374, '321': 3.4255761737823325, '323': 3.582995292577089, '523': 3.8264920814198704, '578': 3.7336870235701647, '589': 3.5362121709147183, '595': 3.477515084562825, '596': 2.9651976935110156, '597': 2.947685572208666, '599': 3.3912376459396496, '603': 3.2556663349464685, '605': 2.992060327703684, '607': 3.818581502230781, '608': 3.0624051360741986, '610': 3.079412227427331, '611': 3.600318674866337, '613': 3.5752260385869468, '614': 3.488656336391808, '615': 3.699803677034347, '618': 3.813786665946778, '622': 3.143378479895245, '624': 3.3912376459396496, '625': 3.3912376459396496, '626': 3.3912376459396496, '627': 2.669877146823934, '630': 3.512756670282797, '631': 3.3632154334216873, '632': 3.341588973466233, '633': 3.3912376459396496, '634': 3.1694137159874716, '635': 3.339603659277735, '636': 3.1632716681147874, '637': 3.1148166616107593, '638': 3.240722648107659, '639': 3.488656336391808, '640': 3.8542930073944937, '641': 3.7955959210425996, '644': 3.1148166616107593, '645': 3.3944501175737622, '646': 3.313786665946778, '647': 3.5752260385869468, '648': 3.198912429814816, '649': 3.413553242483187, '652': 3.3912376459396496, '748': 3.949803677034347, '749': 3.31290961095252, '750': 4.025061175950342, '752': 3.79912747088759, '753': 3.4930075406624366, '754': 3.721233901315839, '757': 2.823620183810228, '758': 3.5502212457306084, '759': 3.3844879598362234, '761': 3.3619192558212725, '762': 3.7710556058609663, '763': 4.084985440524541, '766': 3.72126132546284, '767': 3.405240028951635})
#        self.assertEqual(freq_ic_profile, open(os.path.join(DATA_TEST_PATH,"freq_ic_profile.txt")).read().strip())
#
#
#    def test_get_profiles_mean_size(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        profile_mean_size = self.patient_data.get_profiles_mean_size()
#        self.assertEqual(profile_mean_size, 5.3095)
#
#
#    def test_get_profile_length_at_percentile(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        length_percent = self.patient_data.get_profile_length_at_percentile()
#        self.assertEqual(length_percent, 4)
#
#
#    def test_get_dataset_specifity_index(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        #dsi_uniq = self.patient_data.get_dataset_specifity_index("uniq")
#        #dsi_weigthed = self.patient_data.get_dataset_specifity_index("weigthed")
#        #print(dsi_uniq, dsi_weigthed)
#        #Traceback (most recent call last):
#        #  File "/mnt/home/users/bio_267_uma/jperezg/dev_py/pets/test/test_cohort.py", line 288, in test_get_dataset_specifity_index
#        #    dsi_uniq = self.patient_data.get_dataset_specifity_index("uniq")
#        #  File "/mnt/home/users/bio_267_uma/jperezg/dev_py/pets/pets/cohort.py", line 193, in get_dataset_specifity_index
#        #    dsi = ont.get_dataset_specifity_index(type)
#        #  File "/mnt/home/soft/soft_bio_267/programs/x86_64/pyenv/.pyenv/versions/3.10.8/lib/python3.10/site-packages/py_semtools/ontology.py", line 1147, in get_dataset_specifity_index
#        #    ontology_levels, distribution_percentage = self.get_profile_ontology_distribution_tables()
#        #  File "/mnt/home/soft/soft_bio_267/programs/x86_64/pyenv/.pyenv/versions/3.10.8/lib/python3.10/site-packages/py_semtools/ontology.py", line 1104, in get_profile_ontology_distribution_tables
#        #    cohort_ontology_levels = self.get_ontology_levels_from_profiles(uniq=False)
#        #  File "/mnt/home/soft/soft_bio_267/programs/x86_64/pyenv/.pyenv/versions/3.10.8/lib/python3.10/site-packages/py_semtools/ontology.py", line 1092, in get_ontology_levels_from_profiles
#        #    terms_levels = self.dicts['level']['byValue']
#        #KeyError: 'level'
#
#    def test_compare_profiles(self):
#        self.patient_data.link2ont(Cohort.act_ont)
#        similarities = self.patient_data.compare_profiles()
#        self.assertEqual(f"{similarities}", open(os.path.join(DATA_TEST_PATH,"similarities.txt")).read().strip())


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
        os.makedirs(os.path.join(ROOT_PATH, "tmp"), exist_ok=True)
        tmp_file = os.path.join(ROOT_PATH, "tmp", "test_paco_format.txt")
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

        self.assertEqual(total_files, len(self.patient_data.profiles.keys()))
        
        #for id in self.patient_data.profiles.keys():
        #    pheno_id_file = os.path.join(tmp_file, id + ".json")
        #    os.remove(f"{pheno_id_file}")
