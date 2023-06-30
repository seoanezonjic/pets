import os, unittest 
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.parsers.reference_parser import Reference_parser
from pets.parsers.coord_parser import Coord_Parser
from pets.cohort import Cohort

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
PATIENTS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, '100_test_dataset_with_sex.txt'))


class CoordParserTestSuite(unittest.TestCase):
    def setUp(self):
        self.options = {"input_file": PATIENTS_FILE,
                   "id_col":"patient_id", "chromosome_col": "chr", 
                   "header": True, "separator": "|", "names": True,
                   "start_col":"start", "end_col":"end", "ont_col":"phenotypes",
                   "sex_col": "sex"}


    def test_load(self):
        chroms = set([str(chrm) for chrm in range(1,23)] + ['X', 'Y'])
        fields = set(["start", "stop", "to", "chrm"])

        gen_features = Coord_Parser.load(self.options)
        # It should be a Genomic_Feature object
        self.assertIsInstance(gen_features, Genomic_Feature)

        # It should have 84 elements and 24 chromosomes (1-22, X and Y) and the fields start, stop, to and chrm
        self.assertEqual(len(gen_features.reg_by_to),84)
        self.assertEqual(set(gen_features.regions.keys()), chroms)
        self.assertEqual(set(gen_features.reg_by_to["763"].keys()), fields)

        self.assertEqual(type(gen_features.reg_by_to["763"]["start"]), int)
        self.assertEqual(type(gen_features.reg_by_to["763"]["stop"]), int)
        self.assertEqual(type(gen_features.reg_by_to["763"]["to"]), str)
        self.assertEqual(type(gen_features.reg_by_to["763"]["chrm"]), str)
