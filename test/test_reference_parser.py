import os, unittest 
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.parsers.reference_parser import Reference_parser
from pets.cohort import Cohort

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
PATIENTS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, '100_test_dataset_with_sex.txt'))
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))


class RefParserTestSuite(unittest.TestCase):
    def setUp(self):
        self.gtf_file = os.path.join(DATA_TEST_PATH, "ref_parser", "toy_example.gtf")


    def test_load(self):
        ref_all = Reference_parser.load( file_path=self.gtf_file)
        ref_exons = Reference_parser.load( file_path=self.gtf_file, feature_type="exon")
        ref_transcripts = Reference_parser.load( file_path=self.gtf_file, feature_type="transcript")

        # all variables should be Genomic_Feature objects
        self.assertIsInstance(ref_all, Genomic_Feature)
        self.assertIsInstance(ref_exons, Genomic_Feature)
        self.assertIsInstance(ref_transcripts, Genomic_Feature)

        # Checking that ref_all has 5 elements, ref exon 3 and ref transcript 2
        self.assertEqual(len(ref_all.reg_by_to), 5)
        self.assertEqual(len(ref_exons.reg_by_to), 3)
        self.assertEqual(len(ref_transcripts.reg_by_to), 2)
