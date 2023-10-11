import os, unittest 
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.parsers.reference_parser import Reference_parser
from pets.parsers.coord_parser import Coord_Parser
from pets.cohort import Cohort

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
COORDS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, 'coords_toy_dataset.txt'))


class CoordParserTestSuite(unittest.TestCase):
    def setUp(self):
        self.options = {"input_file": COORDS_FILE, "header": True,
                   "id_col":"patient_id", "chromosome_col": "chr", "start_col":"start", "end_col":"end"}

        self.regions_by_to = {  "624": {"start":100963222, "stop":101153990, "chrm":"5", "to":"624"},
                                "625": {"start":102358320, "stop":105487655, "chrm":"7", "to":"625"},
                                "626": {"start":44083882, "stop":44210195, "chrm":"17", "to":"626"},
                                "627": {"start":154208417, "stop":154309447, "chrm":"X", "to":"627"},
                                "628": {"start": 31923988, "stop": 32092796, "chrm":'6', "to": "628"}}
        
        self.regions = { "5":  [{"start":100963222, "stop":101153990, "chrm":"5", "to":"624"}],
                         "6":  [{"start": 31923988, "stop": 32092796, "chrm":'6', "to": "628"}],
                         "7":  [{"start":102358320, "stop":105487655, "chrm":"7", "to":"625"}],
                         "17": [{"start":44083882, "stop":44210195, "chrm":"17", "to":"626"}],
                         "X":  [{"start":154208417, "stop":154309447, "chrm":"X", "to":"627"}]}
    
    def test_load(self):
        gen_features = Coord_Parser.load(self.options)
        self.assertEqual(self.regions_by_to, gen_features.reg_by_to)
        self.assertEqual(self.regions, gen_features.regions)

          

