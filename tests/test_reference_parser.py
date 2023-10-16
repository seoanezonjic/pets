import os, unittest 
from pets.genomic_features import Genomic_Feature
from pets.parsers.cohort_parser import Cohort_Parser
from pets.parsers.reference_parser import Reference_parser
from pets.cohort import Cohort

ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')
PATIENTS_FILE = os.path.abspath(os.path.join(DATA_TEST_PATH, '100_test_dataset_with_sex.txt'))



class RefParserTestSuite(unittest.TestCase):
    def setUp(self):
        self.gtf_file = os.path.join(DATA_TEST_PATH, "ref_parser", "toy_example.gtf")
        self.gtf_gz = os.path.join(DATA_TEST_PATH, "ref_parser", "toy_example_compr.gtf.gz")

        attrs_dict = {"ENSG00000170006.12": {"feature": "gene", "source": "HAVANA", "gene_id": "ENSG00000170006.12", "gene_type": "protein_coding", "gene_name": "TMEM154", "level": "2", "hgnc_id": "HGNC:26489", "havana_gene": "OTTHUMG00000161463.2"},
                      "ENSG00000248571.1":  {"feature": "gene", "source": "HAVANA", "gene_id": "ENSG00000248571.1", "gene_type": "lncRNA", "gene_name": "AC106882.1", "level": "2", "havana_gene": "OTTHUMG00000161462.1"},
                      "ENSG00000169989.3":  {"feature": "gene", "source": "HAVANA", "gene_id": "ENSG00000169989.3", "gene_type": "protein_coding", "gene_name": "TIGD4", "level": "2", "hgnc_id": "HGNC:18335", "havana_gene": "OTTHUMG00000161464.2"},
                      "ENSG00000164144.16": {"feature": "gene", "source": "HAVANA", "gene_id": "ENSG00000164144.16", "gene_type": "protein_coding", "gene_name": "ARFIP1", "level": "2", "hgnc_id": "HGNC:21496", "havana_gene": "OTTHUMG00000161468.3"},
                      "ENSG00000250980.1":  {"feature": "gene", "source": "HAVANA", "gene_id": "ENSG00000250980.1", "gene_type": "processed_pseudogene", "gene_name": "NSA2P6", "level": "1", "hgnc_id": "HGNC:54587", "tag": "pseudo_consens", "havana_gene": "OTTHUMG00000161467.1"}
        }

        self.regions = {"4": {"ENSG00000170006.12": {"chrm":"4","start":152618628,"stop":152680012,"to":"ENSG00000170006.12", "attrs": attrs_dict["ENSG00000170006.12"]},
                         "ENSG00000248571.1":  {"chrm":"4","start":152666368,"stop":152670107,"to":"ENSG00000248571.1", "attrs": attrs_dict["ENSG00000248571.1"]},
                         "ENSG00000169989.3":  {"chrm":"4","start":152769354,"stop":152779730,"to":"ENSG00000169989.3", "attrs": attrs_dict["ENSG00000169989.3"]},
                         "ENSG00000164144.16": {"chrm":"4","start":152779937,"stop":152918463,"to":"ENSG00000164144.16", "attrs": attrs_dict["ENSG00000164144.16"]},
                         "ENSG00000250980.1":  {"chrm":"4","start":152796063,"stop":152796831,"to":"ENSG00000250980.1", "attrs": attrs_dict["ENSG00000250980.1"]}
                    }}

    def test_load(self):
        reference_chunk = Reference_parser.load( file_path=self.gtf_file, feature_type="gene")

        self.assertIsInstance(reference_chunk, Genomic_Feature) #Asserting that the resulting object is a Genomic_Feature
        self.assertEqual(len(reference_chunk.reg_by_to), 5) #Asserting that we have 5 genes
        self.assertEqual(list(reference_chunk.regions.keys()), ["4"]) #Asserting all the regions are from chromosome 4

        for gene in ["ENSG00000170006.12", "ENSG00000248571.1", "ENSG00000169989.3", "ENSG00000164144.16", "ENSG00000250980.1"]:
            self.assertEqual(self.regions["4"][gene], reference_chunk.reg_by_to[gene])
        

    def test_load_gz(self):
        reference_chunk = Reference_parser.load( file_path=self.gtf_gz, feature_type="gene")

        self.assertIsInstance(reference_chunk, Genomic_Feature) #Asserting that the resulting object is a Genomic_Feature
        self.assertEqual(len(reference_chunk.reg_by_to), 5) #Asserting that we have 5 genes
        self.assertEqual(list(reference_chunk.regions.keys()), ["4"]) #Asserting all the regions are from chromosome 4

        for gene in ["ENSG00000170006.12", "ENSG00000248571.1", "ENSG00000169989.3", "ENSG00000164144.16", "ENSG00000250980.1"]:
            self.assertEqual(self.regions["4"][gene], reference_chunk.reg_by_to[gene])