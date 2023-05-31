import sys 
import unittest
import os
from pets import Genomic_feature
ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')

class BaseNetTestCase(unittest.TestCase):
	def setUp(self):
		features = [
			['chr1', 5000, 10000],
			['chr1', 11000, 12000],
			['chr2', 501000, 600000],
			['chr3', 2010, 2010],
		] 
		self.genomic_feature = Genomic_feature(features)

		named_features = [
			['chr1', 5000, 10000, 'a'],
			['chr1', 11000, 12000, 'b'],
			['chr2', 501000, 600000, 'c'],
			['chr3', 2010, 2010, 'd'],
		] 
		self.named_genomic_feature = Genomic_feature(named_features)

	def test_genomic_feature_attr(self):
		self.assertEqual(self.genomic_feature.len(), 4)
		self.assertEqual(self.named_genomic_feature.len(), 4)
