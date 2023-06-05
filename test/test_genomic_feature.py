import sys 
import unittest
import os
from pets import Genomic_Feature
ROOT_PATH=os.path.dirname(__file__)
DATA_TEST_PATH = os.path.join(ROOT_PATH, 'data')

class BaseNetTestCase(unittest.TestCase):
	def setUp(self):
		self.features = [
			['chr1', 5000, 10000],
			['chr1', 11000, 12000],
			['chr2', 501000, 600000],
			['chr3', 2010, 2010],
		] 
		self.genomic_feature = Genomic_Feature(self.features)

		self.named_features = [
			['chr1', 5000, 10000, 'a'],
			['chr1', 11000, 12000, 'b'],
			['chr2', 501000, 600000, 'c'],
			['chr3', 2010, 2010, 'd'],
		]
		self.annotations = {"a": "status:pathogenic;xref:OMIM4345", "c": "status:benign"}
		self.named_genomic_feature = Genomic_Feature(self.named_features, annotations=self.annotations)

		self.repeated_features_size = [
			['chr1', 5000, 5000, 'a'],
			['chr1', 8000, 8000, 'b'],
			['chr2', 9000, 9000, 'c'],
			['chr3', 2010, 2015, 'd'],
			['chr4', 2015, 2020, 'e'],
			['chr5', 1000, 100000, 'e'],			
		]
		self.repeated_genomic_feature_size = Genomic_Feature(self.repeated_features_size)

	def test_array2genomic_feature(self):
		gen_feature_from_array = Genomic_Feature.array2genomic_feature(self.features, lambda r: [r[0], r[1], r[2]])
		named_gen_feature_from_array = Genomic_Feature.array2genomic_feature(self.named_features, lambda r: [r[0], r[1], r[2], r[3]], annotations=self.annotations)

		self.assertEqual(gen_feature_from_array.regions, self.genomic_feature.regions)
		self.assertEqual(named_gen_feature_from_array.regions, self.named_genomic_feature.regions)

	def test_hash2genomic_feature(self):
		gen_feature_from_hash = Genomic_Feature.hash2genomic_feature(self.genomic_feature.reg_by_to, lambda k,v: [v["chrm"], v["start"], v["stop"], v["to"]] )
		self.assertEqual(gen_feature_from_hash.regions, self.genomic_feature.regions)

	def test_add_reference(self):
		Genomic_Feature.add_reference(self.genomic_feature)
		self.assertEqual(Genomic_Feature.ref, self.genomic_feature)

	def test_genomic_feature_attr(self):
		#Checking that length function is corretly working and returning the
		#total number of regions
		self.assertEqual(self.genomic_feature.len(), 4)
		self.assertEqual(self.named_genomic_feature.len(), 4)

	def test_init_loaded_attributes(self):
		#Checking that the regions are correctly loaded in the object attributes
		#for the named features. Both in regions by chromosome and regions by to
		self.assertEqual(list(self.named_genomic_feature.regions.keys()), ['chr1', 'chr2', 'chr3'])

		#Checking that the the chr1 has the expected number of regions and that
		#a region has the expected number of elements
		self.assertEqual(len(self.named_genomic_feature.regions['chr1']), 2)
		self.assertEqual(len(self.named_genomic_feature.regions['chr1'][0]), 5) #5 because of the to attribute
		self.assertEqual(len(self.named_genomic_feature.regions['chr1'][1]), 4) #4 because this region has no additional attributes

		#Checking that the regions by to has the expected number of elements
		self.assertEqual(len(self.named_genomic_feature.reg_by_to.keys()), 4)

		#Checking that the regions are correctly loaded in the object attributes 
		# for the unnamed features. Both in regions by chromosome and regions by to
		self.assertEqual(list(self.genomic_feature.regions.keys()), ['chr1', 'chr2', 'chr3'])

		#Checking that the the chr1 has the expected number of regions and that
		#a region has the expected number of elements
		self.assertEqual(len(self.genomic_feature.regions['chr2']), 1)
		self.assertEqual(len(self.genomic_feature.regions['chr2'][0]), 4)
		
		#Checking that the regions by to has the expected number of elements
		self.assertEqual(len(self.genomic_feature.reg_by_to.keys()), 4)

	def test_each(self):
		#Checking that each method is correctly iterating over the regions
		regions = 0
		for chrm, region in self.genomic_feature.each():
			regions += 1
			self.assertEqual(len(region), 4) #Checking that each region has 4 elements
		self.assertEqual(regions, self.genomic_feature.len())

	def test_each_chr(self):
		#Checking that each_chr is correctly iterating over the chromosomes
		chrms = 0
		for chrm, regs in self.genomic_feature.each_chr():
			chrms += 1
			self.assertEqual(len(regs), len(self.genomic_feature.regions[chrm])) #Checking that each chromosome has the correct number of regions
		self.assertEqual(chrms, len(self.genomic_feature.regions.keys())) #Checking that the number of chromosomes is correct

	def test_get_chr(self):
		self.assertEqual(list(self.genomic_feature.get_chr()), ["chr1", "chr2", "chr3"])

	def test_get_chr_regs(self):
		#checking that chromosome 1 has the expected number of regions
		self.assertEqual(len(self.genomic_feature.get_chr_regs("chr1")), 2)

		#checking that regions are correct
		self.assertEqual(self.genomic_feature.get_chr_regs("chr1"),
			[{"chrm": 'chr1', "start": 5000, "stop": 10000, "to": 0},
			{"chrm": 'chr1', "start": 11000, "stop": 12000, "to": 1}])

		#checking that trying to access a chromosome that does not exist returns None
		self.assertEqual(self.genomic_feature.get_chr_regs("chr4"), None)

	def test_region_by_to(self):
		#Looking for a region using region id (unnamed case) and no attrs
		self.assertEqual(self.genomic_feature.region_by_to(2),
			{"chrm": 'chr2', "start": 501000, "stop": 600000, "to": 2})
		#Looking for a region using region id (named case) with attrs
		self.assertEqual(self.named_genomic_feature.region_by_to("c"),
			{"chrm": 'chr2', "start": 501000, "stop": 600000, "to": "c", 'attrs': 'status:benign'})
		
	def test_get_sizes(self):
		#Returning the sizes of each of the regions
		self.assertEqual(self.genomic_feature.get_sizes(), [5001, 1001, 99001, 1])

	def get_features(self):
		pass

	def test_match(self):
		#Checking that match is correctly returning the regions that overlap with other genomic_feature	
		expected = {0: ["a"], 1: ["b"], 2: ["c"], 3: ["d"]}
		returned = self.genomic_feature.match(self.named_genomic_feature)
		self.assertEqual(expected, returned)
	
	def test_get_summary_sizes(self):
		#It calculate the sizes of each region and return an ordered list
		#of list in which the first element is the size and the second the number
		#of regions with that size
		expected = [(1, 3), (6, 2), (99001, 1)]
		returned = self.repeated_genomic_feature_size.get_summary_sizes()
		self.assertEqual(expected, returned)

	def test_merge(self):
		self.named_genomic_feature.merge(self.genomic_feature)
		#checking that the length of the two merged genomic features is the expected
		self.assertEqual(self.named_genomic_feature.len(), 8)

		#checking that each of the region "to" has been added
		for chrm, region in self.named_genomic_feature.each():
			self.assertIn(region['to'], [4, 5, 6, 7, "a", "b", "c", "d"])
			
			#checking that each region has all the expected attributes
			if len(region.keys()) == 5:
				self.assertEqual(list(region.keys()), ['chrm', 'start', 'stop', 'to', 'attrs'])
			elif len(region.keys()) == 4:
				self.assertEqual(list(region.keys()), ['chrm', 'start', 'stop', 'to'])
			else:
				self.fail("Unexpected number of keys")

	#TODO: Ask about this method (reference is a list of [start,stop] instead of a genomic_feature?)
	def test_get_reference_overlaps(self):
		reference = [[item[1], item[2]] for item in self.features]
		genomic_ranges = [item[1] for item in self.named_genomic_feature.each()]
		self.assertEqual(
			self.genomic_feature.get_reference_overlaps(genomic_ranges, reference),
		   [["a"], ["b"], ["c"], ["d"]])
	
	def test_generate_cluster_regions(self):
		pass

	def test_compute_windows(self):
		expected = {'chr1': [[5000, 10000], [10000, 11000], [11000, 12000]], 'chr2': [[501000, 600000]], 'chr3': [[2010, 2010]]}
		
		self.genomic_feature.compute_windows("reg_overlap")
		self.assertEqual(self.genomic_feature.windows, expected)
		