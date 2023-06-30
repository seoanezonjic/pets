#! /usr/bin/env python
import argparse, os, re, sys

ROOT_PATH=os.path.dirname(__file__)
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

from pets.genomic_features import Genomic_Feature
from pets.parsers.reference_parser import Reference_parser
from pets.parsers.coord_parser import Coord_Parser

############################################################################################
## OPTPARSE
############################################################################################
parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')

parser.add_argument("-c", "--chromosome_col", dest="chromosome_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the chromosome")

parser.add_argument("-d", "--id_col", dest="id_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the id")

parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")

#chr\tstart\tstop
parser.add_argument("-H", "--header", dest="header", default= True, action="store_false",
                    help="Set if the file has a line header. Default true")

parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                    help="Input file path")

parser.add_argument("-r", "--reference_file", dest="reference_file", default= None,
                    help="Reference file with genome annotation")

parser.add_argument("-o", "--output_file", dest="output_file", default= None,
                    help="Output file with patient data")

parser.add_argument("-s", "--start_col", dest="start_col", default= None,
                    help="Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate")

parser.add_argument("-t", "--feature_type", dest="feature_type", default= None,
                    help="Keep features from reference whose are tagged with this feature type")

parser.add_argument("-n", "--feature_name", dest="feature_name", default= None,
                    help="Use this feature id that is present in attributes/annotation field of reference")

opts = parser.parse_args()
  
options = vars(opts)
regions = Coord_Parser.load(options)
Genomic_Feature.add_reference(
	Reference_parser.load(
		options["reference_file"], 
		feature_type= options["feature_type"]
	)
)
gene_features = regions.get_features(attr_type= options["feature_name"])

with open(options["output_file"], 'w') as f:
	for id, feat_ids in gene_features.items():
		for ft_id in feat_ids:
			f.write(f"{id}\t{ft_id}\n")