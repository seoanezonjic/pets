#! /usr/bin/env python
import argparse, os, re, sys
import pets
from pets import Genomic_Feature
from pets import Reference_parser

ROOT_PATH=os.path.dirname(__file__)

##################################
## METHODS
##################################

def get_data(options):
	fields2extract = get_fields2extract(options)
	field_numbers = fields2extract.values()
	records = read_records(options, fields2extract, field_numbers)
	return records


def read_records(options, fields2extract, field_numbers): # Modified from cohort_parset
	records = []
	count = 0
	with open(options["input_file"]) as f:
		for line in f:
			line = line.strip()
			if options["header"] and count == 0:
				line = re.sub(r"#\s*", "", line) # correct comment like	headers
				field_names = line.split("\t")
				get_field_numbers2extract(field_names, fields2extract)
				field_numbers = fields2extract.values()
			else:
				fields = line.split("\t")
				record = [fields[n] for n in field_numbers]
				if fields2extract.get("id_col") == None:
					id = f"rec_{count}" #generate ids
				else:
					id = record.pop(0)
				record[1] = int(record[1]) 
				record[2] = int(record[2])
				record.append(id)
				records.append(record)
			count +=1
	return records

def get_fields2extract(options):
	fields2extract = {}
	for field in ["id_col", "chromosome_col", "start_col", "end_col"]:
		col = options.get(field)
		if col:
			if not options["header"]: col = int(col) 
			fields2extract[field] = col
	return fields2extract

def get_field_numbers2extract(field_names, fields2extract):
	for field, name in fields2extract.items():
		fields2extract[field] = field_names.index(name)

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

parser.add_argument("-h", "--help", help="Show this message", action="store_true")

opts = parser.parse_args()

if opts.help:
    parser.print_help()
    sys.exit()
    
options = vars(opts)

regions = Genomic_Feature(get_data(options))
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