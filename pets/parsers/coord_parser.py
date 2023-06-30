import re, sys
from pets.genomic_features import Genomic_Feature

class Coord_Parser():

	@classmethod
	def load(cls, options):
		fields2extract = cls.get_fields2extract(options)
		field_numbers = fields2extract.values()
		records = cls.read_records(options, fields2extract, field_numbers)
		genomic_features = Genomic_Feature(records)
		return genomic_features

	@classmethod
	def read_records(cls, options, fields2extract, field_numbers): # Modified from cohort_parset
		records = []
		count = 0
		with open(options["input_file"]) as f:
			for line in f:
				line = line.strip()
				if options["header"] and count == 0:
					line = re.sub(r"#\s*", "", line) # correct comment like	headers
					field_names = line.split("\t")
					cls.get_field_numbers2extract(field_names, fields2extract)
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

	@classmethod
	def get_fields2extract(cls, options):
		fields2extract = {}
		for field in ["id_col", "chromosome_col", "start_col", "end_col"]:
			col = options.get(field)
			if col:
				if not options["header"]: col = int(col) 
				fields2extract[field] = col
		return fields2extract

	@classmethod
	def get_field_numbers2extract(cls, field_names, fields2extract):
		for field, name in fields2extract.items():
			fields2extract[field] = field_names.index(name)