import re, sys
from pets.parsers.file_parser import File_Parser
from pets.genomic_features import Genomic_Feature

class Coord_Parser(File_Parser):

	@classmethod
	def load(cls, options):
		valid_fields = ["id_col", "chromosome_col", "start_col", "end_col"]
		fields2extract, records = cls.get_records(valid_fields, options)
		genomic_features = Genomic_Feature(records)
		return genomic_features

	@classmethod
	def read_records(cls, options, fields2extract, field_numbers):
		records = []
		count = 0
		with open(options["input_file"]) as f:
			for line in f:
				line = line.strip()
				if options["header"] and count == 0:
					field_numbers = cls.get_header(line, fields2extract)
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