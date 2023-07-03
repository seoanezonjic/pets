import re, sys

class File_Parser():
	@classmethod
	def load(cls, options): # To be overwrited in childclass
		return None

	@classmethod
	def get_header(cls, line, fields2extract):
		line = re.sub(r"#\s*", "", line) # correct comment like headers
		field_names = line.split("\t")
		cls.get_field_numbers2extract(field_names, fields2extract)
		field_numbers = fields2extract.values()
		return field_numbers

	@classmethod 
	def read_records(cls, options, fields2extract, field_numbers): # To be overwrited in childclass
		records = []

	@classmethod
	def get_records(cls, valid_fields, options):
		fields2extract = cls.get_fields2extract(valid_fields, options)
		field_numbers = fields2extract.values()
		records = cls.read_records(options, fields2extract, field_numbers)
		return fields2extract, records

	@classmethod
	def get_fields2extract(cls, valid_fields, options):
		fields2extract = {}
		for field in valid_fields:
			col = options.get(field)
			if col:
				if not options["header"]: col = int(col) 
				fields2extract[field] = col
		return fields2extract

	@classmethod
	def get_field_numbers2extract(cls, field_names, fields2extract):
		for field, name in fields2extract.items():
			fields2extract[field] = field_names.index(name)
