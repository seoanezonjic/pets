import re, sys
from pets import Cohort

class Cohort_Parser():
    
    @classmethod
    def load(cls, options):
        fields2extract = Cohort_Parser.get_fields2extract(options)
        field_numbers = fields2extract.values()
        records = Cohort_Parser.read_records(options, fields2extract, field_numbers)
        options["extracted_fields"] = fields2extract.keys()
        cohort, rejected_terms, rejected_recs = Cohort_Parser.create_cohort(records, options)
        return cohort, rejected_terms, rejected_recs

    @classmethod
    def read_records(cls, options, fields2extract, field_numbers):
        records = {}
        count = 0
        with open(options["input_file"]) as f:
            for line in f:
                line = line.strip()
                if options["header"] and count == 0:
                    line = re.sub(r"#\s*", "", line) # correct comment like	headers
                    field_names = line.split("\t")
                    Cohort_Parser.get_field_numbers2extract(field_names, fields2extract)
                    field_numbers = fields2extract.values()
                else:
                    fields = line.split("\t")
                    record = [fields[n] for n in field_numbers]
                    if fields2extract.get("id_col") == None:
                        id = f"rec_{count}" #generate ids
                    else:
                        id = record.pop(0)

                    if len(record) > 0 and record[0] != None:
                        record[0] = record[0].split(options["separator"])
                    else:
                        record[0] = []

                    if options.get("start_col"): record[2] = int(record[2])
                    if options.get("end_col"): record[3] = int(record[3]) 
                    query = records.get(id)
                    if query == None:
                        records[id] = [record]
                    else:
                        query.append(record)
                count +=1
        return records

    @classmethod
    def get_fields2extract(cls, options):
        fields2extract = {}
        for field in ["id_col", "ont_col", "chromosome_col", "start_col", "end_col", "sex_col"]:
            col = options.get(field)
            if col:
                if not options.get("header"): col = int(col)
                fields2extract[field] = col
        return fields2extract

    @classmethod
    def get_field_numbers2extract(cls, field_names, fields2extract):
        for field, name in fields2extract.items():
            fields2extract[field] = field_names.index(name)

    @classmethod
    def create_cohort(cls, records, options):
        ont = Cohort.get_ontology(Cohort.act_ont)
        rejected_terms = []
        rejected_recs = []
        cohort = Cohort()
        for id, record in records.items():
            rec = record[0]
            terms = rec[0]
            if options.get("names"): # Translate hpo names 2 codes
                init_term_number = len(terms)
                terms, rec_rejected_terms = ont.translate_names(terms)
                if rec_rejected_terms and len(rec_rejected_terms) > 0:
                    sys.stderr.write(f"WARNING: record {id} has the unknown term NAMES '{','.join(rec_rejected_terms)}'. Terms removed.\n")
                    rejected_terms.extend(rec_rejected_terms)

                if (not terms or len(terms) == 0) and init_term_number > 0:
                    rejected_recs.append(id)
                    continue

            if len(rec) > 1: # there is genomic region attributes
                variants = [v[1:4] for v in record]
            else:
                variants = [] # Not exists genomic region attributes so we create a empty array

            other_attr = {}
            if "sex_col" in options["extracted_fields"]: # Check for additional attributes. -1 is applied to ignore :id in extracted fields
                other_attr["sex"] = record[0][options["extracted_fields"].index("sex_col") -1]

            cohort.add_record([id, terms, Cohort_Parser.check_variants(variants)], other_attr)

        return cohort, list(set(rejected_terms)), rejected_recs

    @classmethod
    def check_variants(cls, vars):
        checked_vars = []
        for var in vars: #[chr, start, stop]
            if var[0] == '-': # the chr must be defined
                sys.stderr.write(f"WARNING: variant {(',').join(var)} has been removed") 
            else:
                checked_vars.append(var)
        return vars