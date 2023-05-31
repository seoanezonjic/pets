import sys

class Genomic_feature:

	#If any method use gen_fet as name is a Genomic_Feature object
	def __init__(self, feat_list, annotations = None): # [[chr1, start1, stop1],[chr1, start1, stop1]]
		self.regions = {}
		self.reg_by_to = {}
		self.reg_id = -1
		self.load_features(feat_list)
		if annotations != None: self.load_annotations(annotations) 

	def load_features(self, feat_list):
		feat_len = len(feat_list[0])
		for ft_list in feat_list:
			self.reg_id +=1
			if feat_len == 4:
				chrm, start, stop, to = ft_list
				r_id = to
			elif feat_len == 3:
				chrm, start, stop = ft_list
				r_id = self.reg_id
			region = {'chrm': chrm, 'start': start, 'stop': stop, 'to': r_id }
			self.reg_by_to[r_id] = region
			self.add_record(self.regions, chrm, region)


	def load_annotations(self, annotations):
		for r in self.each(): 
			chrm, reg = r
			annot = annotations.get(reg['to'])
			if annot != None: reg['attrs'] = annot 

	def each(self):
		for chrm, regs in self.regions.items():
			for reg in regs: yield(chrm, reg)

	def each_chr(self):
		for chrm, regs in self.regions.items():
			yield(chrm, regs)

	def len(self):
		total_regions = 0
		for chrm, regs in self.each_chr(): total_regions += len(regs)
		return total_regions

	# Private
	def add_record(self, dic, key, record, uniq = False):
		query = dic.get(key)
		if query == None:
			dic[key] = [record]
		elif not uniq: # We not take care by repeated entries
			query.append(record)
		elif not record in query: # We want uniq entries
			query.append(record)