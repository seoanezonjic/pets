import sys
from collections import defaultdict
import copy

class Genomic_Feature:
    ref = None

    @classmethod
    def array2genomic_feature(cls, arr, func, **kwargs):
        return cls([func(r) for r in arr], **kwargs)

    @classmethod
    def hash2genomic_feature(cls, h, func, **kwargs):
        vars = []
        for k, v in h.items():
            vars.append(func(k, v))
        return cls(vars, **kwargs)

    @classmethod
    def add_reference(cls, genomic_regions):
        cls.ref = genomic_regions

    #If any method use gen_fet as name is a Genomic_Feature object
    def __init__(self, feat_list, annotations = None): # [[chr1, start1, stop1],[chr1, start1, stop1]]
        self.regions = {}
        self.reg_by_to = {}
        self.reg_id = -1
        self.load_features(feat_list)
        if annotations != None: self.load_annotations(annotations)
    
    def __eq__(self, other):
        result = False
        if isinstance(other, Genomic_Feature) and self.regions == other.regions and self.reg_by_to == other.reg_by_to:
            result = True
        return result

    def load_features(self, feat_list):
        if not feat_list == None and len(feat_list) > 0:            
            for ft_list in feat_list:
                feat_len = len(ft_list)
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
        for chrm, reg in self.each(): 
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
    
    def get_chr(self):
        return self.regions.keys()
    
    def get_chr_regs(self, chrm):
        return self.regions.get(chrm)
    
    def region_by_to(self, to):
        return self.reg_by_to.get(to)
    
    def get_sizes(self):
        sizes = []
        for chr, region in self.each():
            size = region["stop"] - region["start"] + 1
            sizes.append(size)
        return sizes
  
    def get_features(self, attr_type= None):
        features = self.match(Genomic_Feature.ref) # [{self_id1: [ref_id1, ref_id5], self_id2: [ref_id8]}]
        def _get_attr(fi):
            attrs = Genomic_Feature.ref.region_by_to(fi).get("attrs") #Can be None or value
            if attrs != None: attrs = attrs.get(attr_type) #Can be None or value
            return attrs #Can be None or value
        
        if attr_type:
            for reg_id, feat_ids in features.items():
                new_feat_ids = list(map(_get_attr, feat_ids))            
                features[reg_id] = self.uniq_list([item for item in new_feat_ids if item != None])             
        return features


    def match(self, other_gen_feat):
        all_matches = {}
        for chr, regs in self.each_chr():
            other_regs = other_gen_feat.get_chr_regs(chr)
            if other_regs == None: continue
            for reg in regs: 
                local_matches = []
                start = reg["start"] 
                stop = reg["stop"] 
                for other_reg in other_regs:
                    if self.coor_overlap(start, stop, other_reg): local_matches.append(other_reg["to"])
                if local_matches: all_matches[reg["to"]] = local_matches
        return all_matches # [{reg_id1: [other_id1, other_id5], reg_id2: [other_id8]}]

    def get_summary_sizes(self):
        sizes = defaultdict(lambda: 0)
        for chrm, region in self.each():
            size = region["stop"] - region["start"] + 1
            sizes[size] += 1
        return sorted([list(sublist) for sublist in sizes.items()], key=lambda s: [s[1],s[0]], reverse=True)
 
    def merge(self, gen_fet, to = None): # 'to' the regions must be connected "to" given id
        for chrm, region in gen_fet.each():
            if to == None: 
                self.reg_id +=1
                region["to"] = self.reg_id 
            else:
                region["to"] = to
            self.add_record(self.regions, chrm, region)

    def get_reference_overlaps(self, genomic_ranges, reference):
		#Given a reference (that is the list of genomic windows for a given chromosome)
		#it returns the ids ("to") of the patient regions that fall in each of the windows,
        #ex: [[id1, id2, id3], [id4, id5, id6], [id7, id8, id9]]
        overlaps = []
        for start, stop in reference: 
            reg_ids = []
            for reg in genomic_ranges:
                overlap = self.coor_overlap(start, stop, reg)
                if overlap: reg_ids.append(reg["to"])
            overlaps.append(self.uniq_list(reg_ids))
        return overlaps
    

    def generate_cluster_regions(self, meth, tag, ids_per_reg = 1, obj = False):
        self.compute_windows(meth) # Get putative genome windows
        ids_by_cluster = {}
        annotated_full_ref = [] # All reference windows wit uniq id and chr tagged
        for chrm, regs in self.regions.items():
            reference = self.windows[chrm]
            overlaps = self.get_reference_overlaps(regs, reference) #[[id1, id2], [id4]]
            clust_numb = 0
            for i, ref in enumerate(reference):
                current_ids = overlaps[i] #[id1, id2]
                if len(current_ids) > ids_per_reg:
                    clust_id = f"{chrm}.{clust_numb}.{tag}.{len(current_ids)}"
                    clust_numb +=1
                    for curr_id in current_ids:
                        self.add_record(ids_by_cluster, curr_id, clust_id, True)
                    ref_copy = copy.deepcopy(ref)
                    ref_copy.extend([chrm, clust_id])
                    annotated_full_ref.append(ref_copy)
        if obj: annotated_full_ref = Genomic_Feature.array2genomic_feature(annotated_full_ref, lambda r: [r[2], r[0], r[1], r[3]])
        return ids_by_cluster, annotated_full_ref


    def compute_windows(self, meth):
        self.windows = {}
        for chrm, regs in self.regions.items():
            chr_windows = None
            if meth == "reg_overlap":
                chr_windows = self.compute_region_overlap_windows(regs)
            self.windows[chrm] = chr_windows

    # Private
    def add_record(self, dic, key, record, uniq = False):
        query = dic.get(key)
        if query == None:
            dic[key] = [record]
        elif not uniq: # We not take care by repeated entries
            query.append(record)
        elif not record in query: # We want uniq entries
            query.append(record)

    def compute_region_overlap_windows(self, genomic_ranges):
        reference = []
        single_nt = []
        for gr in genomic_ranges:
            start = gr["start"]
            stop = gr["stop"]
            if stop - start > 0:
                reference.append(start) # get start
                reference.append(stop) # get stop
            else: # Build a window of at least one nt for snv
                single_nt.append(start)

        reference = self.uniq_list(reference)    
        for snt in single_nt: # add start stop for snv
            reference.append(snt)
            reference.append(snt)

        reference.sort()
        #Define overlap ranges
        final_reference = []
        last_len = 1
        for i, coord in enumerate(reference): 
            if i+1 < len(reference):
                next_coord = reference[i + 1]
                current_len = next_coord - coord
                if last_len == 0: coord = coord + 1# Separate SNV window from others
                if current_len == 0 and last_len > 0 and len(final_reference) > 0:
                    final_reference[-1][1] -= 1 # Separate SNV window from others				
                final_reference.append([coord, next_coord])
                last_len = current_len
            
        return final_reference

    def coor_overlap(self, start, stop, reg):
        overlap = False
        reg_start = reg["start"]
        reg_stop = reg["stop"]
        if ((start <= reg_start and stop >= reg_stop) or
            (start > reg_start and stop < reg_stop) or
            (stop > reg_start and stop <= reg_stop) or
            (start >= reg_start and start < reg_stop)):
            overlap = True
        return overlap
    
    ## Accessory methods
    def flatten(self, nested_list):
        flat_list = []
        for item in nested_list:
            if isinstance(item, list):
                flat_list.extend(self.flatten(item))
            else:
                flat_list.append(item)
        return flat_list
    
    def uniq_list(self, list):
        used = set() #Doing this way to preserve order when doing uniq (because set does not preserve order)
        unique_list_with_initial_order = [item for item in list if item not in used and (used.add(item) or True)]
        return unique_list_with_initial_order