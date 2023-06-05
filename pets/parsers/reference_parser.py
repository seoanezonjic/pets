from pets import Genomic_Feature
class Reference_parser():
    
    @classmethod
    def load(cls, file_path, file_format: None, feature_type: None):
        if file_format == None: file_format = file_path.split('.', 2).last
        if file_format == 'gtf':
            regions, all_attrs = Reference_parser.parse_gtf(file_path, feature_type= feature_type)
        return Genomic_Feature(regions, annotations= all_attrs)

    @classmethod
    def parse_gtf(file_path, feature_type= None): # https://www.ensembl.org/info/website/upload/gff.html
        features = []
        all_attrs = {}
        with open(file_path) as f:
            for line in f:
                if line.startswith("#"): continue
            seqname, source, feature, start, stop, score, strand, frame, attribute = line.strip().split("\t")
            if feature_type == None or feature_type == feature:
                attrs = Reference_parser.process_attrs(attribute, ';', ' ')
                attrs['source'] = source
                attrs['feature'] = feature
                id = attrs['gene_id']
                features.append([seqname.replace('chr',''), int(start), int(stop), id])
                all_attrs[id] = attrs
        return features, all_attrs

    #private
    @classmethod
    def process_attrs(attributes, tuple_sep, field_sep):
        def format_tuple(attr_pair):
            tuple = attr_pair.strip().split(field_sep, maxsplit=2)
            tuple[-1] = tuple[-1].replace('"','')
            return tuple 
        
        return dict(list(
            map(lambda attr_pair: format_tuple(attr_pair), 
                attributes.split(tuple_sep))
                ))