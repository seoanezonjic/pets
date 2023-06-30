from pets import Genomic_Feature
import warnings
class Reference_parser():
    
    @classmethod
    def load(cls, file_path, file_format = None, feature_type = None):
        if file_format == None: file_format = file_path.split('.', maxsplit=2)[-1]
        if file_format == 'gtf':
            regions, all_attrs = cls.parse_gtf(file_path, feature_type=feature_type)
        return Genomic_Feature(regions, annotations= all_attrs)

    @classmethod
    def parse_gtf(cls, file_path, feature_type= None): # https://www.ensembl.org/info/website/upload/gff.html
        features = []
        all_attrs = {}
        with open(file_path) as f:
            for line in f:
                if line.startswith("#"): continue
                seqname, source, feature, start, stop, score, strand, frame, attribute = line.strip().split("\t")
                if feature_type == None or feature_type == feature:
                    attrs = cls.process_attrs(attribute, ';', ' ')
                    attrs['source'] = source
                    attrs['feature'] = feature
                    id = attrs['gene_id']
                    features.append([seqname.replace('chr',''), int(start), int(stop), id])
                    all_attrs[id] = attrs
        return features, all_attrs

    #private
    @classmethod
    def process_attrs(cls, attributes, tuple_sep, field_sep):
        attrs_dict = {}
        for attr_pair in attributes.split(tuple_sep):
            if len(attr_pair) == 0: continue
            tuple = attr_pair.strip().split(field_sep, maxsplit=2)
            tuple[-1] = tuple[-1].replace('"','')
            if len(tuple) == 2:
                attrs_dict[tuple[0]] = tuple[1]
            else:
                warnings.warn(f"Attribute {attr_pair} from {attributes} is not a tuple of length 2.")
        return attrs_dict