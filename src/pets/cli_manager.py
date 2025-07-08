import argparse, inspect, sys
from pets.main_modules import * 

## TYPES
def tolist(string): 
  if string == "": return []
  return string.split(',')

def double_split(string, sep1=";", sep2=","):
    return [sublst.split(sep2) for sublst in string.strip().split(sep1)]

def loading_dic(string, sep1=";", sep2=","):
    return {key: value for key, value in double_split(string, sep1=sep1, sep2=sep2)}

##############################################
#Add parser common options
def add_parser_commom_options(parser):
    parser.add_argument("-c", "--chromosome_col", dest="chromosome_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the chromosome")

    parser.add_argument("-d", "--pat_id_col", dest="id_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the patient id")
    
    parser.add_argument("-s", "--start_col", dest="start_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate")

    parser.add_argument("-G", "--genome_assembly", dest="genome_assembly", default= "hg38",
                        help="Genome assembly version. Please choose between hg18, hg19 and hg38. Default hg38")

    #chr\tstart\tstop
    parser.add_argument("-H", "--header", dest="header", default= True, action="store_false",
                        help="Set if the file has a line header. Default true")

    parser.add_argument("-x", "--sex_col", dest="sex_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the patient sex")
    
    parser.add_argument("-S", "--hpo_separator", dest="separator", default='|',
                        help="Set which character must be used to split the HPO profile. Default '|'")

    parser.add_argument("-X", "--excluded_hpo", dest="excluded_hpo", default= None,
                        help="File with excluded HPO terms")
    
    parser.add_argument("--hard_check", dest="hard_check", default= True, action="store_false",
                        help="Set to disable hard check cleaning. Default true")    



###############################################
#CLI Scripts
###############################################
def monarch_entities(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    parser.add_argument('-i', '--input_file', dest='input_file', 
                        help='Input file with one entity id per line')
    parser.add_argument('-d', '--input_id', dest='input_id', 
                        help='Entity id to search with Monarch API')
    parser.add_argument('-o', '--output', dest='output', 
                        help='Output file')
    parser.add_argument('-r', '--relation', dest='relation', 
                        help='Describes which relation must be searched, It has been defines as "input_entity_type-desired_ouput_entity_type where type could be disease, phenotype, gene')
    parser.add_argument("-a","--add_entity2output", dest="add_entity2output", default= False, action="store_true",
                        help="Add queried entity as first column in output file")

    opts =  parser.parse_args(args)
    main_monarch_entities(opts)

def get_gen_features(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                        help="Input file path")

    parser.add_argument("-r", "--reference_file", dest="reference_file", default= None,
                        help="Reference file with genome annotation")

    parser.add_argument("-o", "--output_file", dest="output_file", default= None,
                        help="Output file with patient data")

    parser.add_argument("-t", "--feature_type", dest="feature_type", default= None,
                        help="Keep features from reference whose are tagged with this feature type")

    parser.add_argument("-n", "--feature_name", dest="feature_name", default= None,
                        help="Use this feature id that is present in attributes/annotation field of reference")
    parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")
    opts =  parser.parse_args(args)
    main_get_gen_features(opts)

def paco_translator(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-o", "--output_file", dest="output_file", default= None,
                        help="Output paco file with HPO names")

    parser.add_argument("-P", "--input_file", dest="input_file", default= None,
                        help="Input file with PACO extension")

    parser.add_argument("--n_phens", dest="n_phens", default= None, type=int,
                        help="Remove records with N or less phenotypes")

    parser.add_argument("-m","--save_mode", dest="save_mode", default= "default",
                        help="Set output data mode")

    parser.add_argument("-t","--translate", dest="translate", default= False, action="store_true",
                        help="Set to translate from hpo codes to names. By default, ther is not translation")

    parser.add_argument("-C", "--clean_PACO", dest="clean_PACO", default= False, action="store_true",
                        help="Clean PACO files")

    parser.add_argument("-r", "--removed_path", dest="removed_path", default= None,
                        help="Desired path to write removed profiles file")
    
    parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
                        help="Define if the input HPO are human readable names. Default false")
    
    parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                        help="Column name if header true or 0-based position of the column with the HPO terms")
    parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")
    opts =  parser.parse_args(args)
    main_paco_translator(opts)

def profiles2phenopacket(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                        help="Input file with patient data")
    parser.add_argument("-I", "--vcf_index", dest="vcf_index", default= None,
                        help="VCF file with patient id pointing to vcf path")
    parser.add_argument("-A", "--attr_index", dest="attr_index", default= None,
                        help="File with patient id and attribute in tabulated columns")
    parser.add_argument("--attr_name", dest="attr_name", default= None,
                        help="Attribute name to insert in phenopacket data")                        
    parser.add_argument("-o", "--output_file", dest="output_folder", default= None,
                        help="Output folder")
    parser.add_argument("--v2", dest="v2", default= False, action="store_true",
                        help="Use v2 phenopacket definition")
    parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
                        help="Define if the input HPO are human readable names. Default false")
    
    parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                        help="Column name if header true or 0-based position of the column with the HPO terms")
    parser.add_argument("--neg_hpo_term_col", dest="neg_ont_col", default= None,
                        help="Column name if header true or 0-based position of the column with the EXCLUDED HPO terms fora specific record")
    parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")
    opts =  parser.parse_args(args)
    main_profiles2phenopacket(opts)

def cohort_analyzer(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-a","--coverage_analysis", dest="coverage_analysis", default= True, action="store_false",
      help="Deactivate genome coverage analysis. Default true")
    parser.add_argument("-b", "--bin_size", dest="bin_size", default= 50000, type=int,
      help="Maximum number of bins to plot the coverage")
    parser.add_argument("-C", "--clusters2show", dest="clusters2show_detailed_phen_data", default= 3, type=int,
      help="How many patient clusters are show in detailed phenotype cluster data section. Default 3")
    parser.add_argument("-D","--detailed_clusters", dest="detailed_clusters", default= False, action="store_true",
      help="Show detailed cluster comparation using heatmaps. Default false")
    parser.add_argument("--detailed_cluster_yaxis", dest="detailed_cluster_yaxis", default= "cluster",
      help="What to show in heatmap clusters Y axis. 'cluster' for common phenotype in each cluster, 'cohort' for common phenotype in the whole cohort, 'cohort_sort' for same sorting for all clusters. Default: cluster")
    parser.add_argument("-M", "--minClusterProportion", dest="minClusterProportion", default= 0.01, type=float,
      help="Minimum percentage of patients per cluster")
    parser.add_argument("-f", "--patients_filter", dest="patients_filter", default= 2, type=int,
      help="Minimum number of patients sharing SORs. Default 2")
    parser.add_argument("-g", "--clusters2graph", dest="clusters2graph", default= 30, type=int,
      help="How may patient clusters are plotted in cluster plots. Default 30")
    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
      help="Input file with patient data")
    parser.add_argument("-m", "--similarity_methods", dest="clustering_methods", default=['lin'], type=tolist,
      help="Similarity methods to use in clustering step")
    parser.add_argument("-o", "--output_file", dest="output_file", default= None,
      help="Output file with patient data")
    parser.add_argument("-P", "--hpo_file", dest="hpo_file", default= None,
      help="Input HPO file for extracting HPO codes")
    parser.add_argument("-r", "--root_node", dest="root_node", default= "HP:0000118",
      help="Root node from which the ontology will be represented")
    parser.add_argument("-t", "--ic_stats", dest="ic_stats", default= "freq",
      help="'freq' to compute IC based en hpo frequency in the input cohort. 'freq_internal' to use precomputed internal IC stats. 'onto' to compute ic based on hpo ontology structure.. Default freq")
    parser.add_argument("-T", "--threads", dest="threads", default= 1, type=int,
      help="Number of threads to be used in calculations. Default 1")
    parser.add_argument("--reference_profiles", dest="reference_profiles", default= None,
      help="Path to file tabulated file with first column as id profile and second column with ontology terms separated by separator.")
    parser.add_argument("--sim_thr", dest="sim_thr", default=None, type=float,
      help="Keep pairs with similarity value >= FLOAT.")	
    parser.add_argument("--cl_size_factor", dest="cl_size_factor", default=1.0, type=float,
      help="When using dinamyc clustering weigths the contribution of the cluster size in tree cut. For smaller clusters use values > 1 for greater clusters use values < 1.")
    parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
      help="Define if the input HPO are human readable names. Default false")
    parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
      help="Column name if header true or 0-based position of the column with the HPO terms")
    parser.add_argument("-e", "--end_col", dest="end_col", default= None,
      help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")
    parser.add_argument("--ontoplot_mode", dest="ontoplot_mode", default='static',
      help="Use to set the mode for the ontoplot. Options: 'static','dynamic' or 'canvas'. Default is 'static'.")
    opts =  parser.parse_args(args)
    main_cohort_analyzer(opts)

def evidence_profiler(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-p", "--profiles_file", dest="profiles_file", default= None,
                        help="Path to profiles file. One profile per line and HP terms must be comma separated")

    parser.add_argument("-e", "--evidence_file", dest="evidence_file", default= None,
                        help="Path to evidence file. The file must have a column with the evidence id and a second column with the profile with HP terms comma separated")

    parser.add_argument("-g", "--genomic_coordinates_file", dest="coordinates_file", default= None,
                        help="Path to file with genomic coordinates for each genomic element in evidence file. One genomic element per line with format id, chr and start position.")

    parser.add_argument("-E", "--evidence_folder", dest="evidences", default= None,
                        help="Path to evidence folder.")

    parser.add_argument("-o", "--output_folder", dest="output_folder", default= "evidence_reports",
                        help="Folder to save reports from profiles")

    parser.add_argument("-V", "--variant_data", dest="variant_data", default= None,
                        help="Folder to tables of patient variants")

    parser.add_argument("-P", "--pathogenic_scores", dest="pathogenic_scores", default= None, # TODO: Generalize to a folder with a table per patient
                        help="File with genome features an their pathogenic scores")
    opts =  parser.parse_args(args)
    main_evidence_profiler(opts)

def diseasome_generator(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)


    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                    help="Input file with the ontology term and genes")
    parser.add_argument("-C", "--disorder_class", dest="disorder_class", default= None,
                    help="Input file with the ontology terms and its respective text")
    parser.add_argument("-O", "--ontology", dest="ontology", default= None, 
                    help="Path to the ontology file")
    parser.add_argument("-g", "--generate", dest = "generate", default=False, action="store_true",
        help = "To generate a new diseasome")
    parser.add_argument("-A", "--analyze", dest= "analyze", default=False, action="store_true",
        help = "Analyze the diseasome given some descriptive stats")
    parser.add_argument("-D", "--diseasome", dest="diseasome", default=None,
        help = "a path for a file conatining a diseasome following the next TABULATED format: MONDO_disease_term, disorder_class")
    parser.add_argument("-o", "--output_file", dest="output_file", default= None, 
                    help="Path to the output file to write results")

    opts = parser.parse_args(args)
    main_diseasome_generator(opts)

def filter_omim(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                    help="Input file with the ontology term and genes")
    parser.add_argument("-o", "--output_file", dest="output_file", default= None, 
                    help="Path to the output file to write results")

    opts = parser.parse_args(args)
    main_filter_omim(opts)

def collapse_terms(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
                    help="Input file with the ontology terms")
    parser.add_argument("-n", "--terms2text", dest="terms2text", default= None,
                    help="Input file with the ontology terms and its respective text")
    parser.add_argument("-O", "--ontology", dest="ontology", default= None, 
                    help="Path to the ontology file")
    parser.add_argument("-o", "--output_file", dest="output_file", default= None, 
                    help="Path to the output file to write results")
    parser.add_argument("--sim_algorithm", dest="sim_algorithm", default= 'white',
                    help="Similarity algorithm")
    parser.add_argument("-r", "--remove_chars", dest="rm_char", default="", 
                    help="Chars to be excluded from comparissons.")
    parser.add_argument("-t", "--threshold", dest="threshold", default=0.70, type=float,
                    help="Threshold to consider a pair of terms similar")
    parser.add_argument("-u","--uniq_parent", dest="uniq_parent", default= False, action = "store_true",
                    help="Just add the uniq most representative parent")
    parser.add_argument("-l","--just_leaves", dest="just_leaves", default=False, action = "store_true",
        help= "When just leaves want to be added on the collapse terms")

    opts = parser.parse_args(args)

    main_collapse_terms(opts)

def report_prioritizer(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("--prioritizers", dest="prioritizers", default= None, type=lambda x: loading_dic(x, sep1=";", sep2=","),
                    help="Format prioritizer:path_to_prioritizer_file")
    parser.add_argument("--integrated_report", dest="integrated_report", default= False, action="store_true",
                    help="Select if integration is needed")
    parser.add_argument("--comparing_report", dest="report", default= None, 
                    help="Path to the comparing report file")
    parser.add_argument("--read_tmp",dest="read_tmp",default=False, action="store_true",
                        help="Read processed file for report analysis")
    parser.add_argument("--write_tmp",dest="write_tmp",default=None, type=str,
                    help="Write processed file for report analysis. This flag cancel the posterior report analysis")
    parser.add_argument("--benchmark_type",dest="benchmark_type",default="gene", type=str,
                    help="Choose type of benchmark. Choose between 'gene' and 'variant', or both. Options: 'variant', 'gene'")
    parser.add_argument("-o", "--output_file", dest="output_file", default= "report_prioritizer", 
                    help="Path to the output file to write results")
    opts = parser.parse_args(args)
    main_report_prioritizer(opts)

def meta_prioritizer(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)
    parser.add_argument("--prioritizers", dest="prioritizers", default= None, type=lambda x: loading_dic(x, sep1=";", sep2=","),
                    help="Format prioritizer:path_to_prioritizer_file")
    parser.add_argument("--labels", dest="labels", help="""path to file with format in three columns: 
                        patient\tcandidateId""") 
    parser.add_argument("--mode", dest="mode", default="both", type=str, help="To select if you want to test and or train a model, the options are 'train', 'predict' or 'both'. Default both")
    parser.add_argument("--training_options", dest="training_options", default=None, type=lambda x: loading_dic(x, sep1=";", sep2=","),
                    help="Format option:value for the training model")
    parser.add_argument("--read_tmp",dest="read_tmp",default=False, action="store_true",
                        help="Read processed file for report analysis")
    parser.add_argument("--model_type", dest="model_type", default="xgboost", type=str,
                        help="Type of model to use. Options: 'heuristic', 'xgboost', 'logistic_regression' . Default xgboost")
    parser.add_argument("--model_path", dest="model_path", default=None, type=str,
                    help="Path to the model file. If not provided, a new model will be trained.")
    parser.add_argument("--write_tmp",dest="write_tmp",default=None, type=str,
                    help="Write processed file for report analysis. This flag cancel the posterior report analysis")
    parser.add_argument("--benchmark_type",dest="benchmark_type",default="gene", type=str,
                    help="Choose type of benchmark. Choose between 'gene' and 'variant', or both. Options: 'variant', 'gene'")
    parser.add_argument("-o", "--output_file", dest="output_file", default= "report_prioritizer", 
                    help="Path to the output file to write results")
    parser.add_argument("--dropna", dest="dropna", default= False, action="store_true",
                    help="Select this flag if model do not tolerate NaNs values")
    opts = parser.parse_args(args)
    main_meta_prioritizer(opts)

def phenPatMaster(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    parser.add_argument("-i", "--input_folder", dest="input_folder", default= None,
                    help="Input phenopacket folder.")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default= None,
                    help="Output phenopacket folder.")
    parser.add_argument("--overwrite_id", dest="overwrite_id", default= False, action = 'store_true',
                    help="Generate automatic id for phenopackets. it will used in ALL id fields.")
    parser.add_argument("--overwrite_file_name", dest="overwrite_file_name", default= False, action = 'store_true',
                    help="When used, the output phenopacket uses phenopacket id as file name.")
    parser.add_argument("--clean_phen", dest="clean_phen", default= False, action = 'store_true',
                    help="When used, check HP terms, update obsolete and remove parental relations.")
    parser.add_argument("--output_file_index", dest="output_file_index", default= None,
                    help="Generate a index with phenotypical and genomical features described in phenopackets.")

    opts = parser.parse_args(args)

    main_phenPatMaster(opts)


def vcf2effects(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    parser.add_argument("-i", "--input_vcf", dest="input_vcf",
      help="Vcf file")
    parser.add_argument("-o", "--output", dest="output",
      help="output with variant effects")
    parser.add_argument("-g", "--genome", dest="genome",
      help="Genome name to use as reference")
    opts = parser.parse_args(args)

    main_vcf2effects(opts)