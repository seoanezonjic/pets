import argparse
import sys
import os
import re
import inspect
from importlib.resources import files
import urllib.parse
import requests

import numpy as np
from py_report_html import Py_report_html
import pets
from pets.cohort_analyser_methods import *
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort
from pets.io import load_hpo_ci_values, write_tabulated_data, load_profiles, load_variants, load_evidences, load_index, parse_morbid_omim, list2dic
from pets.genomic_features import Genomic_Feature
from pets.parsers.reference_parser import Reference_parser
from pets.parsers.coord_parser import Coord_Parser
from py_exp_calc.exp_calc import invert_hash, uniq
from py_semtools.ontology import Ontology
from py_semtools.sim_handler import similitude_network
from NetAnalyzer import NetAnalyzer
import pandas as pd
from py_cmdtabs.cmdtabs import CmdTabs

# https://setuptools.pypa.io/en/latest/userguide/datafiles.html
HPO_FILE = str(files('pets.external_data').joinpath('hp.json'))
MONDO_FILE = str(files('pets.external_data').joinpath('mondo.obo'))
IC_FILE = str(files('pets.external_data').joinpath('uniq_hpo_with_CI.txt'))
GENCODE = str(files('pets.external_data').joinpath('gencode.v43.basic.annotation.gtf.gz'))


## TYPES
def tolist(string): 
  if string == "": return []
  return string.split(',')

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

def get_sorted_profs(args=None):
    if args == None: args = sys.argv[1:]
    parser = argparse.ArgumentParser(description=f'Usage: {inspect.stack()[0][3]} [options]')
    add_parser_commom_options(parser)

    parser.add_argument("-o", "--output_file", dest="output_file", default= 'report.html',
                        help="Output paco file with HPO names")

    parser.add_argument("-P", "--input_file", dest="input_file", default= None,
                        help="Input file with PACO extension")

    parser.add_argument("-f", "--general_prof_freq", dest="term_freq", default= 0, type= int,
                        help="When reference profile is not given, a general ine is computed with all profiles. If a freq is defined (0-1), all terms with freq minor than limit are removed")

    parser.add_argument("-L", "--matrix_limits", dest="matrix_limits", default= [20, 40], type= lambda data: [int(i) for i in data.split(",")],
                        help="Number of rows and columns to show in heatmap defined as 'Nrows,Ncols'. Default 20,40")

    parser.add_argument("-r", "--ref_profile", dest="ref_prof", default= None, 
                        type = lambda file: [line.strip() for line in open(file).readlines()],
                        help="Path to reference profile. One term code per line")
    
    parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                        help="Column name if header true or 0-based position of the column with the HPO terms")
    parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")
    opts =  parser.parse_args(args)
    main_get_sorted_profs(opts)

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
                        
    parser.add_argument("-o", "--output_file", dest="output_folder", default= None,
                        help="Output folder")
    
    parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
                        help="Define if the input HPO are human readable names. Default false")
    
    parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                        help="Column name if header true or 0-based position of the column with the HPO terms")
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
    parser.add_argument("-M", "--minClusterProportion", dest="minClusterProportion", default= 0.01, type=float,
      help="Minimum percentage of patients per cluster")
    parser.add_argument("-f", "--patients_filter", dest="patients_filter", default= 2, type=int,
      help="Minimum number of patients sharing SORs. Default 2")
    parser.add_argument("-g", "--clusters2graph", dest="clusters2graph", default= 30, type=int,
      help="How may patient clusters are plotted in cluster plots. Default 30")
    parser.add_argument("-i", "--input_file", dest="input_file", default= None,
      help="Input file with patient data")
    parser.add_argument("-m", "--clustering_methods", dest="clustering_methods", default=['lin'], type=tolist,
      help="Clustering methods")
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
    parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
      help="Define if the input HPO are human readable names. Default false")
    parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
      help="Column name if header true or 0-based position of the column with the HPO terms")
    parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                        help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")
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

###########################################################
# Main functions
###########################################################

def main_collapse_terms(opts):
    options = vars(opts)
    ontology_file = options["ontology"] if options["ontology"] else MONDO_FILE
    ontology = Ontology(file= ontology_file, load_file = True) 

    input_file = CmdTabs.load_input_data(options["input_file"])
    terms_to_annot = None
    if len(input_file[0]) > 1:
        terms_to_annot =  list2dic(input_file)
        terms = list(set(terms_to_annot.keys()))
    else:
        terms = list(set([term[0] for term in input_file]))

    if options["terms2text"]:
        term_to_txt = load_index(options.get("terms2text"))
        #term_to_txt = {key: value[0] for key, value in term_to_txt.items()}
    else: 
        term_to_txt = {term: ontology.translate_id(term) for term in terms}
    txt_to_term = {value: key for key, value in term_to_txt.items()}

    terms2collapse = [term for term in terms if term in term_to_txt.keys()] # selecting just terms with txt
    terms2collapse = list(set(terms2collapse))
    terms2collapse = filter_out_non_leafs_nodes(list(set(terms)), ontology) if options["just_leaves"] else terms2collapse

    parent_to_childs_terms = get_parent_and_childs_nodes_dict(terms2collapse, ontology)
    parent_to_childs_txt = {parent: [term_to_txt[child] for child in childs] for parent, childs in parent_to_childs_terms.items()}
    similarities = get_txt_to_txt_similarities(parent_to_childs_txt, rm_char=options["rm_char"])
    terms_to_parents_collapsed = get_thresholded_childs_to_parents_dict(similarities, threshold = options.get("threshold"),
     ontology=ontology, rm_char = options["rm_char"], txt_to_term = txt_to_term, uniq_parent = options["uniq_parent"])

    with open(options["output_file"], "w") as f:
        # Wath out: This is implying that all term which are not leafs would be added!
        for term in terms:
            term_id = ""
            if terms_to_parents_collapsed.get(term) is not None:
                parents = terms_to_parents_collapsed[term]
            else:
                parents = [term]
            for parent in parents:
                if terms_to_annot:
                    annotations = "\t".join(terms_to_annot[term])
                    f.write(f"{term}\t{annotations}\t{parent}\n")
                else:
                    f.write(f"{term}\t{parent}\n")

def main_filter_omim(opts):
    options = vars(opts)
    parsed_morbid = parse_morbid_omim(options["input_file"])
    with open(options["output_file"], "w") as f:
        for omim_code, annots in parsed_morbid.items():
            genes = ",".join(uniq(annots[1]))
            omim_txt = annots[0]
            f.write("\t".join([omim_code, omim_txt, genes])+"\n")

def main_monarch_entities(opts):
    entities = []
    if opts.input_id != None:
        entities.append(opts.input_id)
    else:
        with open(opts.input_file) as f:
            for line in f: entities.append(line.rstrip())

    base_url='api-v3.monarchinitiative.org/v3/api'
    limit = 500
    with open(opts.output, 'w') as f:
        for entity in entities:
            data = get_monarch_data(entity, opts.relation, base_url, limit)
            for d in data:
                if opts.add_entity2output: d = [entity] + d
                f.write("\t".join(d) + "\n")    

def get_monarch_data(entity, relation, base_url, limit):
    retrieved_data = []
    if(relation == 'phenotype-disease'):
        api_url = f"{base_url}/entity/{entity}/biolink:DiseaseToPhenotypicFeatureAssociation"
    elif(relation == 'disease-gene'):
        api_url = f"{base_url}/entity/{entity}/biolink:CausalGeneToDiseaseAssociation"
    query = 'https://' + urllib.parse.quote(api_url)+ f'?limit={limit}'
    total = 1
    recovered = 0
    while total > recovered:
        if recovered == 0:
            offset = 0
        else:
            offset = recovered - 1
        final_query = query + f"&offset={offset}"
        response = requests.get(final_query)
        data = response.json()
        total = data['total']
        items = data['items']
        recovered += len(items)
        for it in items:
            entity = it['subject']
            if(relation == 'phenotype-disease'):
                original_entity = it['original_subject']
                fields = [entity, original_entity]
            elif(relation == 'disease-gene'):
                gene_name = it['subject_label']
                fields = [entity, gene_name]
            retrieved_data.append(fields)
    return retrieved_data

def main_diseasome_generator(opts):
    # Loading and parsing inputs
    options = vars(opts)
    ontology_file = options["ontology"] if options["ontology"] else MONDO_FILE
    ontology = Ontology(file= options["ontology"], load_file = True, extra_dicts=[['xref', {'select_regex': "OMIM:[0-9]*", 'store_tag': 'tag', 'multiterm': False}]])
    ontology.dicts['tag']["byTerm"] = {value: key for key, value in ontology.dicts['tag']["byValue"].items()}

    if options["disorder_class"]:
        disorder_class_file = options["disorder_class"]
    else:
        disorder_class_file = str(files('pets.external_data').joinpath('disorder_classes'))
        
    disorder_class = {}
    disorder_class = load_index(disorder_class_file, True)
    diseasome = None
    if options["diseasome"]:
        diseasome = load_index(options["diseasome"])
        diseasome = {ontology.dicts["tag"]["byValue"][omim]: group for omim, group in diseasome.items() if ontology.dicts["tag"]["byValue"].get(omim)}

    # Execution modes
    if options["generate"] and not diseasome:
        diseases =  CmdTabs.load_input_data(options["input_file"])
        disease_annotation = {}
        if len(diseases[0]) > 1: disease_annotation = list2dic(diseases)
        omim_diseases = list(set([disease[0] for disease in diseases])) 
        mondo_diseases = [ontology.dicts['tag']["byValue"][omim] for omim in omim_diseases if ontology.dicts['tag']["byValue"].get(omim)]

        disease2disclass = get_dis2dclass(mondo_diseases, disorder_class, ontology)
        dependency_map = get_dependency_map(set(disorder_class.keys()), ontology)
        diseasome = clean_dis2dclass(disease2disclass, ontology, dependency_map, disorder_class)
    
        with open(options["output_file"], "w") as f:
            for disease in mondo_diseases:
                if diseasome.get(disease):
                        # Pasamos a omim -> mondo -> annotations -> mondo parent
                        omim = ontology.dicts['tag']["byTerm"][disease]
                        if disease_annotation:
                            annotations = "\t".join(disease_annotation[omim])
                            f.write(f"{omim}\t{disease}\t{annotations}\t{diseasome[disease]}\n")
                        else:
                            f.write(f"{omim}\t{disease}\t{diseasome[disease]}\n")

    if options["analyze"]:
        tripartite = generate_tripartite_diseasome(list(diseasome.items()), ontology, ["disease","group","parentals"])
        tripartite.compute_autorelations = False
        projection = tripartite.get_association_values(("group", "parentals"), "disease", "counts")
        projection = pd.DataFrame(projection, columns=["disgroup","parental","counts"])
        select_top = lambda group: group.nlargest(20, 'counts')
        # Group by 'disgroup' and apply the function to select the top values within each group
        best_picks = projection.groupby("disgroup").apply(select_top)
        best_picks = best_picks.values.tolist()
        with open(options["output_file"]+"_analysis", "w") as f:
            for line in best_picks:
                line = [line[0],ontology.translate_id(line[1]),str(line[2])]
                f.write("\t".join(line)+"\n")

def main_get_gen_features(opts):
    options = vars(opts)
    regions = Coord_Parser.load(options)
    reference = options["reference_file"] if options.get("reference_file") else GENCODE
    Genomic_Feature.add_reference(
        Reference_parser.load(
            reference, 
            feature_type= options["feature_type"]
        )
    )
    gene_features = regions.get_features(attr_type= options["feature_name"])

    with open(options["output_file"], 'w') as f:
        for id, feat_ids in gene_features.items():
            for ft_id in feat_ids:
                f.write(f"{id}\t{ft_id}\n")


def main_get_sorted_profs(opts):
    options = vars(opts)
    hpo_file = os.environ.get('hpo_file') if os.environ.get('hpo_file') else HPO_FILE
    Cohort.load_ontology("hpo", hpo_file)
    Cohort.act_ont = "hpo"
    hpo = Cohort.get_ontology(Cohort.act_ont)
    patient_data, _, _ = Cohort_Parser.load(options)
    patient_data.check(hard=True)

    clean_profiles = patient_data.profiles

    if options.get("ref_prof"):
      ref_profile = hpo.clean_profile_hard(options["ref_prof"])
    else:
      ref_profile = patient_data.get_general_profile(options["term_freq"])

    hpo.load_profiles({"ref": ref_profile}, reset_stored= True)

    similarities = hpo.compare_profiles(external_profiles= clean_profiles, sim_type= "lin", bidirectional= False)

    candidate_sim_matrix, candidates, candidates_ids = patient_data.get_term2term_similarity_matrix(ref_profile, similarities["ref"], clean_profiles, hpo, options["matrix_limits"][0], options["matrix_limits"][-1])
    candidate_sim_matrix.insert(0, ['HP'] + candidates_ids)

    template = open(str(files('pets.templates').joinpath('similarity_matrix.txt'))).read()
    container = { "similarity_matrix": candidate_sim_matrix }
    report = Py_report_html(container, 'Similarity matrix')
    report.build(template)
    report.write(options["output_file"])

    with open(options["output_file"].replace('.html','') +'.txt', 'w') as f:
      for candidate, value in similarities["ref"].items():
        f.write("\t".join([str(candidate), str(value)])+"\n")



def main_paco_translator(opts):
    options = vars(opts)
    hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
    Cohort.load_ontology("hpo", hpo_file, options.get("excluded_hpo"))
    Cohort.act_ont = "hpo"

    patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(options)

    if options.get("clean_PACO"):
        removed_terms, removed_profiles = patient_data.check(hard=True)
        if options.get("removed_path") and removed_profiles != None and len(removed_profiles) > 0:
            rejected_file = os.path.basename(options["input_file"]).split(".")[0] +'_excluded_patients'
            file = os.path.join(options["removed_path"], rejected_file)
            with open(file, 'w') as f:
                for profile in removed_profiles:
                    f.write(profile+'\n')

    if options.get("n_phens"): rejected_patients_by_phen = patient_data.filter_by_term_number(options["n_phens"])
    patient_data.save(options["output_file"], options["save_mode"], options["translate"])


def main_profiles2phenopacket(opts):
    options = vars(opts)
    hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
    Cohort.load_ontology("hpo", hpo_file, options.get("excluded_hpo"))
    Cohort.act_ont = "hpo"

    patient_data, rejected_hpos_L, rejected_patients_L = Cohort_Parser.load(options)
    rejected_hpos_C, rejected_patients_C = patient_data.check(hard=True)
    patient_data.link2ont(Cohort.act_ont)

    vcf_index = None
    if options.get("vcf_index"): vcf_index = load_index(options["vcf_index"])
    patient_data.export_phenopackets(options["output_folder"], options["genome_assembly"], vcf_index= vcf_index)


def main_cohort_analyzer(options):
    opts = vars(options)
    if opts['genome_assembly'] == 'hg19' or opts['genome_assembly'] == 'hg37':
      CHR_SIZE = str(files('pets.external_data').joinpath('chromosome_sizes_hg19.txt'))
    elif opts['genome_assembly'] == 'hg38':
      CHR_SIZE = str(files('pets.external_data').joinpath('chromosome_sizes_hg38.txt'))
    elif opts['genome_assembly'] == 'hg18':
      CHR_SIZE = str(files('pets.external_data').joinpath('chromosome_sizes_hg18.txt'))
    else:
      raise Exception('Wrong human genome assembly. Please choose between hg19, hg18 or hg38.')
    chr_sizes = dict(map(lambda sublist: [sublist[0], int(sublist[1])], [line.strip().split("\t") for line in open(CHR_SIZE).readlines()]))

    output_folder = os.path.dirname(opts['output_file'])
    detailed_profile_evaluation_file = os.path.join(output_folder, 'detailed_hpo_profile_evaluation.csv')
    rejected_file = os.path.join(output_folder, 'rejected_records.txt')
    temp_folder = os.path.join(output_folder, 'temp')
    hpo_frequency_file = os.path.join(temp_folder, 'hpo_cohort_frequency.txt')

    if not os.path.exists(temp_folder): os.mkdir(temp_folder)

    hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
    Cohort.load_ontology("hpo", hpo_file, opts.get("excluded_hpo"))
    Cohort.act_ont = "hpo"

    opts['check'] = True
    patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(opts)
    with open(rejected_file, 'w') as f: f.write("\n".join(rejected_patients))

    patient_data.link2ont(Cohort.act_ont) # TODO: check if method load should call to this and use the semtools checking methods (take care to only remove invalid terms)

    profile_sizes, parental_hpos_per_profile = patient_data.get_profile_redundancy()
    patient_data.check(hard=True)
    hpo_stats = patient_data.get_profiles_terms_frequency() # hpo NAME, freq
    for stat in hpo_stats: stat[1] = stat[1]*100
    with open(hpo_frequency_file, 'w') as f:
      for hpo_code, freq in patient_data.get_profiles_terms_frequency(translate= False): # hpo CODE, freq
        f.write(f"{hpo_code}\t{freq}\n")

    suggested_childs, fraction_terms_specific_childs = patient_data.compute_term_list_and_childs(file = detailed_profile_evaluation_file)

    ontology_levels, distribution_percentage = patient_data.get_profile_ontology_distribution_tables()
    onto_ic, freq_ic, onto_ic_profile, freq_ic_profile = patient_data.get_ic_analysis()

    if opts['ic_stats'] == 'freq_internal': # TODO: Make semtools to load ci external values
      ic_file = os.environ['ic_file'] if os.environ.get('ic_file') else IC_FILE
      freq_ic = load_hpo_ci_values(ic_file)
      phenotype_ic = freq_ic
      freq_ic_profile = {}
      for pat_id, phenotypes in patient_data.each_profile():
        freq_ic_profile[pat_id] = patient_data.get_profile_ic(phenotypes, phenotype_ic)
    elif opts['ic_stats'] == 'freq':
      phenotype_ic = freq_ic
    elif opts['ic_stats'] == 'onto':
      phenotype_ic = onto_ic

    all_ics, prof_lengths, clust_by_chr, top_clust_phen, multi_chr_clusters = patient_data.process_dummy_clustered_patients(opts, phenotype_ic, temp_folder = temp_folder)

    summary_stats = get_summary_stats(patient_data, rejected_patients, hpo_stats, fraction_terms_specific_childs, rejected_hpos)

    all_cnvs_length = []
    all_sor_length = []
    if not opts.get('chromosome_col') == None:
      summary_stats.append(['Number of clusters with mutations accross > 1 chromosomes', multi_chr_clusters])
      
      patient_data.index_vars()
      all_cnvs_length = patient_data.get_vars_sizes(False)
      summary_stats.append(['Average variant size', round(np.mean(all_cnvs_length), 4)])
      #----------------------------------
      # Prepare data to plot coverage
      #----------------------------------
      if opts.get('coverage_analysis'):
        patients_by_cluster, sors = patient_data.generate_cluster_regions('reg_overlap', 'A', 0)

        ###1. Process CNVs
        raw_coverage, n_cnv, nt, pats_per_region = calculate_coverage(sors)
        summary_stats.extend(
          [['Nucleotides affected by mutations', nt],
          ['Number of genome windows', n_cnv],
          ['Mean patients per genome window', round(pats_per_region, 4)]])
        coverage_to_plot = get_final_coverage(raw_coverage, opts['bin_size'])
        ###2. Process SORs
        raw_sor_coverage, n_sor, nt, pats_per_region = calculate_coverage(sors, opts['patients_filter'] - 1)
        summary_stats.extend(
          [[f"Number of genome window shared by >= {opts['patients_filter']} patients", n_sor],
          ["Number of patients with at least 1 SOR", len(patients_by_cluster)],
          ['Nucleotides affected by mutations', nt]])
        sor_coverage_to_plot = get_final_coverage(raw_sor_coverage, opts['bin_size'])

        all_sor_length = get_sor_length_distribution(raw_sor_coverage)  

    dummy_cluster_chr_data = []
    if not opts.get('chromosome_col') == None:
      dummy_cluster_chr_data = get_cluster_chromosome_data(clust_by_chr, opts['clusters2graph'])

    #----------------------------------
    # CLUSTER COHORT ANALYZER REPORT
    #----------------------------------
    reference_profiles = None
    if opts.get('reference_profiles') != None: reference_profiles = load_profiles(opts['reference_profiles'], Cohort.get_ontology('hpo'))
    template = str(files('pets.templates').joinpath('cluster_report.txt'))
    clustering_data = get_semantic_similarity_clustering(opts, patient_data, reference_profiles, temp_folder, template)

    #----------------------------------
    # GENERAL COHORT ANALYZER REPORT
    #----------------------------------
    new_cluster_phenotypes = get_top_dummy_clusters_stats(top_clust_phen)

    container = {
      'temp_folder' : temp_folder,
      # 'top_clust_phen' : len(top_clust_phen),
      'summary_stats' : summary_stats,
      'clustering_methods' : opts['clustering_methods'],
      'hpo_stats' : hpo_stats,
      'all_cnvs_length' : [ [l] for l in all_cnvs_length ],
      'all_sor_length' : [ [l] for l in all_sor_length ],
      'new_cluster_phenotypes' : len(new_cluster_phenotypes),
      'ontology_levels' : ontology_levels,
      'distribution_percentage' : distribution_percentage,
      'hpo_ic_data': [ list(p) for p in zip(list(onto_ic.values()),list(freq_ic.values())) ],
      'hpo_ic_data_profiles': [ list(p) for p in zip(list(onto_ic_profile.values()), list(freq_ic_profile.values())) ],
      'parents_per_term': [ list(p) for p in zip(profile_sizes, parental_hpos_per_profile) ],
      'dummy_cluster_chr_data' : dummy_cluster_chr_data,
      'dummy_ic_data' : format_cluster_ic_data(all_ics, prof_lengths, opts['clusters2graph']),
      'chr_sizes': chr_sizes,
      #'term_freq_table': dict(Cohort.ont['hpo'].get_profiles_terms_frequency(translate = False)), #TODO: check, it is giving all the frequencies equal
      'term_freq_table': {Cohort.ont['hpo'].translate_name(hpo): value/100 for hpo, value in (dict(hpo_stats)).items()}, 
      'ontology': Cohort.ont['hpo']
    }

    clust_info = []
    for clusterID, info in new_cluster_phenotypes.items():
        phens = ', '.join(info[1])
        freqs = ', '.join([ str(round(a,4)) for a in info[2]])
        container[f"clust_{clusterID}"] = [[info[0], phens, freqs]]

    for meth, data in clustering_data.items(): 
      for item, obj in data.items(): container[f"{meth}_{item}"] = obj

    if opts.get('chromosome_col') != None and opts['coverage_analysis']:
      coverage_to_plot.insert(0, ['Chr', 'Pos', 'Count'])
      container['cnv_coverage'] = coverage_to_plot # chr, start bin, count
      sor_coverage_to_plot.insert(0, ['Chr', 'Pos', 'Count'])
      container['sor_coverage'] = sor_coverage_to_plot

    report = Py_report_html(container)
    report.build(open(str(files('pets.templates').joinpath('cohort_report.txt'))).read())
    report.write(opts["output_file"] + '.html')


def main_evidence_profiler(opts):
    options = vars(opts)
    hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE

    hpo = Ontology(file= hpo_file, load_file= True)
    cohort = Cohort() # needed until get_similarity

    profiles = load_profiles(options["profiles_file"], hpo)
    profile_variants = {} if options.get("variant_data") == None else load_variants(options["variant_data"])
    evidences, genomic_coordinates = load_evidences(options["evidences"], hpo)
    pathogenic_scores = {} if options.get("pathogenic_scores") == None else load_pathogenic_scores(options["pathogenic_scores"])

    hpo.load_profiles(profiles)
    evidences_similarity = {}
    for pair, data in evidences.items():
        entity, profile_type = pair.split('_')
        if profile_type == 'HP':
            evidence_profiles = data["prof"]
            #evidence_profiles.transform_keys!{|prof_id, terms| prof_id.to_sym} #not needed in python
            similarities = hpo.compare_profiles(external_profiles= evidence_profiles, sim_type= "lin", bidirectional= False)
            if len(similarities) > 0: evidences_similarity[pair] = similarities

    template = open(str(files('pets.templates').joinpath('evidence_profile.txt'))).read()
    os.makedirs(options["output_folder"], exist_ok=True)
    for profile_id, reference_prof in profiles.items():
        all_candidates = []
        all_genomic_coordinates = {}
        similarity_matrixs = {}
        for pair, ev_profiles_similarity in evidences_similarity.items():
            entity = pair.split('_')[0]
            similarities = ev_profiles_similarity[profile_id]
            candidate_sim_matrix, candidates, candidates_ids = cohort.get_term2term_similarity_matrix(reference_prof, similarities, evidences[pair]["prof"], hpo, 40, 40)
            coords = get_evidence_coordinates(entity, genomic_coordinates, candidates_ids)
            candidate_sim_matrix.insert(0, ['HP'] + candidates_ids)
            if len(pathogenic_scores) > 0: # priorize by pathogenic scores
                candidate_sim_matrix_patho, candidates_patho, candidates_ids_patho = cohort.get_term2term_similarity_matrix(
                    reference_prof, similarities, 
                    evidences[pair]["prof"], hpo, 40, 40, 
                    other_scores = pathogenic_scores, id2label = evidences[pair]["id2lab"])
                if len(candidate_sim_matrix_patho) > 0:
                    candidate_sim_matrix_patho.insert(0, ['HP'] + candidates_ids_patho)	
                    similarity_matrixs[pair + '_path_vars'] = candidate_sim_matrix_patho
                    evidences[pair + '_path_vars'] = evidences[pair]
            if coords == None: continue
            all_candidates.extend(candidates)
            similarity_matrixs[pair] = candidate_sim_matrix
            all_genomic_coordinates.update(coords)
        prof_vars = profile_variants[profile_id]
        hotspots_with_pat_vars = generate_prediction(similarity_matrixs, all_genomic_coordinates, prof_vars)
        make_report(
            profile_id, 
            all_candidates, 
            all_genomic_coordinates, 
            similarity_matrixs, 
            evidences, prof_vars,
            hotspots_with_pat_vars,
            template, options["output_folder"]
        )

#############################################################################################
## METHODS
############################################################################################
def load_pathogenic_scores(path):
    scores = {}
    with open(path) as f:
        for line in f:
            feature, score = line.split("\t")
            scores[feature] = float(score)
    return scores


def get_evidence_coordinates(entity, genomic_coordinates, candidates_ids):
    coords = None
    all_coordinates = genomic_coordinates.get(entity)
    if all_coordinates: coords = {id: coordinates for id, coordinates in all_coordinates.items() if id in candidates_ids}
    #coords = all_coordinates.select{|id, coordinates| candidates_ids.include?(id.to_sym)} if !all_coordinates.nil?
    return coords


def make_report(profile_id, all_candidates, all_genomic_coordinates, similarity_matrixs, 
                            evidences, prof_vars, hotspots_with_pat_vars, template, output):
    var_ids, var_coors = format_variants4report(prof_vars)
    for c in all_candidates:
        c[0] = str(c[0])

    container = {
        "profile_id": profile_id,
        "candidates": all_candidates,
        "genomic_coordinates": {k: c[:2] for k, c in all_genomic_coordinates.items()},
        "similarity_matrixs": similarity_matrixs,
        "evidences": evidences,
        "var_ids": var_ids,
        "var_coordinates": var_coors,
        "hotspot_table": hotspots_with_pat_vars
    }
    report = Py_report_html(container, 'Evidence profile report')
    report.build(template)
    report.write(os.path.join(output, str(profile_id) + '.html'))


def format_variants4report(var_data):
    if var_data == None:
        var_ids, var_coors = None, None
    else:
        var_ids = []
        var_coors = {}
        count = 0
        for chrm, reg in var_data.each():
            var_id = f"var_{count}"
            var_ids.append([var_id, 0])
            var_coors[var_id] = [str(chrm), reg["start"]]
            count += 1
    return var_ids, var_coors


def generate_prediction(similarity_matrixs, all_genomic_coordinates, prof_vars):
    hotspots_with_pat_vars = []
    if prof_vars:
        phen_regions = Genomic_Feature.hash2genomic_feature( all_genomic_coordinates, lambda k, v: v[0:3]+[k] )
        phen_candidates_by_hotspot, phen_genome_hotspots = phen_regions.generate_cluster_regions("reg_overlap", 'A', 0, True)
        genome_matches = phen_genome_hotspots.match(prof_vars)
        hotspot_with_phen_candidates = invert_hash(phen_candidates_by_hotspot) 
        for hotspot_id, pat_vars in genome_matches.items():
            reg = phen_genome_hotspots.region_by_to(hotspot_id)
            coords = [reg["chrm"], reg["start"], reg["stop"]]
            hotspots_with_pat_vars.append([hotspot_id, coords, hotspot_with_phen_candidates[hotspot_id], pat_vars])
        # TODO: see to use original similarities without use top candidates in similarity_matrixs
        # TODO: COMPLETE UNTIL FULL PREDICTOR
    return hotspots_with_pat_vars

# Diseasome functions
#####################

def generate_tripartite_diseasome(adj_list, ontology, layers=["disease","group","parentals"]):
    # Create network
    net = NetAnalyzer(layers)
    for disease, group in adj_list:
        net.add_node(disease,"disease")
        net.add_node(group,"group")
        parentals = ontology.get_ancestors(disease)
        for parental in parentals:
            net.add_node(parental, "parentals")
            net.add_edge(disease, parental)
        net.add_edge(disease, group)
    return net

def get_dis2dclass(diseases, disorder_class, ontology):
    dis_class = set(disorder_class.keys())
    disease2disclass = {}
    for disease in diseases:
        parents = ontology.get_ancestors(disease)
        parents = set(parents) & dis_class
        disease2disclass[disease] = list(parents)
    return disease2disclass

def clean_dis2dclass(disease2disclass, ontology, dependency_map, disorder_class):
    dis_class = set(disorder_class.keys())
    term2ic = get_term2ic(dis_class, ontology)
    for disease, parents in disease2disclass.items():
        # intersect
        #if "MONDO:0021147" in parents or "MONDO:0005071" in parents: continue
        # Syndromics terms: "MONDO:0002254" Multiple
        # direct classify: "MONDO:0045024" Cancer
        # exclusive clasify: "MONDO:0005071" System nervious

        # syndromic indicator
        nmax = 3
        if "MONDO:0002254" in parents: # syndromic
            parents = parents - {"MONDO:0002254"}
            nmax = 2
        # Direct or not classifier
        if "MONDO:0045024" in parents: # One way terms (e.g. cancer)
            parents = ["MONDO:0045024"]
        else:
            parents = just_child_dependencies(parents, dependency_map)
        # Exlusive term
        if "MONDO:0005071" in parents:
                nmax = 1

        number_classes = len(parents)
        if number_classes == 0:
            disclass = "unclasiffied"
        elif number_classes > nmax:
            disclass = "multiple"
        else:
            max_term = ""
            max_ic = 0
            for parent in parents:
                ic = term2ic[parent]
                if ic > max_ic:
                    max_ic = ic
                    max_term = parent 
            disclass = disorder_class[max_term][0]
        disease2disclass[disease] = disclass
    return disease2disclass

def just_child_dependencies(terms, dependency_map):
    filtered_terms = terms
    for term in terms:
        if dependency_map.get(term):
                terms2remove = dependency_map[term]
                filtered_terms = filtered_terms - terms2remove
        continue
    return filtered_terms

def get_dependency_map(terms, ontology):
    term2dep = {}
    for term in terms:
        parentals = ontology.get_ancestors(term)
        dependencies = set(parentals) & terms
        if dependencies:
                term2dep[term] = dependencies
    return term2dep

def get_term2ic(terms, ontology):
    term2ic = {}
    for term in terms:
        term2ic[term] = ontology.get_IC(term)
    return  term2ic 

# collapse functions
#####################

def filter_out_non_leafs_nodes(mondo_terms, ontology):
    filtered = []
    for term in mondo_terms:
        childs = ontology.get_descendants(term)
        if not childs or len(childs) == 0:
            filtered.append(term)
    return filtered

def get_parent_and_childs_nodes_dict(mondo_terms, ontology):
    parent_to_childs_dict = {}
    for term in mondo_terms:
        parents = ontology.get_direct_related(term, relation="ancestor")
        for parent in parents:
            query = parent_to_childs_dict.get(parent)
            if query is None:
                parent_to_childs_dict[parent] = [term]
            else:
                parent_to_childs_dict[parent].append(term)

    parent_to_childs_dict = {parent: childs for parent, childs in parent_to_childs_dict.items() if len(childs) > 1}
    return parent_to_childs_dict

def get_txt_to_txt_similarities(parent_to_childs_txt, rm_char = ""):
    similarities = {}
    for parent, txts in parent_to_childs_txt.items():
        if len(txts) > 1: similarities[parent] = similitude_network(txts, charsToRemove = rm_char)
    return similarities

def get_thresholded_childs_to_parents_dict(similarities, threshold, ontology, rm_char = "", txt_to_term = None, uniq_parent = False):
    terms_to_parents_collapsed = {}

    for parent, childs in similarities.items():
        for child1, other_childs in childs.items():
            for child2, similarity in other_childs.items():     
                # Watch out: This is selecting for childs with at least one similarity beyond the threshold.  
                if similarity >= threshold:
                    # Just neccesary the parent id. translated_parent = ontology.translate_name(parent)
                    if txt_to_term:
                        child1_term = txt_to_term[child1]
                        child2_term = txt_to_term[child2]

                    if terms_to_parents_collapsed.get(child1_term) == None: 
                        terms_to_parents_collapsed[child1_term] = [parent]
                    else:
                        terms_to_parents_collapsed[child1_term].append(parent)
                
                    if terms_to_parents_collapsed.get(child2_term) == None: 
                        terms_to_parents_collapsed[child2_term] = [parent]
                    else:
                        terms_to_parents_collapsed[child2_term].append(parent)
    if uniq_parent: 
        terms_to_parents_collapsed = get_collapsed_with_unique_parents(ontology, terms_to_parents_collapsed, rm_char)
    else:
        terms_to_parents_collapsed = { child: list(set(parents)) for child, parents in terms_to_parents_collapsed.items() }

    return terms_to_parents_collapsed


def get_collapsed_with_unique_parents(ontology, terms_to_parents_collapsed, rm_char = ""):
    collapsed_with_unique_parents = {}
    for child, parents in terms_to_parents_collapsed.items():
        parents = list(set(parents))
        parents_depth = [ontology.term_paths[parent]["largest_path"] for parent in parents]
        max_depth_indexes = [i for i, x in enumerate(parents_depth) if x == max(parents_depth)]
        if len(max_depth_indexes) == 1:
            collapsed_with_unique_parents[child] = [parents[max_depth_indexes[0]]]
        else:
            deepest_parents = [parents[i] for i in max_depth_indexes]
            translated_child = ontology.translate_id(child)
            translated_parents = [ontology.translate_id(parent) for parent in deepest_parents]
            similarities = similitude_network([translated_child] + translated_parents, charsToRemove = rm_char)
            collapsed_with_unique_parents[child] = [ontology.translate_name(max(similarities[translated_child].items(), key=lambda x: x[1])[0])]
    return collapsed_with_unique_parents