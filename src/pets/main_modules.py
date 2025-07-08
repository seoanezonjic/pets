import os
import glob, json
from importlib.resources import files
import urllib.parse
import requests
from collections import Counter

import numpy as np
import pandas as pd
from py_report_html import Py_report_html
import pets
import pets.report_pets
from pets.cohort_analyser_methods import *
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort
from pets.io import write_tabulated_data, load_profiles, load_variants, load_evidences, load_index, parse_morbid_omim, list2dic
from pets.genomic_features import Genomic_Feature
from pets.parsers.reference_parser import Reference_parser
from pets.parsers.coord_parser import Coord_Parser
from pets.genomic_prioritizer import (
    GenomicPrioritizer,
    AimarrvelPrioritizer,
    LiricalPrioritizer,
    Phen2GenePrioritizer,
    GadoPrioritizer,
    ExomiserPrioritizer,
    PhenogeniusPrioritizer,
    DefaultGenomicPrioritizer,
    MetaGenomicPrioritizer,
    HeuristicModel,
    XGBoostRankerModel,
    LogisticRegressionModel
)
from py_exp_calc.exp_calc import invert_hash, uniq
from py_semtools.ontology import Ontology
from py_semtools.sim_handler import similitude_network
from py_cmdtabs.cmdtabs import CmdTabs
import pickle

# https://setuptools.pypa.io/en/latest/userguide/datafiles.html
HPO_FILE = str(files('pets.external_data').joinpath('hp.json'))
MONDO_FILE = str(files('pets.external_data').joinpath('mondo.obo'))
IC_FILE = str(files('pets.external_data').joinpath('uniq_hpo_with_CI.txt'))
GENCODE = str(files('pets.external_data').joinpath('gencode.v43.basic.annotation.gtf.gz'))


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
    similarities = get_txt_to_txt_similarities(parent_to_childs_txt, rm_char=options["rm_char"], algorithm = options['sim_algorithm'])
    terms_to_parents_collapsed = get_thresholded_childs_to_parents_dict(similarities, threshold = options.get("threshold"),
     ontology=ontology, rm_char = options["rm_char"], txt_to_term = txt_to_term, uniq_parent = options["uniq_parent"], algorithm = options['sim_algorithm'])

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


def main_paco_translator(opts):
    options = vars(opts)
    hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
    Cohort.load_ontology("hpo", hpo_file, options.get("excluded_hpo"))
    Cohort.act_ont = "hpo"

    patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(options)

    if options.get("clean_PACO"):
        removed_terms, removed_profiles = patient_data.check(hard=options["hard_check"])
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
    rejected_hpos_C, rejected_patients_C = patient_data.check(hard=options["hard_check"])
    patient_data.link2ont(Cohort.act_ont)

    vcf_index = None
    attr_index = None
    if options.get("vcf_index"): vcf_index = load_index(options["vcf_index"])
    if options.get("attr_index"): attr_index = load_index(options["attr_index"])
    patient_data.export_phenopackets(options["output_folder"], options["genome_assembly"], 
        vcf_index= vcf_index, attr_index=attr_index, attr_name = options.get("attr_name"), v2=options['v2'])


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
    hpo = Cohort.ont['hpo']

    #Temporal Ontology Object to save the common profile for each cluster when --detailed_clusters is active and detailed_cluster_yaxis is set to 'cluster'
    temporal_hpo = None
    if opts['detailed_clusters'] and opts['detailed_cluster_yaxis'] == 'cluster':
        temporal_hpo = Cohort.load_ontology("hpo", hpo_file, opts.get("excluded_hpo"), inplace=False)

    opts['check'] = True
    patient_data, rejected_hpos, rejected_patients = Cohort_Parser.load(opts)
    with open(rejected_file, 'w') as f: f.write("\n".join(rejected_patients))

    patient_data.link2ont(Cohort.act_ont) # TODO: check if method load should call to this and use the semtools checking methods (take care to only remove invalid terms)

    patient_data.get_profile_redundancy() # GEt term redundancy BEFORE cleaning
    patient_data.check(hard=opts["hard_check"])
    patient_data.link2ont(Cohort.act_ont) #Now that we have calculate profiles redundancy, we synchronize the cleaned profiles from HPO object
    
    hpo.get_profiles_terms_frequency() # hpo CODE, freq
    with open(hpo_frequency_file, 'w') as f:
      for hpo_code, freq in hpo.dicts['term_stats'].items(): f.write(f"{hpo_code}\t{freq}\n")

    suggested_childs, fraction_terms_specific_childs = patient_data.compute_term_list_and_childs(file = detailed_profile_evaluation_file)

    phenotype_ic = patient_data.get_ic_analysis(freq_type = opts['ic_stats'], 
        ic_file = os.environ['ic_file'] if os.environ.get('ic_file') else IC_FILE)

    all_ics, prof_lengths, clust_by_chr, top_clust_phen, multi_chr_clusters = patient_data.process_dummy_clustered_patients(opts, phenotype_ic, temp_folder = temp_folder)

    summary_stats = get_summary_stats(patient_data, rejected_patients, hpo.ics['resnik_observed'], fraction_terms_specific_childs, rejected_hpos)

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
    sortByPhens = None
    if opts['detailed_cluster_yaxis'] == 'cohort_sort':
        phens_ocurrences = Counter()
        for profile in patient_data.profiles.values(): phens_ocurrences.update([hpo.translate_id(term) for term in profile])
        sortByPhens = sortYaxisByPhens(phens_ocurrences)
    reference_profiles = None
    if opts.get('reference_profiles') != None: reference_profiles = load_profiles(opts['reference_profiles'], Cohort.get_ontology('hpo'))
    template = str(files('pets.templates').joinpath('cluster_report.txt'))
    clustering_data = get_semantic_similarity_clustering(opts, patient_data, reference_profiles, temp_folder, template, temporal_hpo, ySortFunc=sortByPhens)

    #----------------------------------
    # GENERAL COHORT ANALYZER REPORT
    #----------------------------------
    new_cluster_phenotypes = get_top_dummy_clusters_stats(top_clust_phen)

    container = {
      'temp_folder' : temp_folder,
      'summary_stats' : summary_stats,
      'clustering_methods' : opts['clustering_methods'],
      'all_cnvs_length' : [ [l] for l in all_cnvs_length ],
      'all_sor_length' : [ [l] for l in all_sor_length ],
      'new_cluster_phenotypes' : len(new_cluster_phenotypes),
      'dummy_cluster_chr_data' : dummy_cluster_chr_data,
      'dummy_ic_data' : format_cluster_ic_data(all_ics, prof_lengths, opts['clusters2graph']),
      'chr_sizes': chr_sizes,
      'ontology': Cohort.ont['hpo'],
      'ontoplot_mode': opts['ontoplot_mode']
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
            candidate_sim_matrix, candidates, candidates_ids = hpo.get_term2term_similarity_matrix(reference_prof, similarities, evidences[pair]["prof"], 40, 40)
            coords = get_evidence_coordinates(entity, genomic_coordinates, candidates_ids)
            candidate_sim_matrix.insert(0, ['HP'] + candidates_ids)
            if len(pathogenic_scores) > 0: # priorize by pathogenic scores
                candidate_sim_matrix_patho, candidates_patho, candidates_ids_patho = hpo.get_term2term_similarity_matrix(
                    reference_prof, similarities, 
                    evidences[pair]["prof"], 40, 40, 
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

def main_report_prioritizer(opts):
    options = vars(opts)
    print("AQUI ESTOY CHAVAL"*30)
    
    prioritizer = {}
    for prioritizer_type, path2folder_results in options["prioritizers"].items():
        if prioritizer_type == "phen2gene":
            prioritizer[(prioritizer_type, path2folder_results)] = Phen2GenePrioritizer()
        elif prioritizer_type == "gado":
            prioritizer[(prioritizer_type, path2folder_results)] = GadoPrioritizer()
        elif prioritizer_type == "phenogenius":
            prioritizer[(prioritizer_type, path2folder_results)] = PhenogeniusPrioritizer()
        elif prioritizer_type == "exomiser":
            prioritizer[(prioritizer_type, path2folder_results)] = ExomiserPrioritizer()
        elif prioritizer_type == "aimarrvel":
            prioritizer[(prioritizer_type, path2folder_results)] = AimarrvelPrioritizer()
        elif prioritizer_type == "lirical":
            prioritizer[(prioritizer_type, path2folder_results)] = LiricalPrioritizer()
        elif prioritizer_type == "default":
            prioritizer[(prioritizer_type, path2folder_results)] = DefaultGenomicPrioritizer()
        else:
            raise Exception(f"Unknown prioritizer: {prioritizer}")
        if options["benchmark_type"] == "gene" or options["benchmark_type"] == "both":
            prioritizer[(prioritizer_type, path2folder_results)].post_process_results_genes(path2folder_results, 
                             write_tmp=options["write_tmp"], read_tmp=options["read_tmp"])
        elif options["benchmark_type"] == "variant" or options["benchmark_type"] == "both":
            prioritizer[(prioritizer_type, path2folder_results)].post_process_results_variants(path2folder_results, 
                             write_tmp=options["write_tmp"], read_tmp=options["read_tmp"])

    if not options["write_tmp"]:
        if options["integrated_report"]:
            if len(prioritizer.keys()) > 1:
                metaprioritizer = MetaGenomicPrioritizer(prioritizer)
                metaprioritizer.get_features(type=options["benchmark_type"])
                metaprioritizer.test_patients = metaprioritizer.get_all_patients()
                print(metaprioritizer.test_patients)
                metaprioritizer.model = HeuristicModel()
                metaprioritizer.predict_test(type=options["benchmark_type"])
                prio_table, quantitative_feature, qualitative_feature = metaprioritizer.get_combined_results(type=options["benchmark_type"])

            else:
                first_prioritizer = list(prioritizer.values())[0]
                prio_table, quantitative_feature, qualitative_feature = first_prioritizer.get_combined_results(options["benchmark_type"])

            container = {
                "quantitative": quantitative_feature,
                "qualitative": qualitative_feature,
                "prio_table": prio_table
            }
            template="integrated_by_patient_prioreport.txt"
        else:
            first_prioritizer = list(prioritizer.values())[0]

            quantitative_feature = list(first_prioritizer.quant_features_idx.values())[0]
            qualitative_feature = list(first_prioritizer.qual_features_idx.values())[0]
            if options["benchmark_type"] == "gene" or options["benchmark_type"] == "both":
                prio_table = list(first_prioritizer.patient2gene_results.values())[0] 
            elif options["benchmark_type"] == "variant" or options["benchmark_type"] == "both":
                prio_table = list(first_prioritizer.patient2variant_results.values())[0] 

            container = {
                "quantitative": quantitative_feature,
                "qualitative": qualitative_feature,
                "prio_table": prio_table
            }
            template="individual_prioreport.txt"

        report = Py_report_html(container)
        report.build(open(str(files('pets.templates').joinpath(template))).read())
        report.write(options["output_file"] + '.html')

def main_meta_prioritizer(opts):
    options = vars(opts)
    
    # loading every prioritizer
    prioritizer = {}
    for prioritizer_type, path2folder_results in options["prioritizers"].items():
        print("ST 1 - Pass to process step\n--------------------------\n---------------------------")
        if prioritizer_type == "phen2gene":
            prioritizer[(prioritizer_type, path2folder_results)] = Phen2GenePrioritizer()
        elif prioritizer_type == "gado":
            prioritizer[(prioritizer_type, path2folder_results)] = GadoPrioritizer()
        elif prioritizer_type == "phenogenius":
            prioritizer[(prioritizer_type, path2folder_results)] = PhenogeniusPrioritizer()
        elif prioritizer_type == "exomiser":
            prioritizer[(prioritizer_type, path2folder_results)] = ExomiserPrioritizer()
        elif prioritizer_type == "aimarrvel":
            prioritizer[(prioritizer_type, path2folder_results)] = AimarrvelPrioritizer()
        elif prioritizer_type == "lirical":
            prioritizer[(prioritizer_type, path2folder_results)] = LiricalPrioritizer()
        elif prioritizer_type == "default":
            prioritizer[(prioritizer_type, path2folder_results)] = DefaultGenomicPrioritizer()
        else:
            raise Exception(f"Unknown prioritizer: {prioritizer}")
        if options["benchmark_type"] == "gene" or options["benchmark_type"] == "both":
            prioritizer[(prioritizer_type, path2folder_results)].post_process_results_genes(path2folder_results, 
                             write_tmp=options["write_tmp"], read_tmp=options["read_tmp"])
        elif options["benchmark_type"] == "variant" or options["benchmark_type"] == "both":
            prioritizer[(prioritizer_type, path2folder_results)].post_process_results_variants(path2folder_results, 
                             write_tmp=options["write_tmp"], read_tmp=options["read_tmp"])


    # Create MetaPrioritizer and getting the features
    metaprioritizer = MetaGenomicPrioritizer(prioritizer)
    #metaprioritizer.get_features(type=options["benchmark_type"])

    if options["labels"]:
        patient_labels = CmdTabs.load_input_data(options["labels"])
        metaprioritizer.load_patient_labels(patient_labels)
    
    # load the model
    if options["model_type"] == "heuristic":
        metaprioritizer.model = HeuristicModel()
    elif options["model_type"] == "xgboost":
        metaprioritizer.model = XGBoostRankerModel()
    elif options["model_type"] == "logistic_regression":
        metaprioritizer.model = LogisticRegressionModel()
    else:
        if options["model_type"] is not None:
            with open(options["model_type"], "rb") as f:
                metaprioritizer.model = pickle.load(f)

    metaprioritizer.get_features(type=options["benchmark_type"], dropna=options["dropna"])
    
    if options["mode"] == "predict":
        metaprioritizer.test_patients = metaprioritizer.get_all_patients(type=options["benchmark_type"])
        metaprioritizer.predict_test(type=options["benchmark_type"])
    elif options["mode"] == "train":
        metaprioritizer.train_patients = metaprioritizer.get_all_patients(type=options["benchmark_type"])
        metaprioritizer.train_model(type=options["benchmark_type"])
    elif options["mode"] == "train_predict":
        metaprioritizer.split_patients(type="gene", test_size=0.3, random_state=42)
        metaprioritizer.train_model(type=options["benchmark_type"])
        metaprioritizer.predict_test(type=options["benchmark_type"])

    if options["mode"] == "predict" or options["mode"] == "train_predict":
        results_id = "patient2gene_results" if options["benchmark_type"] == "gene" else "patient2variant_results"
        os.makedirs(options["output_file"], exist_ok=True)
        for patient, results in getattr(metaprioritizer, results_id).items():
            results.to_csv(os.path.join(options["output_file"],patient), index=False, sep="\t")

    if options["mode"] == "train":
        with open(os.path.join(options["output_file"],f"trained_{options["model_type"]}_model.pkl"), "wb") as f:
            pickle.dump(metaprioritizer.model, f)

#############################################################################################
## METHODS
############################################################################################
def sortYaxisByPhens(phens_ocurrences):
    def sortByPhens(phen):
        return phens_ocurrences[phen]
    return sortByPhens


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
    from NetAnalyzer import NetAnalyzer

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

def get_txt_to_txt_similarities(parent_to_childs_txt, rm_char = "", algorithm = 'white'):
    similarities = {}
    for parent, txts in parent_to_childs_txt.items():
        if len(txts) > 1: similarities[parent] = similitude_network(txts, charsToRemove = rm_char, algorithm = algorithm)
    return similarities

def get_thresholded_childs_to_parents_dict(similarities, threshold, ontology, rm_char = "", txt_to_term = None, uniq_parent = False, algorithm = 'white'):
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
        terms_to_parents_collapsed = get_collapsed_with_unique_parents(ontology, terms_to_parents_collapsed, rm_char, algorithm = algorithm)
    else:
        terms_to_parents_collapsed = { child: list(set(parents)) for child, parents in terms_to_parents_collapsed.items() }

    return terms_to_parents_collapsed


def get_collapsed_with_unique_parents(ontology, terms_to_parents_collapsed, rm_char = "", algorithm = 'white'):
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
            similarities = similitude_network([translated_child] + translated_parents, charsToRemove = rm_char, algorithm = algorithm)
            collapsed_with_unique_parents[child] = [ontology.translate_name(max(similarities[translated_child].items(), key=lambda x: x[1])[0])]
    return collapsed_with_unique_parents

def main_phenPatMaster(opts):
    # TODO: Refactor this creating a parser and writer of phenopackets and integrate Cohort class to perform this operations
    hpo_file = os.environ['hpo_file'] if os.environ.get('hpo_file') else HPO_FILE
    ontology = Ontology(file= hpo_file, load_file = True) 

    phenopacket_files=glob.glob(os.path.join(opts.input_folder, '*.json'))
    count = 1
    pp_dict = dict()
    index = []
    for pp_path in phenopacket_files:
        new_id = "pp"+str(count)
        phenopacket = json.loads(open(pp_path).read().encode("utf-8"))
        if opts.overwrite_id:
            pp_dict[phenopacket['id']] = new_id
            phenopacket['id'] = new_id
            phenopacket['subject']['id'] = new_id
            for inter in phenopacket['interpretations']:
                inter['id'] = new_id
                genomic_inter= inter['diagnosis']['genomicInterpretations']
                for gen_inter in genomic_inter:
                    gen_inter['subjectOrBiosampleId'] = new_id
        count += 1
        
        if opts.clean_phen or opts.output_file_index != None:
            pp_phens = phenopacket['phenotypicFeatures']
            phens = []
            neg_phens = [] # HPOs that has NOT present the patient
            for ph in pp_phens:
                hp_code = ph["type"]['id']
                if ph.get("excluded"):
                    neg_phens.append(hp_code)
                else:
                    phens.append(hp_code)

            if opts.clean_phen:
                phens = ontology.clean_profile_hard(phens)
                neg_phen_defs, _ = ontology.check_ids(neg_phens)
                pp_phens = []
                for hp in phens: pp_phens.append({'type':{'id': hp , 'label': ontology.translate_id(hp) }})
                for hp in neg_phen_defs: pp_phens.append({'type':{'id': hp , 'label': ontology.translate_id(hp) }, 'excluded': True })
                phenopacket['phenotypicFeatures'] = pp_phens

            if opts.output_file_index != None:
                for interp in phenopacket["interpretations"]:
                    genomic_inter = interp['diagnosis']['genomicInterpretations']
                    for gen_interp in genomic_inter:
                        variant = gen_interp['variantInterpretation']['variationDescriptor'].get('vcfRecord')
                        if variant != None:
                            index.append([phenopacket['id'], phens, variant['chrom'], variant['pos'], variant['pos']])
                        else:
                            index.append([phenopacket['id'], phens, "", "", ""])

        if opts.output_file_index != None:
            with open(opts.output_file_index, "w") as outfile:
                for p_id, phens, chrom, start,stop in index: outfile.write(f"{p_id}\t{','.join(phens)}\t{chrom}\t{start}\t{stop}\n")

        json_object = json.dumps(phenopacket, indent=4)
        if opts.overwrite_file_name:
            pp_name = phenopacket['id'] + '.json'
        else:
            pp_name = os.path.basename(pp_path)
        with open(os.path.join(opts.output_folder, pp_name), "w") as outfile: outfile.write(json_object)

    if opts.overwrite_id:
        with open(os.path.join(opts.output_folder, 'pp_dict.txt'), "w") as outfile:
            for old_id, new_id in pp_dict.items():
                outfile.write(old_id+"\t"+new_id+"\n")


def main_vcf2effects(opts):
    import logging
    from varcode import load_vcf

    for _ in logging.root.manager.loggerDict: # to disable info logger in pyensembl and varcode
        logging.getLogger(_).setLevel(logging.CRITICAL)

    vcfVariants = load_vcf(opts.input_vcf, allow_extended_nucleotides= True, genome=opts.genome)
    var_effects = []
    for variant in vcfVariants:
        effects = variant.effects().drop_silent_and_noncoding()
        if len(effects) > 0:
            top_effect = effects.top_priority_effect()
            var = top_effect.variant
            var_effects.append([var.contig, var.start , var.ref, var.alt, top_effect.transcript_name, top_effect.short_description])

    with open(opts.output, 'w') as f:
        for contig, start, ref, alt, t_name, effect in var_effects:
            f.write(f"{contig}\t{start}\t{ref}\t{alt}\t{t_name}\t{effect}\n")