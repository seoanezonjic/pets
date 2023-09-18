#! /usr/bin/env python

#require 'fileutils'
import argparse, os, sys
from py_semtools import Ontology

ROOT_PATH = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(ROOT_PATH, '..'))

sys.path.insert(0, os.path.join(ROOT_PATH, '..', '..', 'py_report_html'))

from py_report_html import Py_report_html
import pets
from pets.cohort import Cohort
from pets.cohort_analyser_methods import load_profiles
from pets.evidence_profiler_methods import load_variants, load_evidences
from pets.genomic_features import Genomic_Feature
#from pets.io import load_profiles, load_variants, load_evidences



ROOT_PATH = os.path.dirname(__file__)
CONSTANTS_PATH = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'constants.py'))
with open(CONSTANTS_PATH) as infile:
    exec(infile.read())

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


def invert_hash(h):
	new_h = {}
	for k, vals in h.items():
		for v in vals:
			query = new_h.get(v)
			if query == None:
				new_h[v] = [k]
			else:
				query.append(k)
	return new_h

############################################################################################
## OPTPARSE
############################################################################################
parser = argparse.ArgumentParser(description=f'Usage: {os.path.basename(__file__)} [options]')

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

opts = parser.parse_args()
options = vars(opts)


#############################################################################################
## MAIN
############################################################################################

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


template = open(os.path.join(REPORT_FOLDER, 'evidence_profile.txt')).read()
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
