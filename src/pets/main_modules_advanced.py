import os, glob, json, requests, pickle
from importlib.resources import files
import urllib.parse
import warnings
from collections import Counter

import pandas as pd
from py_cmdtabs.cmdtabs import CmdTabs
from py_exp_calc.exp_calc import invert_hash, uniq
from py_report_html.py_report_html import Py_report_html
from py_semtools.ontology import Ontology
from py_semtools.sim_handler import similitude_network
import pets.report_pets
from pets.parsers.cohort_parser import Cohort_Parser
from pets.cohort import Cohort
from pets.io import load_profiles, load_variants, load_evidences, load_index, parse_morbid_omim, list2dic
from pets.genomic_features import Genomic_Feature

# https://setuptools.pypa.io/en/latest/userguide/datafiles.html
HPO_FILE = str(files('pets.external_data').joinpath('hp.json'))
MONDO_FILE = str(files('pets.external_data').joinpath('mondo.obo'))
IC_FILE = str(files('pets.external_data').joinpath('uniq_hpo_with_CI.txt'))
GENCODE = str(files('pets.external_data').joinpath('gencode.v43.basic.annotation.gtf.gz'))

def main_pedigree_analysis(opts):
    from pets.pedigree_analysis import PedigreeAnalyzer
    options = opts
    analyzer = PedigreeAnalyzer()
    #analyzer.de_novo_tolerant = False
    analyzer.load_pedigree(options["pedigree_file"])
    vcf_ref = analyzer.load_vcf_merged(options["merged_vcf"])

    if options.get("variant_prefilter"):
        print(f"Applying variant prefilters: {options['variant_prefilter']}")
        for variant_filter in options["variant_prefilter"]:
            if variant_filter not in ["de_novo", "homozygous", "compound_het"]:
                raise ValueError(f"Invalid variant filter: {variant_filter}. Must be one of 'de_novo', 'homozygous', 'compound_het'.")
            else:
                print(f"Applying variant filter: {variant_filter}")
                analyzer.filter_variants_by_type(variant_filter)
    
    analyzer.build_matrixes()
    if options.get("variant_attrs"):
        analyzer.load_variant_attributes(options["variant_attrs"])
        analyzer.build_var2attr_matrix()

    ref_vars = set(vcf_ref.keys())
    analyzer.filter_by_inheritance("ARC", "S1132-2")#options["desired_moi"])
    merged_vars = set(analyzer.variant_ids)
    print(f"Reference variants: {len(ref_vars)}")
    print(f"Variants after filtering by inheritance: {len(merged_vars)}")
    print(f"Variants after filtering by inheritance and reference: {len(merged_vars - ref_vars)}")
    print(f"Example of variants after filtering by inheritance and reference: {list(merged_vars - ref_vars)[:10]}")
    #if options.get("de_novo_tolerant"): analyzer.de_novo_tolerant = True
    #analyzer.remove_de_novo_variants()
    #analyzer.analyze()
    #analyzer.generate_report()