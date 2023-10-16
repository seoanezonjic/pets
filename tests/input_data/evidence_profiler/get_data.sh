#! /usr/bin/env bash
. ~soft_bio_267/initializes/init_R

mkdir -p processed_data
zgrep 'HP:' raw_data/gene_phenotype.all.tsv.gz | grep 'NCBITaxon:9606' | aggregate_column_data.rb -i - -x 1 -a 4 | sed 's/ (human)//g' | head -n 10000 > processed_data/gene_phenotype
cut -f 1 processed_data/gene_phenotype > processed_data/gene_list
get_genomic_coordinates.R -i processed_data/gene_list -l region -o processed_data/gene_coordinates
awk 'FS="\t" {if($5 != "NA") print $1 "\t" $3 "\t" $5}' processed_data/gene_coordinates > processed_data/gene_coordinates_clean
