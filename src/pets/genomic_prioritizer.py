
import os
import pandas as pd
import pickle
import logging
import re
import json
from py_exp_calc.exp_calc import get_rank_metrics

logger = logging.getLogger(__name__)

def parse_hgnc_data() -> pd.DataFrame:
    base_path = os.path.dirname(__file__) 
    hgnc_path = os.path.join(base_path, "external_data", "hgnc_complete_set.txt")

    df = pd.read_csv(hgnc_path, sep="\t", dtype={"omim_id": str, "entrez_id":str})

    df = df.rename(columns={
        "hgnc_id": "hgnc_id",
        "symbol": "gene_symbol",
        "ensembl_gene_id": "ensembl_id",
        "entrez_id": "entrez_id",
        "refseq_accession": "refseq_accession",
        "prev_symbol": "previous_symbol_raw",
    })

    df["prev_symbols"] = df["previous_symbol_raw"].fillna("").apply(
        lambda x: [s.strip('"') for s in x.split("|")] if x else []
    )

    return df[["hgnc_id", "gene_symbol", "ensembl_id", "entrez_id", "refseq_accession", "prev_symbols"]]



class GenomicPrioritizer:
    #identifier_map = create_gene_identifier_map()
    identifier_map = parse_hgnc_data()
    print(identifier_map)
    gene_identifiers = ["gene_symbol", "ensembl_id"]

    def __init__(self):
        # For variant: rank|score|chromosome|start|end|ref|alt
        self.patient2variant_results = {}
        # For Gene: rank|score|gene_symbol|gene_identifier(ensbl)
        self.patient2gene_results = {}
        self.priot_type = ["gene", "variant"]
        self.identifier_translator = {}
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def load_identifier_tanslator(self, source_id, target_id):
        self.identifier_translator[(source_id, target_id)] = {}
        for _, row in type(self).identifier_map[[source_id, target_id]].iterrows():
            self.identifier_translator[(source_id, target_id)][row[source_id]] = row[target_id]
    
    # def post_process(self, raw_results_dir, write_tmp=None, read_tmp=False):
    #     # Placeholder for the implementation
    #     # This function will read the raw results, process them, and save the standardized output.
    #     if "gene" in self.priot_type:
    #         # Process gene results
    #         self.post_process_results_genes(raw_results_dir, write_tmp, read_tmp)
    #     if "variant" in self.priot_type:
    #         # Process variant results
    #         self.post_process_results_variants(raw_results_dir, write_tmp, read_tmp)

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        raise NotImplementedError("To implement this method")

    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        raise NotImplementedError("To implement this method")

class DefaultGenomicPrioritizer(GenomicPrioritizer):

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        print("")
        pass

    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        print("")
        pass

class Phen2GenePrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene"] 
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f] = pickle.load(pf)
            else:
                df = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                self.patient2gene_results[f] = df
                # Process the data
                self.patient2gene_results[f] = self.get_result_gene(df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
            self.quant_features_idx[f] = None
            self.qual_features_idx[f] = [4]
    
    def get_result_gene(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[]):
        """
        Process the data from the raw results directory.
        """
        processed_data = pd.DataFrame()
        processed_data["rank"] = df.iloc[:, rank]
        processed_data["score"] = df.iloc[:, score]
        for key, value in gene.items():
            processed_data[key] = df.iloc[:, value]
        self.load_identifier_tanslator("gene_symbol", "ensembl_id")
        processed_data["ensembl_id"] = [
            self.identifier_translator[("gene_symbol", "ensembl_id")].get(row, None)
            for row in processed_data["gene_symbol"]
        ]
        colnames = df.columns.values.tolist()
        for feature in quali_feature:
            processed_data[colnames[feature]] = df.iloc[:, feature]
        pass
    
        return processed_data
    

class GadoPrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene"] 
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f] = pickle.load(pf)
            else:
                df = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                self.patient2gene_results[f] = df
                # Process the data
                self.patient2gene_results[f] = self.get_result_gene(df, rank=2, gene={"gene_symbol": 1, "ensembl_id":0}, score=3, quali_feature=list(range(4,len(df.columns.values))), quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
            self.quant_features_idx[f] = None
            self.qual_features_idx[f] = list(range(4,self.patient2gene_results[f].shape[1]))

    def get_result_gene(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[]):
        """
        Process the data from the raw results directory.
        """
        processed_data = pd.DataFrame()
        processed_data["rank"] = df.iloc[:, rank]
        processed_data["score"] = df.iloc[:, score]
        for key, value in gene.items():
            processed_data[key] = df.iloc[:, value]
        colnames = df.columns.values.tolist()
        for feature in quali_feature:
            processed_data[colnames[feature]] = df.iloc[:, feature]
        pass
    
        return processed_data
    
class PhenogeniusPrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene"] 
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f] = pickle.load(pf)
            else:
                df = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                self.patient2gene_results[f] = df
                # Process the data
                self.patient2gene_results[f] = self.get_result_gene(df, rank=2, gene={"gene_symbol": 1}, score=3)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
            self.quant_features_idx[f] = list(range(4,self.patient2gene_results[f].shape[1]-1))
            self.qual_features_idx[f] = [self.patient2gene_results[f].shape[1]]
    
    def _get_hpos_scores(self, row, col=4):
        hpos_scores = re.sub("\[|\]|{|}","",row[col]).split(",")
        hpos = []
        scores = []
        if hpos_scores == ['']:
            return hpos, scores
        else:
            for hpo_score in hpos_scores:
                hpo, score = hpo_score.split(": ")
                hpo = hpo.replace("'", "")
                hpo = hpo.strip()
                hpos.append(hpo)
                scores.append(score)
        return hpos, scores

    def get_hpo_table(self, row):
        hpos, scores = self._get_hpos_scores(row, "hpo_implicated")
        hpo_description, _ = self._get_hpos_scores(row, "hpo_description_implicated")
        hpo2desc_score = {}
        for i in range(len(hpos)):
            hpo2desc_score[hpos[i]] = (hpo_description[i], scores[i])
        return hpo2desc_score, set(hpos)

    def get_result_gene(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[], quant_feature=[]):
        """
        Process the data from the raw results directory.
        """
        processed_data = pd.DataFrame()
        processed_data["rank"] = df.iloc[:, rank]
        processed_data["score"] = df.iloc[:, score]
        for key, value in gene.items():
            processed_data[key] = df.iloc[:, value]
        self.load_identifier_tanslator("gene_symbol", "ensembl_id")
        processed_data["ensembl_id"] = [
            self.identifier_translator[("gene_symbol", "ensembl_id")].get(row, None)
            for row in processed_data["gene_symbol"]
        ]
        rows = df.to_dict(orient="records")
        row2hpo_desc_score = {}
        all_hpos = set()
        for idx, row in enumerate(rows):
            hpo2desc_score, hpos = self.get_hpo_table(row)
            row2hpo_desc_score[idx] = hpo2desc_score
            all_hpos = all_hpos.union(hpos)
        all_hpos = sorted(list(all_hpos))

        hpo_score_table = []
        for idx in range(0,len(rows)):
            hpo_score_row = []
            hpo2desc_score = row2hpo_desc_score[idx]
            for hpo in all_hpos:
                desc_score = hpo2desc_score.get(hpo)
                if desc_score:
                    hpo_score_row.append(desc_score[1])
                else:
                    hpo_score_row.append("0")
            hpo_score_table.append(hpo_score_row)
        
        hpo_features_df = pd.DataFrame(hpo_score_table, columns=all_hpos)
        processed_data = pd.concat([processed_data, hpo_features_df], axis=1)
        processed_data["phenotype_specificity"] = df.iloc[:, 6]
        return processed_data

class AimarrvelPrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene", "variant"] 
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f, "prediction", "conf_4Model", "integrated", f"{f}_integrated.csv")
                df = pd.read_csv(file_path, sep=",")
                print(df)
                self.patient2gene_results[f] = df
                # Process the data
                self.patient2gene_results[f] = self.get_result_gene(df)
                print(self.patient2gene_results[f])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
            self.quant_features_idx[f] = None
            self.qual_features_idx[f] = None

    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2variant_results[f] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f, "prediction", "conf_4Model", "integrated", f"{f}_integrated.csv")
                df = pd.read_csv(file_path, sep=",")
                self.patient2variant_results[f] = df
                # Process the data
                self.patient2variant_results[f] = self.get_result_variant(df)
                print(self.patient2variant_results[f])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2variant_results[f], pf)
            self.quant_features_idx[f] = None
            self.qual_features_idx[f] = None

    def get_result_gene(self, df, rank=107, gene={"gene_symbol": 111, "ensemble_id": 112}, score=105, quali_feature=[], quant_feature=[]):
        processed_data = pd.DataFrame()
        processed_data["rank"] = df.iloc[:, rank]
        processed_data["score"] = df.iloc[:, score]
        for key, value in gene.items():
            processed_data[key] = df.iloc[:, value]
        processed_data = processed_data[processed_data["gene_symbol"] != "-"]
        processed_data = processed_data.drop_duplicates() # This part could be controversial
        return processed_data

    def get_result_variant(self, df):
        processed_data = pd.DataFrame()
        processed_data["rank"] = df.iloc[:, 107]
        processed_data["score"] = df.iloc[:, 104]

        # Split variant string into parts
        split_var = df.iloc[:, 0].str.split("-", expand=True)
        chrom = split_var[0]
        pos = split_var[1].astype(int)
        ref = split_var[2]
        alt = split_var[3]
        processed_data["contigName"] = chrom
        processed_data["start"] = pos
        processed_data["end"] = pos + ref.str.len() - 1
        processed_data["ref"] = ref
        processed_data["alt"] = alt
        processed_data["varId"] = (
            chrom + "-" +
            pos.astype(str) + "-" +
            processed_data["end"].astype(str) + "-" +
            ref + "-" +
            alt
        )
        processed_data = processed_data[["rank", "score", "varId", "contigName", "start", "end", "ref", "alt"]]
        return processed_data
    

class LiricalPrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene", "variant"] 
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f)
                header = True
                df = {}
                col_names = []
                with open(file_path) as f:
                    for line in f:
                        if line.startswith("!"): continue
                        line = line.strip().split("\t")
                        if header: 
                            for col in line:
                                col_names.append(col)
                                df[col] = []
                                header = False
                        else:
                            print("here are the colnames")
                            print(col_names)
                            for idx, col in enumerate(col_names):
                                df[col].append(line[idx])
                df = pd.DataFrame.from_dict(df)
                self.patient2gene_results[f] = self.get_result_gene(df)
                print(self.patient2gene_results[f])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)

            self.quant_features_idx[f] = [4,5]
            self.qual_features_idx[f] = [6,7]

    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2variant_results[f] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f)
                header = True
                df = {}
                col_names = []
                with open(file_path) as f:
                    for line in f:
                        if line.startswith("!"): continue
                        line = line.strip().split("\t")
                        if header: 
                            for col in line:
                                col_names.append(col)
                                df[col] = []
                                header = False
                        else:
                            print("here are the colnames")
                            print(col_names)
                            for idx, col in enumerate(col_names):
                                df[col].append(line[idx])
                df = pd.DataFrame.from_dict(df)
                self.patient2variant_results[f] = self.get_result_variant(df)
                print(self.patient2variant_results[f])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2variant_results[f], pf)
            self.quant_features_idx[f] = [8,9,10]
            self.qual_features_idx[f] = [11,12]

    def get_result_gene(self, df):
        processed_data = pd.DataFrame()
        processed_data["rank"] = df["rank"]
        processed_data["score"] = df["compositeLR"]
        for target in ["gene_symbol", "ensembl_id"]:
            self.load_identifier_tanslator("entrez_id",target)
            processed_data[target] = [
                self.identifier_translator[("entrez_id",target)][row.split(":")[1]]
                for row in df["entrezGeneId"]
            ]
        for feature in ["pretestprob", "posttestprob", "diseaseName", "diseaseCurie"]:
            processed_data[feature] = df[feature]
        return processed_data
    
    def get_result_variant(self, df):
        cols_of_interset = ["rank", "score", "varId", "contigName", 
                                             "start", "end", "ref", "alt", 
                                             "pathogenicityScore", "pretestprob", 
                                             "posttestprob", "diseaseName", "diseaseCurie"]
        processed_data ={col: [] for col in cols_of_interset}
        for row in df.to_dict('records'):
            for variant in row["variants"].split(";"):
                variant = variant.strip()
                print(variant)
                # Here I decided to match coding and noncoding variants, but of course that could be no a desire
                # For change: (?P<transcript>NM_[\d\.]+)
                # Instead of: (?P<transcript>[A-Z]{2}_[\d\.]+)
                match = re.match(
                    r"(?P<chr>[^\s:]+):(?P<pos>\d+)(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+|[A-Z]+)?\s+"
                    r"(?P<transcript>[A-Z]{2}_[\d\.]+)(?::(?P<annotation>.*?))?\s+"
                    r"pathogenicity:(?P<pathogenicity>[\d\.]+)\s+\[(?P<genotype>[0-9/]+)\]",
                    variant
                )
                # if not match:
                #     continue  # o print(f"No match: {variant}")
                var_data = match.groupdict()
                processed_data["rank"].append(row["rank"])
                processed_data["score"].append(row["compositeLR"])
                start = int(var_data["pos"])
                ref = var_data["ref"]
                end = start + len(ref) - 1
                var_id = f"{var_data['chr']}-{start}-{end}-{ref}-{var_data['alt']}"
                processed_data["varId"].append(var_id)
                processed_data["contigName"].append(var_data['chr'])
                processed_data["start"].append(start)
                processed_data["end"].append(end)
                processed_data["ref"].append(ref)
                processed_data["alt"].append(var_data["alt"])
                processed_data["pathogenicityScore"].append(float(var_data["pathogenicity"]))
                for feature in ["pretestprob", "posttestprob", "diseaseName", "diseaseCurie"]:
                    processed_data[feature].append(row[feature])

                # quanti : [7,8,9]
                # quali: [10,11,12]
        processed_data = pd.DataFrame.from_dict(processed_data)
        processed_data = processed_data[[*cols_of_interset]]
        return processed_data

class ExomiserPrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene", "variant"] 
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    def post_process_results_genes(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f] = pickle.load(pf)
            else:
                json_results = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                with open(os.path.join(raw_results_dir, f)) as f:
                    json_results = json.load(f)
                self.patient2gene_results[f] = self.get_result_gene(json_results)
                print(self.patient2gene_results[f])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
            self.quant_features_idx[f] = [4,5,6,7]
            self.qual_features_idx[f] = None
    
    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2variant_results[f] = pickle.load(pf)
            else:
                json_results = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                with open(os.path.join(raw_results_dir, f)) as f:
                    json_results = json.load(f)
                self.patient2variant_results[f] = self.get_result_variant(json_results)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2variant_results[f], pf)
            self.quant_features_idx[f] =  [8,9,10,11]
            self.qual_features_idx[f] = None  

    def get_result_variant(self, json_results):
        data = {"gene_symbol": [], "ensembl_id": [], "score": [], "priorityScore": [], "pValue": [], "hiphive_score": [], "omim_score": [], 
                "gene_variant_score": [], "gene_phenotype_score": [], 
                "varId": [], "varIdUniq": [], "pathogenicityScore": [], "contigName": [], "start": [], "end": [], "ref": [], "alt": []}
        for rank, row in enumerate(json_results):
            if not row.get("combinedScore"): continue
            for variants in row["variantEvaluations"]:
                for var_feature in ["pathogenicityScore", "contigName", "start", "end", "ref", "alt"]: # 'phredScore',  "contributingInheritanceModes"
                    data[var_feature].append(variants[var_feature])
                varId = f"{variants['contigName']}:{variants['start']}-{variants['end']}:{variants['ref']}/{variants['alt']}"
                data["varId"].append(varId)
                data["varIdUniq"].append(varId+str(rank))
                # data['clinVarDataInterpretation'].append(variants['pathogenicityData']['clinVarData']['interpretation'])
                # 'pathogenicityData': {'clinVarData': {'primaryInterpretation': 'LIKELY_BENIGN', 'variantEffect': 'MISSENSE_VARIANT'}}
                # data["gene_variant_score"].append(row["geneScores"][2]["variantScore"])
                # data["gene_phenotype_score"].append(row["geneScores"]["phenotypeScore"])
                data["gene_symbol"].append(row["geneIdentifier"]["geneSymbol"])
                data["ensembl_id"].append(row["geneIdentifier"]["geneId"])
                data["score"].append(row["combinedScore"])
                data["pValue"].append(row["pValue"])
                data["hiphive_score"].append(row["priorityResults"]["HIPHIVE_PRIORITY"].get("score"))
                data["omim_score"].append(row["priorityResults"]["OMIM_PRIORITY"]["score"])
        ranking = get_rank_metrics(data["score"], data["varIdUniq"])
        genes = [row[0] for row in ranking]
        print("AHHHHHHHHHHHHHH:",len(genes))
        print("OHHHHHHHHHHH",len(set(data["varIdUniq"])))
        ranking = {row[0]:row[3] for row in ranking}
        data["rank"] = [ranking[gene] for gene in data["varIdUniq"]]

        processed_data = pd.DataFrame()
        for key in ["rank", "score", "varId", "contigName", "start", "end", "ref", "alt", "pathogenicityScore", "pValue", "hiphive_score", "omim_score"]:
            print(key)
            processed_data[key] = data[key]
        return processed_data

    def get_result_gene(self, json_results):
        data = {"gene_symbol": [], "ensembl_id": [], "score": [], "priorityScore": [], "pValue": [], "hiphive_score": [], "omim_score": []}
        for rank, row in enumerate(json_results):
            if not row.get("combinedScore"): continue
            data["gene_symbol"].append(row["geneIdentifier"]["geneSymbol"])
            data["ensembl_id"].append(row["geneIdentifier"]["geneId"])
            data["score"].append(row["combinedScore"])
            data["priorityScore"].append(row["priorityScore"])
            data["pValue"].append(row["pValue"])
            data["hiphive_score"].append(row["priorityResults"]["HIPHIVE_PRIORITY"]["score"])
            data["omim_score"].append(row["priorityResults"]["OMIM_PRIORITY"]["score"])
        ranking = get_rank_metrics(data["score"], data["ensembl_id"])
        genes = [row[0] for row in ranking]
        ranking = {row[0]:row[3] for row in ranking}
        data["rank"] = [ranking[gene] for gene in data["ensembl_id"]]

        processed_data = pd.DataFrame()
        for key in ["rank","score","gene_symbol", "ensembl_id", "priorityScore", "pValue", "hiphive_score", "omim_score"]:
            processed_data[key] = data[key]

        return processed_data




