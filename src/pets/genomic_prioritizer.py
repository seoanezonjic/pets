
import os
import pandas as pd
import pickle
import logging
import re
import json
import numpy as np
from xgboost import XGBRanker
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from py_exp_calc.exp_calc import get_rank_metrics

logger = logging.getLogger(__name__)

def parse_hgnc_data():
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
    identifier_map = parse_hgnc_data()
    gene_identifiers = ["gene_symbol", "ensembl_id"]

    def __init__(self):
        self.patient2variant_results = {} # rank    score   chromosome  start   end ref alt
        self.patient2gene_results = {} # rank   score   gene_symbol gene_identifier(ensbl)
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

    def get_combined_results(self, type="gene"):
        if type not in self.priot_type:
            raise ValueError(f"Invalid type: {type}. Must be one of {self.priot_type}")

        common_results, quant_idx_per_patient, qual_idx_per_patient = self.get_common_results(type)
        
        for i, df in enumerate(common_results.values()):
            # add a column with the key
            df["pat_number"] = i

        # Concatenate only common columns
        combined_df = pd.concat([df for df in common_results.values()], ignore_index=True)
        quant_idx = next(iter(quant_idx_per_patient.values()))
        qual_idx = next(iter(qual_idx_per_patient.values()))
        return combined_df, quant_idx, qual_idx
    
    def get_common_results(self, type):
        if type == "gene":
            results_dict = self.patient2gene_results
            first_feat = ["rank", "score", "gene_symbol", "ensembl_id"]
        elif type == "variant":
            results_dict = self.patient2variant_results
            first_feat = ["rank", "score", "varId", "contigName", "start", "end", "ref", "alt"]

        dfs = list(results_dict.values())
        common_cols = set(dfs[0].columns)
        for df in dfs[1:]:
            common_cols &= set(df.columns)
        common_cols = list(common_cols)

        quant_features = list(self.quant_features_idx.values())[0]
        # obtaining the col names
        quant_features = [dfs[0].columns[i] for i in quant_features if dfs[0].columns[i] in common_cols] if quant_features else []
        qual_features = list(self.qual_features_idx.values())[0]
        # obtaining the col names
        qual_features = [dfs[0].columns[i] for i in qual_features if dfs[0].columns[i] in common_cols] if qual_features else []
        common_results = {patient: result[first_feat+quant_features+qual_features] for patient, result in results_dict.items()}
        # Now we obtain the indexes
        final_idx = len(first_feat)
        quant_idx = None
        qual_idx = None
        if quant_features:
            quant_idx = list(range(final_idx, final_idx + len(quant_features)))
            final_idx += len(quant_features)
        if qual_features:
            qual_idx = list(range(final_idx, final_idx + len(qual_features)))
        quant_idx_per_patient = {patient: quant_idx for patient in results_dict.keys()}
        qual_idx_per_patient = {patient: qual_idx for patient in results_dict.keys()}

        return common_results, quant_idx_per_patient, qual_idx_per_patient


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
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f_name] = pickle.load(pf)
            else:
                df = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                self.patient2gene_results[f_name] = df
                # Process the data
                self.patient2gene_results[f_name] = self.get_result_gene(df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f_name], pf)
            self.quant_features_idx[f_name] = None
            self.qual_features_idx[f_name] = [4]
    
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
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f_name] = pickle.load(pf)
            else:
                df = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                self.patient2gene_results[f_name] = df
                # Process the data
                self.patient2gene_results[f_name] = self.get_result_gene(df, rank=2, gene={"gene_symbol": 1, "ensembl_id":0}, score=3, quali_feature=list(range(4,len(df.columns.values))), quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f_name], pf)
            self.qual_features_idx[f_name] = None
            self.quant_features_idx[f_name] = list(range(4,self.patient2gene_results[f_name].shape[1]))

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
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f_name] = pickle.load(pf)
            else:
                df = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                self.patient2gene_results[f_name] = df
                # Process the data
                self.patient2gene_results[f_name] = self.get_result_gene(df, rank=2, gene={"gene_symbol": 1}, score=3)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f_name], pf)
            self.quant_features_idx[f_name] = list(range(4,self.patient2gene_results[f_name].shape[1]-1))
            self.qual_features_idx[f_name] = [self.patient2gene_results[f_name].shape[1]-1]
    
    def _get_hpos_scores(self, row, col=4):
        pattern = r"'\s*([^']*?)\s*'\s*:\s*([0-9.]+)"
        matches = re.findall(pattern, row[col])
        hpos = []
        scores = []
        if not matches:
            return hpos, scores
        else:
            for hpo, score in matches:
                hpo = hpo.strip()
                scores.append(score)
                hpos.append(hpo)
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
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f_name] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f, "prediction", "conf_4Model", "integrated", f"{f}_integrated.csv")
                df = pd.read_csv(file_path, sep=",")
                self.patient2gene_results[f_name] = df
                # Process the data
                self.patient2gene_results[f_name] = self.get_result_gene(df)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f_name], pf)
            self.quant_features_idx[f_name] = None
            self.qual_features_idx[f_name] = None

    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2variant_results[f_name] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f, "prediction", "conf_4Model", "integrated", f"{f}_integrated.csv")
                df = pd.read_csv(file_path, sep=",")
                self.patient2variant_results[f_name] = df
                # Process the data
                self.patient2variant_results[f_name] = self.get_result_variant(df)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2variant_results[f_name], pf)
            self.quant_features_idx[f_name] = None
            self.qual_features_idx[f_name] = None

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
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f_name] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f)
                header = True
                df = {}
                col_names = []
                with open(file_path) as file:
                    for line in file:
                        if line.startswith("!"): continue
                        line = line.strip().split("\t")
                        if header: 
                            for col in line:
                                col_names.append(col)
                                df[col] = []
                                header = False
                        else:
                            for idx, col in enumerate(col_names):
                                df[col].append(line[idx])
                df = pd.DataFrame.from_dict(df)
                self.patient2gene_results[f_name] = self.get_result_gene(df)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f_name], pf)

            self.quant_features_idx[f_name] = [4,5]
            self.qual_features_idx[f_name] = [6,7]

    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        """
        Post-process the gene results from the raw results directory.
        """
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2variant_results[f_name] = pickle.load(pf)
            else:
                file_path = os.path.join(raw_results_dir, f)
                header = True
                df = {}
                col_names = []
                with open(file_path) as file:
                    for line in file:
                        if line.startswith("!"): continue
                        line = line.strip().split("\t")
                        if header: 
                            for col in line:
                                col_names.append(col)
                                df[col] = []
                                header = False
                        else:
                            for idx, col in enumerate(col_names):
                                df[col].append(line[idx])
                df = pd.DataFrame.from_dict(df)
                self.patient2variant_results[f_name] = self.get_result_variant(df)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2variant_results[f_name], pf)
            self.quant_features_idx[f_name] = [8,9,10]
            self.qual_features_idx[f_name] = [11,12]

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
        # pretestprob is in a format 1/900 and i need to apply to convet to float
        processed_data["pretestprob"] = df["pretestprob"].apply(lambda x: float(x.split("/")[0])/float(x.split("/")[1]))
        processed_data["posttestprob"] = df["posttestprob"].apply(lambda x: float(x.replace("%","")))
        for feature in ["diseaseName", "diseaseCurie"]:
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
                # Here I decided to match coding and noncoding variants, but of course that could be no a desire
                # For change: (?P<transcript>NM_[\d\.]+)
                # Instead of: (?P<transcript>[A-Z]{2}_[\d\.]+)
                match = re.match(
                    r"(?P<chr>[^\s:]+):(?P<pos>\d+)(?P<ref>[ACGT]+)>(?P<alt>[A-Z*]+)?\s+"
                    r"(?P<transcript>[A-Z]{2}_[\d\.]+)(?::(?P<annotation>.*?))?\s+"
                    r"pathogenicity:(?P<pathogenicity>[\d\.]+)\s+\[(?P<genotype>[0-9./]+)\]",
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
                processed_data["pretestprob"].append(float(row["pretestprob"].split("/")[0])/float(row["pretestprob"].split("/")[1]))
                processed_data["posttestprob"].append(float(row["posttestprob"].replace("%","")))
                for feature in ["diseaseName", "diseaseCurie"]:
                    processed_data[feature].append(row[feature])

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
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2gene_results[f_name] = pickle.load(pf)
            else:
                json_results = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                with open(os.path.join(raw_results_dir, f)) as file:
                    json_results = json.load(file)
                self.patient2gene_results[f_name] = self.get_result_gene(json_results)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f_name], pf)
            self.quant_features_idx[f_name] = [4,5,6,7]
            self.qual_features_idx[f_name] = None
    
    def post_process_results_variants(self, raw_results_dir, write_tmp=None, read_tmp=False):
        files = os.listdir(raw_results_dir)
        if write_tmp: os.makedirs(write_tmp, exist_ok=True)
        for f in files:
            # Load the data
            f = os.path.basename(f)
            f_name = os.path.splitext(f)[0]
            if read_tmp:
                with open(os.path.join(raw_results_dir, f), 'rb') as pf:
                    self.patient2variant_results[f_name] = pickle.load(pf)
            else:
                json_results = pd.read_csv(os.path.join(raw_results_dir, f), sep="\t")
                with open(os.path.join(raw_results_dir, f)) as file:
                    json_results = json.load(file)
                self.patient2variant_results[f_name] = self.get_result_variant(json_results)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_processed"), 'wb') as pf:
                        pickle.dump(self.patient2variant_results[f_name], pf)
            self.quant_features_idx[f_name] =  [8,9,10,11]
            self.qual_features_idx[f_name] = None  

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
                data["omim_score"].append(row["priorityResults"]["OMIM_PRIORITY"].get("score"))
        ranking = get_rank_metrics(data["score"], data["varIdUniq"])
        ranking = {row[0]:row[3] for row in ranking}
        data["rank"] = [ranking[gene] for gene in data["varIdUniq"]]

        processed_data = pd.DataFrame()
        for key in ["rank", "score", "varId", "contigName", "start", "end", "ref", "alt", "pathogenicityScore", "pValue", "hiphive_score", "omim_score"]:
            processed_data[key] = data[key]
        return processed_data

    def get_result_gene(self, json_results):
        data = {"gene_symbol": [], "ensembl_id": [], "score": [], "priorityScore": [], "pValue": [], "hiphive_score": [], "omim_score": []}
        for rank, row in enumerate(json_results):
            if not row.get("combinedScore"): continue
            data["gene_symbol"].append(row["geneIdentifier"]["geneSymbol"])
            data["ensembl_id"].append(row["geneIdentifier"]["geneId"])
            data["score"].append(row["combinedScore"])
            data["priorityScore"].append(row.get("priorityScore",None))
            data["pValue"].append(row["pValue"])
            data["hiphive_score"].append(row["priorityResults"]["HIPHIVE_PRIORITY"].get("score"))
            data["omim_score"].append(row["priorityResults"]["OMIM_PRIORITY"].get("score"))
        ranking = get_rank_metrics(data["score"], data["ensembl_id"])
        genes = [row[0] for row in ranking]
        ranking = {row[0]:row[3] for row in ranking}
        data["rank"] = [ranking[gene] for gene in data["ensembl_id"]]

        processed_data = pd.DataFrame()
        for key in ["rank","score","gene_symbol", "ensembl_id", "priorityScore", "pValue", "hiphive_score", "omim_score"]:
            processed_data[key] = data[key]

        return processed_data

class MetaGenomicPrioritizer:
    def __init__(self, prioritizers):
        self.prioritizers = prioritizers # prioritizer_name -> prioritizer_instance
        # features
        self.feature_gene = {}
        self.feature_variant = {}
        self.feature_quant_idx = {}
        self.feature_qual_idx = {}
        self.feature_columns = []  # columns used as features
        # label
        self.label = {} # Patient -> Positives ids
        # corpus split
        self.train_patients = None
        self.test_patients = None
        # model
        self.model = None
        # Results
        self.patient2variant_results = {} # rank    score  chromosome   start end ref alt
        self.patient2gene_results = {} # rank   score gene_symbol gene_identifier(ensbl)
        self.quant_features_idx = {}
        self.qual_features_idx = {}

    # Feature extraction
    ####################

    def get_features(self, type="gene", dropna=False):
        # Mapping type to attributes
        id_candidates = {"gene": "gene_symbol", "variant": "varId"}
        results_attr = {"gene": "patient2gene_results", "variant": "patient2variant_results"}
        merged_attr = {"gene": "feature_gene", "variant": "feature_variant"}
        features_to_remove = {"gene": ["ensembl_id"], "variant": ["contigName", "start", "end", "ref", "alt"]}

        if type not in id_candidates:
            raise ValueError(f"Invalid type: {type}")

        id_candidate = id_candidates[type]
        results_key = results_attr[type]
        merged_key = merged_attr[type]
        feature_to_remove = features_to_remove[type]

        self.clean_features(results_key, type)
        all_patients = list(getattr(self.prioritizers[list(self.prioritizers.keys())[0]], results_key).keys())

        # Merge per patient
        merged_results = {}
        for patient in all_patients:
            dfs = []
            quantitative_idx = []
            qualitative_idx = []
            n_total_cols = 0
            for name, prioritizer in self.prioritizers.items():
                df = getattr(prioritizer, results_key).get(patient)
                name = name[0]
                if df is None or df.empty:
                    continue

                # Keep only one column for ID + rename the rest with the prioritizer name
                renamed = df.copy()
                cols_to_rename = [col for col in df.columns if col != id_candidate and col not in feature_to_remove]
                renamed = renamed[[id_candidate] + cols_to_rename]
                renamed.columns = [id_candidate] + [f"{col}_{name}" for col in cols_to_rename]
                quantitative_idx.extend([n_total_cols + 1, n_total_cols + 2])
                if prioritizer.quant_features_idx.get(patient): 
                    quantitative_idx.extend([ idx + n_total_cols -len(feature_to_remove) for idx in prioritizer.quant_features_idx[patient]])
                if prioritizer.qual_features_idx.get(patient):
                    qualitative_idx.extend([val + n_total_cols - len(feature_to_remove) for val in prioritizer.qual_features_idx[patient]])
                n_total_cols += len(cols_to_rename)
                dfs.append(renamed)

            if not dfs:
                merged_results[patient] = pd.DataFrame()
                continue

            # Merge all DataFrames on the id_candidate
            merged = dfs[0]
            
            for df in dfs[1:]:
                merged = pd.merge(merged, df, on=id_candidate, how="outer")
            
            if dropna:
                merged_results[patient] = merged.dropna() # merged
            else:
                merged_results[patient] = merged

            self.feature_qual_idx[patient] = qualitative_idx
            self.feature_quant_idx[patient] = quantitative_idx


        # assign to attribute
        setattr(self, merged_key, merged_results)
    
    def clean_features(self, results_key, type):
        # Obtain results with common columns inside each prioritizer
        for prioritizer in self.prioritizers.values():
            prio_table, quantitative_feature, qualitative_feature = prioritizer.get_common_results(type)
            prioritizer.quant_features_idx = quantitative_feature
            prioritizer.qual_features_idx = qualitative_feature
            setattr(prioritizer, results_key, prio_table)

        # Select just patients with results in all prioritizers.
        prioritizer_names = list(self.prioritizers.keys())
        all_patients = set(getattr(self.prioritizers[prioritizer_names[0]], results_key).keys())
        for prioritizer_name in prioritizer_names[1:]:
            all_patients_for_priotizer = set(getattr(self.prioritizers[prioritizer_name], results_key).keys())
            all_patients = all_patients.intersection(all_patients_for_priotizer)
        all_patients = sorted(list(all_patients))

        # Update with just the common patients
        for prioritizer in self.prioritizers.values():
            clean_features = {}
            results_dict = getattr(prioritizer, results_key)
            for patient in all_patients:
                clean_features[patient] = results_dict[patient]
            setattr(prioritizer, results_key, clean_features)

    # labels
    ##################

    def load_patient_labels(self, patient_labels):
        for patient_label in patient_labels:
            patient = patient_label[0]
            candidateID = patient_label[1]
            if not self.label.get(patient):
                self.label[patient] = [candidateID]
            else:
                self.label[patient].append(candidateID)

    # Split and preparing corpus
    ######################

    def get_all_patients(self, type="gene"):
        merged_key = "feature_gene" if type == "gene" else "feature_variant"
        all_patients = list(getattr(self, merged_key).keys())
        return all_patients

    def split_patients(self, type="gene", test_size=0.3, random_state=42):
        all_patients = self.get_all_patients(type)
        self.train_patients, self.test_patients = train_test_split(
            all_patients, test_size=1, random_state=random_state
        )

    # Training 
    ######################

    def prepare_training_data(self, type="gene", label_col_substring="score"):
        merged_key = "feature_gene" if type == "gene" else "feature_variant"
        merged_results = getattr(self, merged_key)

        train_dfs = [merged_results[p] for p in self.train_patients if not merged_results[p].empty]
        df = pd.concat(train_dfs, keys=self.train_patients, names=["patient_id"])
        df = df.reset_index()

        # identifier by type
        id_col = "gene_symbol" if type == "gene" else "varId"

        # numerical features
        #feature_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        #feature_cols = [col for col in feature_cols if label_col_substring not in col]
        feature_cols = df.iloc[:, self.quant_features_idx.values()[0]]#.select_dtypes(include=[np.number]).columns.tolist()
        #feature_cols = [col for col in feature_cols if label_col_substring not in col]
        self.feature_columns = feature_cols[1:] # TODo check this

        X = df[self.feature_columns]

        # Obtain labels
        y = []
        for row in zip(df["patient_id"], df[id_col]):
            if row[1] in self.label[row[0]]:
                y.append(1)
            else:
                y.append(0)

        groups = df.groupby("patient_id").size().to_numpy()
        return X, y, groups

    def train_model(self, type="gene"):
        X, y, groups = self.prepare_training_data(type="gene")
        self.model.train(X,y, groups)

    # Prediction
    ######################

    ## Predict the test set
    def predict_test(self, type="gene"):
        merged_key = "feature_gene" if type == "gene" else "feature_variant"
        id_col = "gene_symbol" if type == "gene" else "varId"
        merged_results = getattr(self, merged_key)
        predict_results = getattr(self, f"patient2{type}_results")

        if self.test_patients is None:
            self.test_patients = list(merged_results.keys())

        for patient in self.test_patients:
            df = merged_results[patient]
            if not self.feature_columns:
                feature_cols = df.select_dtypes(include=[np.number]).columns.tolist()
                feature_cols = [col for col in feature_cols if "score" not in col]
                self.feature_columns = feature_cols #[1:] # TODo check this Aqui esta pasando algo incoherente con respecto al training
            X_test = df[self.feature_columns]
            scores = self.model.predict(X_test)

            df = df.copy()
            df["score"] = scores
            df = df.sort_values("score", ascending=False)
            # adding ranking pos
            ranking = get_rank_metrics(df["score"].tolist(), df[id_col].tolist())
            ranking = {row[0]:row[3] for row in ranking}
            df["rank"] = [ranking[gene] for gene in df[id_col]]
            cols = ["rank", "score"] + [col for col in df.columns if col not in ["rank", "score"]]
            predict_results[patient] = df[cols]
            self.quant_features_idx[patient] = [val + 2 for val in self.feature_quant_idx[patient]]
            self.qual_features_idx[patient] = [val + 2 for val in self.feature_qual_idx[patient]]

    
    def get_combined_results(self, type="gene"):
        common_results = getattr(self, f"patient2{type}_results")
        for i, df in enumerate(common_results.values()):
            # add a column with the key
            df["pat_number"] = i

        # Concatenate only common columns
        combined_df = pd.concat([df for df in common_results.values()], ignore_index=True)
        quant_idx = next(iter(self.quant_features_idx.values()))
        qual_idx = next(iter(self.qual_features_idx.values()))
        return combined_df, quant_idx, qual_idx

class HeuristicModel():
    """
    A heuristic model that ranks items based on the minimum rank
    """

    def train(self, X, y=None, groups=None):
        pass  # no training needed

    def predict(self, X):
        rank_cols = [col for col in X.columns if col.startswith("rank_")]
        X[rank_cols] = X[rank_cols].apply(pd.to_numeric, errors='coerce')
        return -1 * X[rank_cols].min(axis=1, skipna=True).to_numpy()

class XGBoostRankerModel:

    def __init__(self, **params):
        default_params = {
            "objective": "rank:pairwise", 
            "learning_rate": 0.1,
            "n_estimators": 100,
            "max_depth": 6,
            "random_state": 42
        }
        default_params.update(params)
        self.model = XGBRanker(**default_params)

    def train(self, X, y, group):
        self.model.fit(X, y, group=group)

    def predict(self, X):
        return self.model.predict(X) #Predicted scores (higher means more relevant)

class LogisticRegressionModel:

    def __init__(self, **params):
        default_params = {}
        default_params.update(params)
        self.model = LogisticRegression()

    def train(self, X, y, group):
        self.model.fit(X, y)

    def predict(self,X):
        return self.model.predict_proba(X)[:,1]