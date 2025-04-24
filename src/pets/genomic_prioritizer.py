
import os
import pandas as pd
import pickle
import logging
import re

logger = logging.getLogger(__name__)

def parse_hgnc_data() -> pd.DataFrame:
    base_path = os.path.dirname(__file__) 
    hgnc_path = os.path.join(base_path, "external_data", "hgnc_complete_set.txt")

    df = pd.read_csv(hgnc_path, sep="\t", dtype={"omim_id": str})

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


# def create_gene_identifier_map() -> pd.DataFrame:
#     """
#     Crea un DataFrame con el mapeo entre identificadores de genes y símbolos, usando los datos HGNC.
#     """
#     logger.info("Creando mapeo de identificadores de genes.")
#     df = parse_hgnc_data()

#     # 'Unpivot' las columnas de ID (melt)
#     melted = df.melt(
#         id_vars=["gene_symbol", "prev_symbols"],
#         value_vars=["ensembl_id", "hgnc_id", "entrez_id", "refseq_accession"],
#         var_name="identifier_type",
#         value_name="identifier"
#     )

#     # Añadir prefijos
#     prefix_map = {
#         "ensembl_id": "ensembl:",
#         "hgnc_id": "",
#         "entrez_id": "ncbigene:",
#         "refseq_accession": ""
#     }
#     melted["prefix"] = melted["identifier_type"].map(prefix_map).fillna("")

#     return melted


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

    def load_identifier_tanslator(self, source_id, target_id):
        self.identifier_translator[(source_id, target_id)] = {}
        for _, row in type(self).identifier_map[[source_id, target_id]].iterrows():
            self.identifier_translator[(source_id, target_id)][row[source_id]] = row[target_id]
    
    def post_process(self, raw_results_dir, write_tmp=None, read_tmp=False):
        # Placeholder for the implementation
        # This function will read the raw results, process them, and save the standardized output.
        if "gene" in self.priot_type:
            # Process gene results
            self.post_process_results_genes(raw_results_dir, write_tmp, read_tmp)
        if "variant" in self.priot_type:
            # Process variant results
            self.post_process_results_variants(raw_results_dir, write_tmp, read_tmp)

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
                self.patient2gene_results[f] = self.process_data(df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
    
    def process_data(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[]):
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
                self.patient2gene_results[f] = self.process_data(df, rank=2, gene={"gene_symbol": 1, "ensembl_id":0}, score=3, quali_feature=list(range(4,len(df.columns.values))), quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
    
    def process_data(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[]):
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
                self.patient2gene_results[f] = self.process_data(df, rank=2, gene={"gene_symbol": 1}, score=3)
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
    
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

    def process_data(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[], quant_feature=[]):
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

class ExomiserPrioritizer(GenomicPrioritizer):

    def __init__(self):
        super().__init__() 
        self.priot_type = ["gene", "variant"] 

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
                self.patient2gene_results[f] = self.process_data(df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[])
                # Save the processed data
                if write_tmp:
                    with open(os.path.join(write_tmp, f+"_proccesed"), 'wb') as pf:
                        pickle.dump(self.patient2gene_results[f], pf)
    
    def process_data(self, df, rank=0, gene={"gene_symbol": 1}, score=3, quali_feature=[4], quant_feature=[]):
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




