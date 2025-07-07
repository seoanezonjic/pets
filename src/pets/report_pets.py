from py_report_html import Py_report_html
from importlib.resources import files
import numpy as np
from collections import defaultdict
import py_exp_calc.exp_calc as pxc
import networkx as nx
import pandas as pd

########################################
## Monkey Patching Methods
########################################

Py_report_html.additional_templates.append(str(files('py_semtools').joinpath('templates')))

## UTILS FUNCTIONS
def _transform_value(self, value, method):
    if method == "cubic_root":
        return value**(1/3)
    elif method == "bins":
        return self._get_alpha_bin(value)
    elif method == "none":
        return value
    elif callable(method):
        return method(value)
    
def make_title(self, type, id, sentence):
        if type == "table":
                key = f"tab:{id}"
                html_title = f"<p style='text-align:center;'> <b> {type.capitalize()} {self.add_table(key)} </b> {sentence} </p>"
        elif type == "figure":
                key = id
                html_title = f"<p style='text-align:center;'> <b> {type.capitalize()} {self.add_figure(key)} </b> {sentence} </p>"
        return html_title

def get_top_rank(self, table, top):
    top_genes = table[table["rank"] <= top]
    return top_genes

def get_corr_table(self, table, columns, method="spearman"):
    df = table.iloc[:,columns]
    corr = df.corr(method=method)
    corr = corr.dropna(axis=1, how='all')
    corr = corr.dropna(axis=0, how='all')
    corr_table = corr.values.tolist()
    corr_table = [[corr.index[i]] + row for i, row in enumerate(corr_table)]
    corr_table.insert(0, [" "] + corr.columns.tolist())
    return corr_table

def df_to_numeric(self, table, numeric_cols):
    table = table.copy()
    for idx in numeric_cols:
        col_name = table.columns[idx]
        table[col_name] = pd.to_numeric(table[col_name], errors='coerce')
    return table

def df_to_list(self, df):
    l = df.values.tolist()
    l.insert(0, df.columns.tolist())
    return l
     

#### METHODS FOR SIMILARITY MATRIX HEATMAP

#### LOADING ALL MONKEYPATCHED METHODS
Py_report_html._transform_value = _transform_value
Py_report_html.get_top_rank = get_top_rank
Py_report_html.get_corr_table = get_corr_table
Py_report_html.make_title = make_title
Py_report_html.df_to_numeric = df_to_numeric
Py_report_html.df_to_list = df_to_list