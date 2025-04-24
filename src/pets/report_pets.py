from py_report_html import Py_report_html
from importlib.resources import files
import numpy as np
from collections import defaultdict
from py_semtools.cons import Cons
import py_exp_calc.exp_calc as pxc
import networkx as nx

########################################
## Monkey Patching Methods
########################################

Py_report_html.additional_templates.append(str(files(Cons.TEMPLATES).joinpath('')))



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


#### METHODS FOR SIMILARITY MATRIX HEATMAP

#### LOADING ALL MONKEYPATCHED METHODS
Py_report_html._transform_value = _transform_value