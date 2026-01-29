import os
from py_semtools.ontology import Ontology

file_Hierarchical = {"file": os.path.join("./hp.obo"), "name": "hp"}
hierarchical = Ontology(file= file_Hierarchical["file"],load_file= True)
hierarchical.precompute()
hierarchical.write(os.path.join("./hp.json"))
