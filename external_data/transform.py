from py_semtools import Ontology, JsonParser
import os

file_Hierarchical = {"file": os.path.join("./hp.obo"), "name": "hp"}
hierarchical = Ontology(file= file_Hierarchical["file"],load_file= True)
hierarchical.precompute()
hierarchical.write(os.path.join("./hp.json"))
