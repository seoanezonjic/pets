.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/pets.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/pets
    .. image:: https://readthedocs.org/projects/pets/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://pets.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/pets/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/pets
    .. image:: https://img.shields.io/pypi/v/pets.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/pets/
    .. image:: https://img.shields.io/conda/vn/conda-forge/pets.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/pets
    .. image:: https://pepy.tech/badge/pets/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/pets
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/pets

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

====
pets
====


    Patient Exploration Tools Suite (PETS)


The Patient Exploration Tools Suite (PETS) is a suite of tools programmed in Python for the analysis of phenotypic and genotypic data from patients with genetic diseases. It is structured into three different functionalities:

* Cohort Analyzer: a tool for evaluating the quality of data in patient cohorts or disease groups. It uses phenotypic information encoded with the Human Phenotype Ontology and calculates semantic similarity of phenotypes with configurable methods (such as Lin and Resnik). Patients or diseases are then clustered hierarchically. Results are presented in HTML files for easy interpretation.
* Evidence Profiler: a tool for prioritizing genomic variants associated with those characterized in patient cohorts. It leverages data from sources like the MONARCH Initiative, which includes information on genetic diseases, associated phenotypes, and characterized variants to perform this prioritization.
* Implementation of PhenoPackets: for the standardization of phenotypic and genotypic information of the patients or diseases analyzed, enhancing the accuracy and interoperability of variant prioritization.

PETS also includes auxiliary scripts for file manipulation, so that the input information is transformed into a single type with which the rest of the tools will work. The library also includes example files to work on.

Please, cite us as: Rojano E., Cordoba-Caballero J., Moreno-Jabato F., Gallego D., Serrano M., Perez B., Pares-Aguilar A., Perkins JR., Ranea JAG., Seoane-Zonjic P. Evaluating, Filtering and Clustering Genetic Disease Cohorts Based on Human Phenotype Ontology Data with Cohort Analyzer. J. Pers. Med. 2021, 11(8), 730; https://doi.org/10.3390/jpm11080730.