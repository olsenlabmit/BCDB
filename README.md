# BCDB
The Block Copolymer Database (BCDB) is a database platform that allows a user to search, submit, visualize, and download experimental phase measurements and their associated characterization information for di- and multi-block copolymers. This database template can accommodate any number of blocks and at the time of publication contains over 5,300 diblock copolymer melt phase measurements mined from literature and manually collated. The chemical structure of the polymer is encoded in BigSMILES, an extension of the Simplified Molecular-Input Line-Entry System (SMILES) into the macromolecular domain, and the user can search repeat units and functional groups using the search syntax SMARTS (SMILES Arbitrary Target Specification). The user can also query characterization and phase information using the Structured Query Language (SQL). This collection of data facilitates the development of data-driven models, benchmarking and visualization of new data against existing samples, and compilation of custom sets of block copolymer data from the literature.  

# Licenses
The visualization software for the database coded in Python is available under the MIT License (https://opensource.org/licenses/MIT). The dataset is released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).

# Requirements
This platform requires the installation of Anaconda https://www.anaconda.com/products/individual.

# Installation
Launch the Anaconda prompt and navigate to the BCDB folder. Locate the environment.yml file. Perform the following steps in the Anaconda prompt:

1. conda env create -f environment.yml
2. conda activate bcdb
3. cd platform
4. python tool.py

Refer to the video Tutorial.mp4 for a brief demonstration of the search, download, and visualization capabilities.
