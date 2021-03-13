# BCDB
The Block Copolymer Database (BCDB) is a database platform that allows a user to search, submit, visualize, and download experimental phase measurements and their associated characterization information for di- and multi-block copolymers. This database template can accommodate any number of blocks and at the time of publication contains over 5,300 diblock copolymer melt phase measurements mined from literature and manually collated. The chemical structure of the polymer is encoded in BigSMILES, an extension of the Simplified Molecular-Input Line-Entry System (SMILES) into the macromolecular domain, and the user can search repeat units and functional groups using the search syntax SMARTS (SMILES Arbitrary Target Specification). The user can also query characterization and phase information using the Structured Query Language (SQL). This collection of data facilitates the development of data-driven models, benchmarking and visualization of new data against existing samples, and compilation of custom sets of block copolymer data from the literature.  

# Installation
This platform requires the installation of Anaconda: https://docs.anaconda.com/anaconda/install/.
Launch the Anaconda prompt and navigate to the download BCDB folder. 
Locate the environment.yml file. 

1. conda env create -f environment.yml
2. conda activate bcdb
3. cd platform
4. python tool.py