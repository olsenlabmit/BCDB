# Block Copolymer Phase Behavior Database (BCDB)
# Creators (ORCID)
- Nathan J. Rebello (https://orcid.org/0000-0002-0178-7701)
- Akash Arora (https://orcid.org/0000-0002-2260-269X)
- Hidenobu Mochigase
- Tzyy-Shyang Lin
- Debra J. Audus (https://orcid.org/0000-0002-5937-7721)
- Bradley D. Olsen (https://orcid.org/0000-0002-7272-7140)
# Dates
Project started on March 2020. Data released on March 2021. Data extracted from peer-reviewed literature from the 1980s to 2010s. 
# Subject 
The Block Copolymer Database (BCDB) is a database platform that allows a user to search, submit, visualize, and download experimental phase measurements and their associated characterization information for di- and multi-block copolymers. This database template can accommodate any number of blocks and at the time of publication contains over 5,300 block copolymer melt phase measurements mined from literature and manually curated. The chemical structure of the polymer is encoded in BigSMILES, an extension of the Simplified Molecular-Input Line-Entry System (SMILES) into the macromolecular domain, and the user can search repeat units and functional groups using the search syntax SMARTS (SMILES Arbitrary Target Specification). The user can also query characterization and phase information using the Structured Query Language (SQL). This platform is an important step in making polymer data more accessible to the broader community, facilitates benchmarking and visualization of new data against existing samples, enables to user to search and download custom sets of block copolymer data from the literature, and will drive the development of data-driven models. We note a model trained on this data significantly outperforms self-consistent field theory. 
# Organization
Here are some files you will encounter in this repository:
- Tutorial.mp4: brief demonstration of the search, download, and visualization capabilities.
environment.yml: packages for installing a new environment on Anaconda to run the software.
- names_dict.csv: block chemistry dictionary
- phases_dict.csv: phase dictionary
- data folder contains all data in separate CSV files so that the data is open and accessible for each n-block chemistry. 
  - column_meanings.csv: column headers with their meanings. The only difference between n-block columns is the Individual Block category. For example, the diblock has two sets of columns from this category describing each block, a triblock has three, and a tetrablock has four.
- platform folder contains the code for running the data (see Run the Tool to run the platform). 
  - plot.py and tool.py contain the code for the user interface for search, visualization, and download
  - All other files are part of the BigSMILES parser, forked from https://github.com/olsenlabmit/BigSMILES_parser. The BigSMILES parser is utilized for stochastic graph search.
# Funders 
This work is funded by the Community Resource for Innovation in Polymer Technology, a project supported by the National Science Foundation (NSF) Convergence Accelerator program (NSF Convergence Accelerator Research-2040636). 
# Rights
The visualization software for the database coded in Python is available under the MIT License (https://opensource.org/licenses/MIT) in Zenodo (http://doi.org/10.5281/zenodo.4780309). The dataset is released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) in Zenodo (http://doi.org/10.5281/zenodo.4780309).
# Location
http://doi.org/10.5281/zenodo.4780309
# Methodology
- Manual curation from tables, text, and figures from peer-reviewed scientific literature. 
- A DOI is provided for each data point. A paper on the full methodology will be published within a week or two. 
- Data from figures was extracted using WebPlotDigitizer.
  - Citation: Rohatgi, A. (2015). WebPlotDigitizer (Version 4.4) [Computer software]. Retrieved from https://automeris.io/WebPlotDigitizer/index.html)
# Requirements
This platform requires the installation of Anaconda https://www.anaconda.com/products/individual.
# How to Use
Launch the Anaconda prompt and navigate to the BCDB folder. Locate the environment.yml file. Perform the following steps in the prompt:

1. conda env create -f environment.yml
2. conda activate bcdb
3. cd platform
4. python tool.py
