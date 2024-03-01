# The Block Copolymer Phase Behavior Database

# Description
The Block Copolymer Database (BCDB) is a platform that allows users to search, submit, visualize, benchmark, and download experimental phase measurements and their associated characterization information for di- and multi-block copolymers. To our knowledge, there is no widely accepted data model for publishing experimental and simulation data on block copolymer self-assembly. This proposed data schema with traceable information can accommodate any number of blocks and at the time of publication contains over 5,400 block copolymer total melt phase measurements mined from literature and manually curated and simulation data points of the phase diagram generated from self-consistent field theory that can rapidly be augmented. This database can be accessed on the Community Resource for Innovation in Polymer Technology (CRIPT) web application and the Materials Data Facility. The chemical structure of the polymer is encoded in BigSMILES, an extension of the Simplified Molecular-Input Line-Entry System (SMILES) into the macromolecular domain, and the user can search repeat units and functional groups using SMARTS search syntax (SMILES Arbitrary Target Specification). The user can also query characterization and phase information using the Structured Query Language (SQL) and download custom sets of block copolymer data to train machine learning models. Finally, a protocol is presented in which GPT-4, an AI-powered large language model, can be used to rapidly screen and identify block copolymer papers from literature using only the abstract text and determine whether they have BCDB data, allowing the database to grow as the number of published papers on the world wide web increases. The F1-score for this model is 0.74. This platform is an important step in making polymer data more accessible to the broader community.

# Organization
Here are some files you will encounter in this repository:
- data folder contains an Excel sheet of all block copolymer data divided into tabs like diblock and triblock
- platform folder contains the code for running the data
  - plot.py and tool.py contain the code that powers the user interface for search, visualization, and download
  - All other files are part of the BigSMILES parser, forked from https://github.com/olsenlabmit/BigSMILES_parser. The BigSMILES parser is utilized for molecular search
- environment.yml: packages for installing a new environment on Anaconda to run the software
- tutorial.mp4: brief demonstration of the search, download, and visualization capabilities

# Data Collection Methodology
- Manual curation from tables, text, and figures from peer-reviewed scientific literature. 
- A DOI is provided for each data point. A paper on the full methodology will be published. 
- Data from figures was extracted using WebPlotDigitizer.
  - Citation: Rohatgi, A. (2015). WebPlotDigitizer (Version 4.4) [Computer software]. Retrieved from https://automeris.io/WebPlotDigitizer/index.html.

# Repository Location
http://doi.org/10.5281/zenodo.4780309

# How to Search, Visualize, and Download Data
This platform requires the installation of Anaconda https://www.anaconda.com/products/individual. 
Launch the Anaconda prompt and navigate to the BCDB folder. Locate the environment.yml file. Perform the following steps in the prompt:

1. conda env create -f environment.yml
2. conda activate bcdb
3. cd platform
4. python tool.py

# Creators
- Nathan J. Rebello (https://orcid.org/0000-0002-0178-7701)
- Akash Arora (https://orcid.org/0000-0002-2260-269X)
- Hidenobu Mochigase (hidenobu.mochigase@furukawaelectric.com)
- Tzyy-Shyang Lin (https://orcid.org/0000-0002-8265-6702)
- Jiale Shi (https://orcid.org/0000-0002-5447-3925)
- Debra J. Audus (https://orcid.org/0000-0002-5937-7721)
- Eric S. Muckley (https://orcid.org/0000-0001-7114-5424)
- Ardiana Osmani (https://orcid.org/0000-0002-5532-0285)
- Bradley D. Olsen (https://orcid.org/0000-0002-7272-7140)

# Rights
The visualization software for the database coded in Python is available under the MIT License (https://opensource.org/licenses/MIT) in Zenodo (http://doi.org/10.5281/zenodo.4780309). The dataset is released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) in Zenodo (http://doi.org/10.5281/zenodo.4780309).
