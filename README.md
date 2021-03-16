# Block Copolymer Phase Behavior Database (BCDB)
# Creators: 
Nathan J. Rebello, Akash Arora, Hidenobu M. Mochigase, Tzyy-Shyang Lin, Bradley D. Olsen
# Dates
Project started in January 2020. Data released in March 2021. Data extracted from literature articles from the 1980s - 2010s. 
# Subject 
The Block Copolymer Database (BCDB) is a database platform that allows a user to search, submit, visualize, and download experimental phase measurements and their associated characterization information for di- and multi-block copolymers. This database template can accommodate any number of blocks and at the time of publication contains over 5,300 block copolymer melt phase measurements mined from literature and manually collated. The chemical structure of the polymer is encoded in BigSMILES, an extension of the Simplified Molecular-Input Line-Entry System (SMILES) into the macromolecular domain, and the user can search repeat units and functional groups using the search syntax SMARTS (SMILES Arbitrary Target Specification). The user can also query characterization and phase information using the Structured Query Language (SQL). This collection of data facilitates the development of data-driven models, benchmarking and visualization of new data against existing samples, and compilation of custom sets of block copolymer data from the literature. 
# Contents
Before running the platform, the user can watch Tutorial.mp4. The environment.yml is for installing a new environment with the appropriate packages to run the platform. The nblocks_csv folder contains all of the data in separate CSV files (UTF-8) so that the data is open and accessible. The platform folder contains the code for running the data (see Run the Tool to run the platform).
# Funders 
This work is funded by the Community Resource for Innovation in Polymer Technology, a project supported by the National Science Foundation (NSF) Convergence Accelerator program (NSF Convergence Accelerator Research-2040636). 
# Rights
The visualization software for the database coded in Python is available under the MIT License (https://opensource.org/licenses/MIT) in Zenodo (http://doi.org/10.5281/zenodo.XXXXXXX). The dataset is released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) in Zenodo (http://doi.org/10.5281/zenodo.XXXXXXX).
# Languages
Python
# Location
http://doi.org/10.5281/zenodo.XXXXXXX
# Methodology
Manual curation from tables, text, and figures from peer-reviewed scientific literature. Data from figures was extracted using WebPlotDigitizer (Citation: Rohatgi, A. (2015). WebPlotDigitizer (Version 4.4) [Computer software]. Retrieved from https://automeris.io/WebPlotDigitizer/index.html).
# Requirements
This platform requires the installation of Anaconda https://www.anaconda.com/products/individual.
# Run the Tool
Launch the Anaconda prompt and navigate to the BCDB folder. Locate the environment.yml file. Perform the following steps in the Anaconda prompt:

1. conda env create -f environment.yml
2. conda activate bcdb
3. cd platform
4. python tool.py

Refer to the video Tutorial.mp4 for a brief demonstration of the search, download, and visualization capabilities.
