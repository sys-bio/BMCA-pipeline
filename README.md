# BMCA-pipeline

This repository accompanies [paper citation](). It contains model files, data files, and Jupyter notebooks designed to run Bayesian Metabolic Control Analysis (BMCA) and perform subsequent analyses.

## Installation

To set up the environment:

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).
2. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/BMCA-pipeline.git
   cd BMCA-pipeline
3. Create a conda environment and install the dependencies listed in `env.yml`

## Usage
These notebooks can be used to reproduce the figures in the paper. 
Each Jupyter notebook is set up with specific models and data files. 
Some notebooks omit certain types of data. 
All notebooks contain analysis following the BMCA run. 
Analysis regarding aggregate data are found in notebooks containing "-stats" in their title.


## License

This code is distributed under the MIT License.

## Acknowledgements

The `emll` folder was forked from [https://github.com/pnnl-predictive-phenomics/emll](https://github.com/pnnl-predictive-phenomics/emll). A copy of it is here to preserve the correct dependencies. 

This project is developed by the [UW Sauro Lab](https://www.sys-bio.org) at the Department of Bioengineering, University of Washington, Seattle.

The information, data, or work presented herein was funded in part by the The Bioenergy Technologies Office (BETO) within the U.S. Department of Energy’s (DOE) Office of Energy Efficiency and Renewable Energy under Award Number DE-EE0008927 and the U.S. Department of Energy, Office of Science, Biological and Environmental Research Program, under Award Number DE-SC0023091. This report was prepared as an account of work sponsored by an agency of the United States Government. Neither the United States Government nor any agency thereof, nor any of their employees, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately owned rights. Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof. The views and opinions of the authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

HMS was supported by NIH Biomedical Imaging and Bioengineering award P41 EB023912 through the Center for Reproducible Biomedical Modeling.

