# NLCD : A method to discover nonlinear causal relations among genes

This repo contains the scripts and data used in the paper: "NLCD: A method to discover nonlinear causal relations among genes. Aravind Easwar and Manikandan Narayanan". 

NLCD stands for Non-Linear Causal Discovery. It is a method to discover causal relations between two genes or between two traits in general (including gene expression or other clinical traits) from matched trait and genetic data. Please find more information about the method and its evaluation and application in our paper.

## Repo Folder Overview

-> `data`:  Contains simulation data and code to generate the simulation data \
-> `figures`: Contains figures and code to reproduce the plots, tables \
-> `results`: \
-> `scripts`: \
-> `yeast`: Contains yeast datasets, and codes to process the raw yeast data. 

## Installation
```
git clone https://github.com/BIRDSgroup/NLCD.git
cd NLCD
pip install -r requirements.txt
```
## Quick Usage (aka Getting Started)

To test NLCD on a single triplet (L, A, B), where L is the genetic instrument variable (SNP) and A and B are the gene pairs (or traits in general), you can run the following commands:  
* Open the user_run.ipynb file, and call the nlcd_single() function as follows: 
* `nlcd_single(L,A,B,shuffles,algo,sample_seed=None,verbose=True,reverse=False)`
  * In the above call, shuffle is the number of times to do the permutation testing, and algo can be SVR, KRR and ANN.
  * The argument seed can also be set to enable reproducible results (if no sample seed is provided, the program will generate a sample seed). 
  * Above function tests if A is causal for B ($A \rightarrow B$), but to ensure B is not causal for A, then we can call nlcd_single in the reverse direction like so:
  * `nlcd_single(L,A,B,shuffles,algo,sample_seed=None,verbose=True,reverse=*True*)`
* Above function is to test a single trio. If you have a list of trios, you can also run NLCD on each trio in the batch mode like so:  
  * `nlcd_batch(data,shuffles,algo,reverse=False,sample_seed=None)`
  * The data variable contains the set of trios to be tested.
  * The sample_seed here refers to a parent seed -- for each triplet, a separate seed (child seed) is produced and output (so that they can be used later with `nlcd_single()` function if needed for reproducibility).     

## Reproduction of Simulated data, Results and Figures
Instructions are provided in the data folder,results, figures and yeast folder on how to reproduce the results,data and figures.    
Files/Results that are large in size can be found here : https://drive.google.com/drive/folders/1W2uq8H3PmJKspsceDLFLkwFMoIW63jgf?usp=drive_link

## Contact 
If you are facing any issues/queries, please contact Aravind Easwar (aeaswar81@gmail.com), Dr. Manikandan Narayanan (nmanik@cse.iitm.ac.in)

## License preamble 

Copyright 2024 BIRDS Group, IIT Madras

NLCD is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

NLCD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with NLCD.  If not, see <https://www.gnu.org/licenses/>.


