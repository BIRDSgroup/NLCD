# NLCD : A method to discover nonlinear causal relations among genes

Non Linear Causal Discovery (NLCD) 

This repo contains the scripts and data used in the paper : cite our paper.

## Repo Folder Overview

-> `data`:  Contains simulation data and code to generate the simulation data \
-> `figures`: Contains figures and code to reproduce the plots, tables \
-> `results`: \
-> `scripts`: \
-> `yeast`: Contains yeast datasets, and code to process the yeast datasets. 

## Installation
```
git clone https://github.com/BIRDSgroup/NLCD.git
cd NLCD
pip install -r requirements.txt
```
## Usage 

Open the user_run.ipynb file, and based on the need call nlcd_single() function or nlcd_batch()  
L is the instrument variable, A and B are the gene pairs, shuffle is the number of times to do the permutation testing, and algorithm can be SVR, KRR and ANN, seed can also be set to reproduce the results  
`nlcd_single(L,A,B,shuffles,algo,sample_seed=None,verbose=True,reverse=False)`: This function is to test a single trio, reverse parameter checks the relation from $B \rightarrow A$, the default is testing from $A \rightarrow B$. If no sample seed is provided, the program will generate a sample seed.  
`nlcd_batch(data,shuffles,algo,reverse=False,sample_seed=None)`: This function can test for a set of trios. data variable contains the set of trios to be tested. This function has a parent seed which we can pass as sample_seed, for each of the trios separate seed (child seed) is produced so that their values can be reproduced using `nlcd_single()` function.     

## Reproduction of Simulated data, Results and Figures
Instructions are provided in the data folder,results, figures and yeast folder on how to reproduce the results,data and figures.    
Files/Results that are large in size can be found here : https://drive.google.com/drive/folders/1W2uq8H3PmJKspsceDLFLkwFMoIW63jgf?usp=drive_link

## Contact 
If you are facing any issues, please contact Aravind Easwar (aeaswar81@gmail.com), Dr. Manikandan Narayanan (nmanik@cse.iitm.ac.in)

## License preamble 

Copyright 2024 BIRDS Group, IIT Madras

NLCD is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

NLCD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with NLCD.  If not, see <https://www.gnu.org/licenses/>.


