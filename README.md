# NLCD : A method to discover nonlinear causal relations among genes
The repo contains the scripts and data used in the paper : cite our paper.  
The scripts folder contains :  
-> `nlcd_main.py`: This file has all the functions defined  
-> `nlcd_simulation.py` : This file is used to run the simulation data and yeast data results  
-> `nlcd_user.py`: This file imports the main file and adds two more functions for the user , add seeds  
-> `user_run.ipynb`: This notebook imports the above nlcd_user.py file and gives examples on how a user can run NLCD  
-> `cit_run.r`: This file can be used to reproduce the results obtained from cit  

## Usage 

Open the user_run.ipynb file, and based on the need call nlcd_single() function or nlcd_batch()  
L is the instrument variable, A and B are the gene pairs, shuffle is the number of times to do the permutation testing, and algorithm can be SVR, KRR and ANN, seed can also be set to reproduce the results  
`nlcd_single(L,A,B,shuffles,algo,sample_seed=None,verbose=True,reverse=False)`: This function is to test a single trio, reverse parameter checks the relation from $B \rightarrow A$, the default is testing from $A \rightarrow B$. If no sample seed is provided, the program will generate a sample seed.  
`nlcd_batch(data,shuffles,algo,reverse=False,sample_seed=None)`: This function can test for a set of trios. data variable contains the set of trios to be tested. This function has a parent seed which we can pass as sample_seed, for each of the trios separate seed (child seed) is produced so that their values can be reproduced using `nlcd_single()` function.     
All the yeast data and files that are big can be found this in drive link : https://drive.google.com/drive/folders/1kuUakCynxg145uJckRI7spIvj4JaGgJO?usp=sharing  
## Reproduction of Results and Figures  
Instructions are provided in the results folder and figures folder on how to reproduce the results and figures for both cit and nlcd.    
