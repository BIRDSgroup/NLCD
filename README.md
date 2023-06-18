# NLCD : A method to discover nonlinear causal relations among genes
The repo contains the scripts and data used in the paper : cite our paper.  
The scripts folder contains :  
-> nlcd_main.py : This file has all the functions defined  
-> nlcd_simulation.py : This file is used to generate the simulation data results  
-> nlcd_user.py : This file imports the main file and adds two more functions for the user   
-> user_run.ipynb : This notebook imports the above nlcd_user.py file and gives examples on how a user can run NLCD  

## Usage 

Open the user_run.ipynb file, and based on the need call nlcd_user() function or nlcd_group()  
L is the instrument variable, A and B are the gene pairs, shuffle is the number of times to do the permutation testing, and algorithm can be SVR, KRR and ANN 
nlcd_user(L,A,B,shuffles,algo) : This function is to test a single trio 
nlcd_group(data,shuffles,algo) : This function can test for a set of trios. data variable contains the set of trios to be tested.  
