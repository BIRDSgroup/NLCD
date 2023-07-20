Internal message: Code ran in ibse actinium server
## How to reproduce the results?
To reproduce AtoB.csv. 
Open the file or just view the first few lines to know the parent seed
![This image shows how to reproduce NLCD results](https://drive.google.com/uc?export=view&id=1WOmDOoz18-iyGBZb8gO2d_TAUgUbMjce)
Note down the parent seed marked in red here
Run the below command:
```
python3 nlcd_simulation.py -a=ANN -s=500 -i=inputfiledirectory/human_muscle.txt -o=AtoB.csv --seed=196824294903844221866610020071140119398
```
To reproduce AtoB_rev.csv. 
Open the file or just view the first few lines to know the parent seed
![This image shows how to reproduce NLCD results](https://drive.google.com/uc?export=view&id=1S0sHtyqPx0dcliOlyLf9VaVo3Y200Olt)
Note down the parent seed marked in red here
Run the below command:
```
python3 nlcd_simulation.py -a=KRR -s=100 -i=inputfiledirectory/human_muscle.txt -o=BtoA_rev.csv --seed=233362693158198143346700831088686120739 --reverse
```
Make sure that the input file path is correct, this file is present in the drive link provided. The outputfilename can be of user's preference.  
Make sure that `nlcd_main.py`, `nlcd_user.py` are also in the same directory as `nlcd_simulation.py`  
Now from the AtoB.csv and BtoA.csv extract the cis-trans and trans-cis trios. For this run the `extract_trios.ipynb`  

Run the `convert_hgnc.r` file to get the gene names from the ensembl names.  
