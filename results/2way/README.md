
## How to reproduce the results?
For example, lets reproduce LinearANN1000s500perm.csv. Here the Linear in the file name indicates we are using Linear dataset, ANN as the algorithm, 1000 samples and 500 permutations.
Open the file or just view the first few lines to know the parent seed
![This image shows how to reproduce NLCD results](https://drive.google.com/uc?export=view&id=1xcAosXJzU2EmrDInBuS7N5ZOlstxsbOT)
Note down the parent seed marked in red here
Run the below command:
```
python3 nlcd_simulation.py -a=ANN -s=500 -i=inputfiledirectory/Linear1000.txt -o=outputfilename.csv --seed=8712540936262824741492729377897620802 --reverse
```
Note that we are giving the parameter reverse to say that the testing is done from $B \rightarrow A$   
Make sure that the input file path is correct and the file is the linear dataset with 1000 samples. The outputfilename can be of user's preference.
Make sure that nlcd_main.py, nlcd_user.py are also in the same directory as nlcd_simulation.py

Based on the filename the algorithm can be ANN, KRR, or SVR. The sample size can be either 1000, 500, or 300. The datasets can be
Linear: Linear dataset
Sine: Sine dataset
Saw: Saw dataset
Para: Parabola dataset
Paravar: Parabola with unequal variance dataset
Linearvar: Linear with unequal varaince dataset
