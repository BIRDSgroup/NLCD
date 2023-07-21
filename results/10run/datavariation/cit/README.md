## How to reproduce the results?
For example, lets reproduce Linear500cit100permrun8.csv. Here the Linear in the file name indicates we are using Linear dataset, 500 samples and 100 permutations and this is running cit on the 8th linear dataset out of the 10 datasets generated. 
Open the file or just view the first few lines to know the parent seed  
![This image shows how to reproduce NLCD results](https://drive.google.com/uc?export=view&id=1TgmfOzrUM2_RoTfzWt9MZhjSusc6x_Da)
Note down the parent seed marked in red here
Open the `cit_run.r`. Compile `read_data()` function and `cit_process()` function.
Call `cit_process()` with the following parameters:
```
cit_process("./data/10rundata/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/filename.csv',seed=1291353810)
```
Make sure that the input file path is correct and the file is among the 10 generated linear datasets with 500 samples. The outputfilename can be of user's preference.

Based on the filename the algorithm can be ANN, KRR, or SVR. The sample size can be either 1000, 500, or 300. The datasets can be
Linear: Linear dataset
Sine: Sine dataset
Saw: Saw dataset
Para: Parabola dataset
Paravar: Parabola with unequal variance dataset
Linearvar: Linear with unequal varaince dataset
