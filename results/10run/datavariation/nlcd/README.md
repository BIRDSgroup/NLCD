## How to reproduce the results?
For example, lets reproduce SineKRR500s100perm1run.csv. Here the Sine in the file name indicates we are using Sine dataset, KRR as the algorithm, 500 samples and 100 permutations and this is on the first dataset among the 10 generated datsets  
Open the file or just view the first few lines to know the parent seed  
![This image shows how to reproduce NLCD results](https://drive.google.com/uc?export=view&id=1PVBdl1Je6oB4lHbcmqAymdzWteRDnQ9P)
Note down the parent seed marked in red here
Run the below command:
```
python3 nlcd_simulation.py -a=ANN -s=500 -i=inputfiledirectory/10rundata/Sine500run1.txt -o=outputfilename.csv --seed=314585293344357982809773176541730115694
```
Make sure that the input file path is correct and the file is the first Sine data in the 10rundata folder with 1000 samples. The outputfilename can be of user's preference.
Make sure that `nlcd_main.py`, `nlcd_user.py` are also in the same directory as `nlcd_simulation.py`

Based on the filename the algorithm can be ANN, KRR, or SVR. The sample size can be either 1000, 500, or 300. The datasets can be
Linear: Linear dataset
Sine: Sine dataset
Saw: Saw dataset
Para: Parabola dataset
Paravar: Parabola with unequal variance dataset
Linearvar: Linear with unequal varaince dataset
