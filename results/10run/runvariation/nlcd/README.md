## How to reproduce the results?
For example, lets reproduce IndpKRR500s500perm6run.csv. Here Independent in the file name indicates we are using Independent dataset, KRR as the algorithm, 500 samples and 500 permutations.
Open the file or just view the first few lines to know the parent seed  
![This image shows how to reproduce NLCD results](https://drive.google.com/uc?export=view&id=1DiR767GqP7xjHP88slQoJvpMGfPlIBkN)
Note down the parent seed marked in red here
Run the below command:
```
python3 nlcd_simulation.py -a=ANN -s=500 -i=inputfiledirectory/Linear1000.txt -o=outputfilename.csv --seed=198324164150862344355037048063459674793
```
Make sure that the input file path is correct and the file is the independent dataset with 500 samples. The outputfilename can be of user's preference.
Make sure that `nlcd_main.py`, `nlcd_user.py` are also in the same directory as `nlcd_simulation.py`

Based on the filename the algorithm can be ANN, KRR, or SVR. The sample size can be either 1000, 500, or 300. The datasets can be
Linear: Linear dataset
Sine: Sine dataset
Saw: Saw dataset
Para: Parabola dataset
Paravar: Parabola with unequal variance dataset
Linearvar: Linear with unequal varaince dataset
