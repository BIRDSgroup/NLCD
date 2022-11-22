#The first argument is for the total runs and the second is batches of runs
for ((i = 0; i < $1; i+=$2 ));
do
  qsub -V -v runs=$2,ini=$i run_indmodel_0.sh 
done
