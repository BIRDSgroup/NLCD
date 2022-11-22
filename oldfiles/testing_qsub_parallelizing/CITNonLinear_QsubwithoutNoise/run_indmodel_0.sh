#!/bin/bash
#$ -N indmodel
#$ -cwd
#$ -o outputfile
#$ -e errorfile
python3 --version
python3 CITNonLinear_QsubwithNormalization.py $runs $ini 
