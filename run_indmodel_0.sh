#!/bin/bash
#$ -N indmodel
#$ -cwd
#$ -o outputfile
#$ -e errorfile
python3 --version
python3 CITNonLinear_wandb_gndtruth1.py > output.log
