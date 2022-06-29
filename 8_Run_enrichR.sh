#!/bin/bash
#SBATCH --job-name=Run_enrichR
#SBATCH --time=11:50:59
#SBATCH --output=Run_enrichR.out
#SBATCH--account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100GB
#SBATCH--gpus-per-node=1

set -e
cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/


module load R/4.1.0-gnu9.1
Rscript 8_Run_enrichR.R
