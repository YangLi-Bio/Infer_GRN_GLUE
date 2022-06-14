#!/bin/bash
#SBATCH --job-name=Preprocess_and_model_training_GLUE
#SBATCH --time=24:20:59
#SBATCH --output=Preprocess_and_model_training_GLUE.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=300GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/


module load python/3.7-2019.10
source activate GLUE_env
#python 1_GLUE_preprocess_example.py
#python 2_Model_training_example.py
python 4_Infer_CRE_gene_links_example.py
conda deactivate
