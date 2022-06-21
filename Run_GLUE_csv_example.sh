#!/bin/bash
#SBATCH --job-name=GLUE_csv_example
#SBATCH --time=12:20:59
#SBATCH --output=GLUE_csv_example.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/


module load python/3.7-2019.10
source activate GLUE_env
python Codes/1_GLUE_preprocessing.py /fs/ess/scratch/PCON0022/liyang/stream/csv_files/SNARE_seq_CellLineMixture_RNA.csv /fs/ess/scratch/PCON0022/liyang/stream/csv_files/SNARE_seq_CellLineMixture_ATAC.csv hg38 csv_example_out

python Codes/2_Model_training.py csv_example_out
python Codes/4_Infer_CRE_gene_links.py csv_example_out
conda deactivate
