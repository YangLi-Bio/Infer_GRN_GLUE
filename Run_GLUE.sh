#!/bin/bash

# Parameters
rna_file=$1
atac_file=$2
org=$3
out_dir=$4
n_pcs=$5


# Main program
python /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/1_GLUE_preprocessing.py $rna_file $atac_file $org $out_dir $n_pcs
python /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/2_Model_training.py $out_dir
python /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/4_Infer_CRE_gene_links.py $out_dir
