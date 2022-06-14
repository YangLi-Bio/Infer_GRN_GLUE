#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/


for file in *_Seurat_to_csv.pbs
do
    echo $file
    sbatch $file
    sleep 0.1s
done
