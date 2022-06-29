#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/


for file in *_GLUE_plus_SCENIC_default_50.pbs
do
    echo $file
    sbatch $file
    sleep 0.1s
done


for file in *_GLUE_plus_SCENIC_100.pbs
do
    echo $file
    sbatch $file
    sleep 0.1s
done