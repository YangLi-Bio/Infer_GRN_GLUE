#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/


for file in *.GLUE_HM.pbs
do
    echo $file
    sbatch $file
    sleep 0.1s
done
