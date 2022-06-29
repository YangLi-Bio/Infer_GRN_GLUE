#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs

data_list="../../main_text_data.txt"

cat $data_list | while read line
do
    array=(${line})
    dir=${array[0]}
    echo "Directory: $dir"
    
    job=$dir
    echo "Job: $job"
    
    # rds=${dir/_dir/.RDS}
    
    org=${array[1]}
    echo "Organism: $org"
    
    # cd $dir
    echo -e "#!/bin/bash\n#SBATCH --job-name=GLUE_EP_${job}\n#SBATCH --time=11:20:59\n#SBATCH --output="GLUE_${job}_EP.out"\n#SBATCH --account=PCON0022\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=80GB\n#SBATCH --gpus-per-node=1\n\nset -e\n\nmodule load R/4.1.0-gnu9.1\n\ncd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/\nstart=\$(date +%s)\nsleep 5;\n\nRscript ../Codes/10_Evaluate_EP_enrich.R ${dir} ${org}" > "${job}.GLUE_EP.pbs"
    # cd ..
    sleep 0.1s
done
