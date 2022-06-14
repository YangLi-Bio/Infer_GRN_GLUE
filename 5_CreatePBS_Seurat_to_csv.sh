#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/

data_list="main_text_data.txt"

cat $data_list | while read line
do
    array=(${line})
    dir=${array[0]}
    echo "Directory: $dir"
    
    job=$dir
    echo "Job: $job"

    echo -e "#!/bin/bash\n#SBATCH --job-name=Seurat_to_csv_${job}\n#SBATCH --time=01:50:59\n#SBATCH --output="Seurat_to_csv_${job}.out"\n#SBATCH--account=PCON0022\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=8\n#SBATCH --mem=40GB\n#SBATCH--gpus-per-node=1\n\nset -e\n\nmodule load R/4.0.2-gnu9.1\n\ncd /fs/ess/PCON0022/liyang/STREAM/benchmarking/\nstart=$(date +%s)\nsleep 5;\n\nRscript /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/5_Seurat_to_csv.R ${dir}.RDS" > "${job}_Seurat_to_csv.pbs"
    # cd ..
done