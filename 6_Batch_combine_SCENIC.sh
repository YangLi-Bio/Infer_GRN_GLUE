#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/
data_list="../../main_text_data.txt"
module load R/4.1.0-gnu9.1


cat $data_list | while read line
do
    array=(${line})
    dir=${array[0]}
    echo "Directory: $dir"
    job=$dir
    # echo "Job: $job"
    
    Rscript /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/6_Combine_SCENIC.R /fs/ess/scratch/PCON0022/liyang/stream/benchmarking/GLUE/${dir}_50_default_out/gene_peak_conn.csv /fs/ess/PCON0022/liyang/STREAM/benchmarking/SCENIC/${dir}/TM_RF_Top_50_Thr_0.03_minJI_0.8.regulon.tsv ${dir}_GLUE_eGRNs_default_50
    
    sleep 0.5s
    
    Rscript /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/6_Combine_SCENIC.R /fs/ess/scratch/PCON0022/liyang/stream/benchmarking/GLUE/${dir}_50_default_out/gene_peak_conn.csv /fs/ess/PCON0022/liyang/STREAM/benchmarking/SCENIC/${dir}/TM_RF_Top_50_Thr_0.03_minJI_0.8.regulon.tsv ${dir}_GLUE_eGRNs_100
    
    sleep 0.5s

    # echo -e "#!/bin/bash\n#SBATCH --job-name=GLUE_plus_SCENIC_default_50_${job}\n#SBATCH --time=00:50:59\n#SBATCH --output="GLUE_plus_SCENIC_default_50_${job}.out"\n#SBATCH--account=PCON0022\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=8\n#SBATCH --mem=100GB\n#SBATCH--gpus-per-node=1\n\nset -e\n\nmodule load R/4.0.2-gnu9.1\n\ncd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/\nstart=$(date +%s)\nsleep 5;\n\nRscript /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/6_Combine_SCENIC.R /fs/ess/scratch/PCON0022/liyang/stream/benchmarking/GLUE/${dir}_50_default_out/gene_peak_conn.csv /fs/ess/PCON0022/liyang/STREAM/benchmarking/SCENIC/${dir}/TM_RF_Top_50_Thr_0.03_minJI_0.8.regulon.tsv ${dir}_GLUE_eGRNs_default_50" > "${job}_GLUE_plus_SCENIC_default_50.pbs"
    
    # echo -e "#!/bin/bash\n#SBATCH --job-name=GLUE_plus_SCENIC_100_${job}\n#SBATCH --time=00:50:59\n#SBATCH --output="GLUE_plus_SCENIC_100_${job}.out"\n#SBATCH--account=PCON0022\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=8\n#SBATCH --mem=100GB\n#SBATCH--gpus-per-node=1\n\nset -e\n\nmodule load R/4.0.2-gnu9.1\n\ncd /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/\nstart=$(date +%s)\nsleep 5;\n\nRscript /fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Codes/6_Combine_SCENIC.R /fs/ess/scratch/PCON0022/liyang/stream/benchmarking/GLUE/${dir}_50_default_out/gene_peak_conn.csv /fs/ess/PCON0022/liyang/STREAM/benchmarking/SCENIC/${dir}/TM_RF_Top_50_Thr_0.03_minJI_0.8.regulon.tsv ${dir}_GLUE_eGRNs_100" > "${job}_GLUE_plus_SCENIC_100.pbs"
done