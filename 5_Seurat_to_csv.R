# Parameters and libraries
args <- commandArgs(T) # 10x_human_brain_3k.RDS
obj.path <- args[1]
prefix <- paste0("/fs/ess/scratch/PCON0022/liyang/stream/csv_files/", 
                 gsub(".RDS", "", obj.path))
source("/fs/ess/PCON0022/liyang/r_utilities/functions/IO_tools.R")


# Write files
Seurat_to_csv(obj = readRDS(obj.path), prefix = prefix, 
               slot = "data")
