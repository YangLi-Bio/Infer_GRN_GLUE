############################################################
#                                                          #
#      Combine the CRE-gene lists (GLUE) with TF-gene      #
#      relations (SCENIC)                                  #
#                                                          #
############################################################


# Load parameters
args <- commandArgs(T)
p.GLUE <- args[1]
# p.GLUE <- "/fs/ess/scratch/PCON0022/liyang/stream/benchmarking/GLUE/SNARE_seq_CellLineMixture_50_default_out/gene_peak_conn.csv"

p.SCENIC <- args[2]
# p.SCENIC <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/SCENIC/SNARE_seq_CellLineMixture/TM_RF_Top_50_Thr_0.03_minJI_0.8.regulon.tsv"

p.output <- args[3] # output file path

message ("Loading CRE-gene linkages from ", p.GLUE, 
         " and TF-gene relations from ", p.SCENIC, " ...\n")


# Load data
library(dplyr)
CRE.GeneLinks <- read.csv(p.GLUE, header = T)
TF.GenePairs <- read.table(p.SCENIC, header = T)[, c(1, 2)]
message ("Loaded ", nrow(CRE.GeneLinks), " CRE-gene linkages.\n")
message ("Loaded ", nrow(TF.GenePairs), " TF-gene pairs.\n")


# Combine the CRE-gene linkages with TF-gene pairs
TFGene.list <- split(TF.GenePairs, f = TF.GenePairs$TF) %>% lapply(., "[[", 2)
TF.links <- lapply(TFGene.list, function(x) {
  CRE.GeneLinks[CRE.GeneLinks$gene %in% x,]
})
message ("Identified CRE-gene linkages associated with ", 
         length(TF.links), " TFs.\n")


# Save the qsave file
qs::qsave(TF.links, paste0(p.output, ".qsave"))
