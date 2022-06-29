###############################################################
#                                                             #
#               Run pathway enrichment analysis               #
#                                                             #
###############################################################


# Set parameters
work.dir <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
data.file <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/main_text_data.txt"


# Get the list of files
source(paste0(tool.dir, "IO_tools.R"))
data.list <- input_txt(data.file, sep = " ")
data.org <- Reduce("rbind", lapply(data.list, function(x) {
  c(x[1], x[2])
}))
data.vec <- data.org[, 2]
names(data.vec) <- data.org[, 1]

setwd(work.dir)
dir.list <- list.dirs()
length(dir.list)
dir.list <- dir.list[2:6]
head(dir.list)
dir.list <- gsub("^\\./", "", dir.list)
library(pbapply)
library(dplyr)
source(paste0(tool.dir, "transcriptome_tools.R"))
pblapply(dir.list, function(dd) {
  org <- data.vec[[dd]]
  ifelse (grepl("^mm", org), species <- "mouse", 
          species <- "human")
  file.list <- setdiff(list.files(dd, pattern = ".qsave"), 
                       list.files(dd, pattern = "_KEGG.qsave")) %>% 
    gsub("\\.qsave", "", .) %>% 
    setdiff(., gsub("_KEGG.qsave", "", 
                    list.files(dd, pattern = "_KEGG.qsave"))) %>% 
    paste0(., ".qsave")
  lapply(file.list, function(ff) {
    res <- qs::qread(paste0(work.dir, dd, "/", ff))
    if (is.null(res) | length(res) < 1) {
      return(NULL)
    }
    genes.list <- lapply(res, function(x) {
      unique(x$gene)
    })
    enriched <- run_GO_and_KEGG(genes.ll = genes.list, org = species, 
                                dbs = "KEGG")
    file.enriched <- paste0(work.dir, dd, "/", 
                            gsub("\\.qsave", "_KEGG.qsave", ff))
    qs::qsave(enriched, file.enriched)
  })
})