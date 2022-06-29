##################################################################
#                                                                #
#         Calculate the precision, recall, and f-score           #
#                                                                #
##################################################################


# Parameters
library(dplyr)
args <- commandArgs(T)
data <- args[1] # data <- "SNARE_seq_CellLineMixture" ; data <- "SHARE_seq_Brain"
org <- args[2] # org <- "hg38" ; org <- "mm10"
files.ll <- setdiff(gsub(".qsave", "", list.files(path = data, pattern = ".qsave")), 
                    gsub(".qsave", "", list.files(path = data, 
                                                  pattern = "_KEGG.qsave"))) %>% 
  paste0(., ".qsave") # param.ll <- param.ll[1:3]
if (length(files.ll) < 1) {
  files.ll <- list.files(path = data, pattern = ".regulon.tsv")
}
if (length(files.ll) < 1) {
  stop ("There is no successful parameter setting!\n")
}


# Collect enhancer-target pairs in databases
atlas.dir <- "/fs/ess/scratch/PCON0022/liyang/stream/comparison/EnhancerAtlas2/"
atlas.files <- list.files(paste0(atlas.dir, org), pattern = ".qsave")
message ("There are in total ", length(atlas.files), " EnhancerAtlas 2.0 files.\n")


# Load the enhancer-target pairs
library(pbmcapply)
library(qs)
ep.ll <- pbmclapply(atlas.files, function(ee) {
  qs::qread(paste0(atlas.dir, org, "/", ee))
}, mc.cores = detectCores())


# Calculate peak-gene linkages in eGRNs
source("/fs/ess/PCON0022/liyang/r_utilities/functions/cistrome_tools.R")
library(pbapply)
library(pbmcapply)
file.scores <- Reduce("rbind", pbmclapply(files.ll, function(ff) {
  message ("----> Processing the file: ", ff, " ...\n")
  egrns.ll <- qs::qread(paste0(data, "/", ff))
  egrns.ll <- lapply(egrns.ll, function(x) {
    list(genes = unique(x$gene), peaks = unique(x$peak))
  })
  egrns.ll <- egrns.ll[sapply(egrns.ll, function(ee) {
    if (length(ee$genes) < 1 | length(ee$peaks) < 1) {
      return(F)
    }
    return(T)
  })]
  if (length(egrns.ll) < 1) {
    return(c(precision = NA, recall = NA, fscore = NA))
  }
  if ("peak.status" %in% names(egrns.ll[[1]])) {
    egrns.ll <- lapply(egrns.ll, function(xxx) {
      list(genes = xxx$genes, peaks = xxx$peaks[unlist(xxx$peak.status)])
    })
  }
  links <- link_peaks_to_genes_in_batch(pairs.ll = egrns.ll, org = org)
  enrich.scores <- Reduce("rbind", pblapply(seq_along(links), function(i) {
    message ("--------> Processing the ", i, "-th eGRN ...\n")
    xx <- links[[i]]
    score.dt <- enhancer_target_precision_recall_fscore(link.pairs = xx, ep.ll = ep.ll, 
                                                        only.overlap = T)
    added.dt <- cbind(link = rep(i, nrow(score.dt)), score.dt)
    added.dt
  }))
  enrich.scores <- na.omit(enrich.scores)
  if (nrow(enrich.scores) < 1) {
    return(c(precision = 0, recall = 0, fscore = 0))
  } else {
    return(apply(enrich.scores[, c("precision", "recall", "fscore")], 2, mean))
  }
}, mc.cores = detectCores()))
rownames(file.scores) <- gsub(".qsave", "", files.ll)
file.scores <- na.omit(file.scores)
# qs::qsave(file.scores, paste0(data, ".EP_scores.qsave"))
qs::qsave(file.scores, paste0(data, ".enhancer_precision.qsave"))
message ("Finished calculating precision, recall, and f-score for ", data, ".\n")