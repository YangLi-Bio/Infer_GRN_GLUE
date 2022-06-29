############################################################################
#                                                                          #
#          Calculate the average precision, recall, and f1 score           #
#                                                                          #
############################################################################


# Parameters

# library(enrichR)
source("/fs/ess/PCON0022/liyang/STREAM/param_tuning/codes/library.R") # load packages
source("/fs/ess/PCON0022/liyang/STREAM/param_tuning/codes/input.R")
# data loading, dimension reduction, and clustering

source("/fs/ess/PCON0022/liyang/STREAM/param_tuning/codes/cisAnalysis.R") # cis-regulatory analysis
source("/fs/ess/PCON0022/liyang/STREAM/param_tuning/codes/globalVar.R") # define global variables
getGlobalVar() # build global variables
# setwd("/fs/ess/PCON0022/liyang/STREAM/benchmarking/STREAM")
# pval.cutoff <- 0.05
if.merge <- T


############################################################################
#                                                                          #
#                Function to evaluate the CREs in each eGRN                #
#                                                                          #
############################################################################


# Evaluate each dataset
assess_similarity <- function(data.name = "10x_e18_mouse_brain_fresh_5k", 
                              param.name = "hbc_15_var_3000_top_3000_cons_0.80.extended_egrns.qsave", 
                              bed.ll, if.merge = T) {
  
  
  # data.name : the name of a dataset, e.g., data.name <- "10x_e18_mouse_brain_fresh_5k"
  # param.name : the name of a parameter setting, e.g., param.name <- "hbc_15_var_3000_top_3000_cons_0.80.extended_egrns.qsave"
  # bed.ll : the list of ChIP-Seq datasets with histone marks
  # if.merge : whether merge all the eGRNs, T by default

  
  # Load the eGRNs
  egrns.ll <- qs::qread(paste0(data.name, "/", param.name))
  message ("----> There are in total ", length(egrns.ll), " eGRNs.\n")
  
  
  # # Information for sampling
  # chr.ll <- org.gs@seqinfo@seqlengths
  # names(chr.ll) <- org.gs@seqinfo@seqnames
  
  
  # Evaluate the enrichment of CREs of each eGRN against H3K27ac histone marks
  if ("peak.status" %in% names(egrns.ll[[1]])) {
    peaks.ll <- lapply(egrns.ll, function(x) {
      x[["peaks"]][unlist(x[["peak.status"]])]
    })
  } else if ("peaks" %in% names(egrns.ll[[1]])) {
    peaks.ll <- lapply(egrns.ll, "[[", "peaks")
  } else if ("data.frame" %in% class(egrns.ll[[1]])) {
    peaks.ll <- lapply(seq_along(egrns.ll), function(i) {
      unique(egrns.ll[[i]]$peak)
    })
  }

  peaks.ll <- peaks.ll[!sapply(peaks.ll, is.null)] # peaks.ll <- peaks.ll[1:5]
  peaks.ll <- peaks.ll[sapply(peaks.ll, length) > 0]
  if (if.merge) {
    peaks.temp <- Reduce("union", peaks.ll)
    peaks.ll <- c()
    peaks.ll[[1]] <- peaks.temp
  }
  # peaks.ll <- peaks.ll[1:5]
  rm(egrns.ll)
  message ("----> Began evaluating each eGRN against H3K27ac histone marks ...\n")
  precision.recall.dt <- rbindlist(pbmclapply(seq_along(peaks.ll), function(i) { # pp <- peaks.ll[[1]]
    pp <- peaks.ll[[i]]
    message ("--------> Processing the ", i, "-th eGRN ...\n")
    peaks.gr <- StringToGRanges(pp)
    overlap.dt <- rbindlist(lapply(seq_along(bed.ll), function(j) { # j <- 1
      # message ("-----------> Comparing the eGRN with the ", j, "-th peak list ...\n")
      if ("data.frame" %in% class(bed.ll[[j]])) {
        bb <- makeGRangesFromDataFrame(bed.ll[[j]])
      } else {
        bb <- bed.ll[[j]]
      }
      query.subject <- findOverlaps(peaks.gr, bb)

      # random.ll <- vector()
      # set.seed(123)
      # for (i in 1:n.samples) {
      #   message (i, "\n")
      #   random.gr <- sample_granges(original.gr = peaks.gr, chr.ll = chr.ll)
      #   random.ll <- c(random.ll, findOverlaps(random.gr, bb) %>% length)
      # }
      # length(which(random.ll < query.subject))
      
      # overlaps <- min(query.subject %>% subjectHits %>% unique %>% length, 
      #                 query.subject %>% queryHits %>% unique %>% length)
      # pval <- phyper(overlaps, length(bb), n.bed - length(bb), length(peaks.gr), lower.tail = F)
      # return(list(eGRN = i, bed = j, pval = pval, overlap = overlaps, 
      #             l.query = length(peaks.gr), l.subject = length(bb)))
      
      overlap.query <- query.subject %>% queryHits %>% unique
      overlap.subject <- query.subject %>% subjectHits %>% unique
      list(query = length(overlap.query) / length(pp), subject = length(overlap.subject) / length(bb))
    }))
    overlap.dt[which.max(overlap.dt$query),]
  }, mc.cores = detectCores()))
  message ("----> Finished evaluating ", nrow(precision.recall.dt), " eGRN against H3K27ac histone marks.\n")
  if (nrow(precision.recall.dt) > 0) {
    return(apply(precision.recall.dt, 2, mean))
  } else {
    return(data.frame(precision = 0, recall = 0, f1.score = 0))
  }
}


############################################################################
#                                                                          #
#                                Main program                              #
#                                                                          #
############################################################################


# Get parameters
args <- commandArgs(T)
data.name <- args[1] # the name of a dataset, e.g., data.name <- "10x_e18_mouse_brain_fresh_5k"
# data.name <- "10x_pbmc_unsorted_3k"
# data.name <- "10x_pbmc_granulocyte_sorted_3k"


# # Annotation databases
# org.pathway <- orgPathway.hash[[org]] # organism ID in KEGG
# org.gs <- orgGS.hash[[org]] # genome sequences
# org.anno <- orgAnno.hash[[org]] # the annotation
# org.db <- orgDB.hash[[org]] # database


# # Convert characters into objects
# specLib(org, org.anno, org.gs, org.db) # load specific libraries
# org.gs <- get(org.gs)
# org.anno <- get(org.anno)
# org.db <- get(org.db)


# Load the Cistrome databases
cistrome.dir <- "/fs/ess/scratch/PCON0022/liyang/stream/comparison/H3K27ac_chip_seq/processed_bed/"
file.ll <- list.files(paste0(cistrome.dir, data.name))
bed.ll <- pblapply(seq_along(file.ll), function(i) {
  ff <- file.ll[[i]]
  bed <- qs::qread(paste0(cistrome.dir, data.name, "/", ff))
  bed
}) # bed.ll <- bed.ll[1:10]
# tail(bed.ll)


# Get the list of parameter settings
param.ll <- setdiff(gsub(".qsave", "", list.files(path = data.name, pattern = ".qsave")), 
                    gsub(".qsave", "", list.files(path = data.name, 
                                                                pattern = "_KEGG.qsave"))) %>% 
            paste0(., ".qsave") # param.ll <- param.ll[1:3]


# Evaluate the eGRNs
data.precision.recall <- Reduce("rbind", lapply(param.ll, function(x) { # x <- param.ll[[1]]
  message ("----> Dataset: ", x, " ...\n")
  assess_similarity(data.name = data.name, param.name = x, bed.ll = bed.ll, if.merge = if.merge) # process one parameter setting
})) %>% as.data.frame
rownames(data.precision.recall) <- gsub(".qsave", "", param.ll)
colnames(data.precision.recall) <- c("precision", "recall")
if (if.merge) {
  qs::qsave(data.precision.recall, paste0(data.name, ".merged_hubs.qsave"))
} else {
  qs::qsave(data.precision.recall, paste0(data.name, ".extended_hubs.qsave")) 
}