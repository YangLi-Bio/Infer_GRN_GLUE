###############################################################
#                                                             #
#     Distribute the results according to three parameters    #
#                                                             #
###############################################################


# Set parameters
work.dir <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Outputs/"


# Get the list of datasets
setwd(work.dir)
qsave.list <- list.files(pattern = "*.qsave")
length(qsave.list)


# Divide result file
library(pbapply)
library(dplyr)
pblapply(qsave.list, function(x) {
  
  # npcs = 50, 100
  orig.res <- qs::qread(x)
  name.ar <- strsplit(x, split = "_GLUE_eGRNs_")[[1]]
  data.dir <- name.ar[1]
  npcs <- strsplit(name.ar[2], split = "[_|.]")[[1]] %>% 
    tail(n = 2) %>% "[" (1)
  dir.create(data.dir)
  
  
  # Other three parameters
  for (qval in c(0.01, 0.05)) {
    for (dist in c(100000, 150000)) {
      for (score in c(0.30, 0.50)) {
        res <- orig.res[!sapply(orig.res, is.null)]
        res.filtered <- lapply(res, function(y) {
          y[y$dist < dist & y$glue > score & y$qval < qval,]
        })
        res.filtered <- res.filtered[sapply(res.filtered, function(y) {
          nrow(y) > 2
        })]
        file.path <- paste0(data.dir, "/NPC_", npcs, "_qval_", 
                           qval, "_dist_", dist, "_score_", 
                           score, ".qsave")
        qs::qsave(res.filtered, file.path)
      }
    }
  }
})
