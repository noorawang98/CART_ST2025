library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(homologene)
library(Scissor)
source('~/scripts/Scissor5.R')
source('~/LJX/code/schard-main/R/functions.R')
source('~/LJX/code/schard-main/R/h5ad_util.R')

setwd('/home/Data/ST_20241110/spatial_input/h5ad/recluster_obj/')

# library(dplyr)
library(purrr)
# library(tibble)

# ## 读取之前保存的 infos2_hm2mm_list
load("./Scissor_rdata/RData/Scissor_5bulk_all_samples.RData")
# 
# df = read.csv('./Scissor_rdata/NR1/meta_GSE248835_NR1.Scissor_260520_alpha0.005.csv')
# 
eval_dir <- "./Scissor_rdata/real_cell_eval"
dir.create(eval_dir, showWarnings = FALSE, recursive = TRUE)

auc_dir <- "./Scissor_rdata/real_cell_auc"
dir.create(auc_dir, showWarnings = FALSE, recursive = TRUE)

source('~/LJX/code/run_one_reliability.R')

# # test single samples
# test_res <- run_one_reliability(
#   scissor_dir = filename,
#   cohort_name = "GSE197977",
#   sample_name = "NR1",
#   infos2_hm2mm_list = infos2_hm2mm_list,
#   n = 100,
#   nfold = 10,
#   skip_if_exists = FALSE
# )

## =========================
## 1. 批量 reliability.test
## =========================
log_file <- file.path(eval_dir, "rstudio_background_job_progress.log")

write_log <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste0(..., collapse = "")
  )
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

all_summary <- list()

for (cohortname in names(infos2_hm2mm_list)) {
  
  for (samplename in names(infos2_hm2mm_list[[cohortname]])) {
    
    one_key <- paste(cohortname, samplename, sep = "__")
    
    write_log("Start: ", one_key)
    
    one_res <- tryCatch({
      
      run_one_reliability(
        scissor_dir = "./Scissor_rdata/RData",
        cohort_name = cohortname,
        sample_name = samplename,
        infos2_hm2mm_list = infos2_hm2mm_list,
        n = 100,
        nfold = 10,
        skip_if_exists = TRUE
      )
      
    }, error = function(e) {
      
      write_log("ERROR: ", one_key, " | ", conditionMessage(e))
      
      data.frame(
        cohort = cohortname,
        sample = samplename,
        cell_num = NA_integer_,
        n_scissor_pos = length(infos2_hm2mm_list[[cohortname]][[samplename]]$Scissor_pos),
        n_scissor_neg = length(infos2_hm2mm_list[[cohortname]][[samplename]]$Scissor_neg),
        alpha = infos2_hm2mm_list[[cohortname]][[samplename]]$para$alpha,
        lambda = infos2_hm2mm_list[[cohortname]][[samplename]]$para$lambda,
        family = infos2_hm2mm_list[[cohortname]][[samplename]]$para$family,
        auc = NA_real_,
        status = paste0("error: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    })
    
    all_summary[[one_key]] <- one_res
    
    write.csv(
      dplyr::bind_rows(all_summary),
      file.path(eval_dir, "real_cell_eval_summary_running.csv"),
      row.names = FALSE
    )
    
    write_log("Done: ", one_key)
    
    gc()
  }
}

eval_summary <- dplyr::bind_rows(all_summary)

write.csv(
  eval_summary,
  file.path(eval_dir, "real_cell_eval_summary_all.csv"),
  row.names = TRUE
)

write_log("All finished.")