###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To estimate missing data and prepare SNP datasets for genetic clustering analyses.

## PART 1: Getting ready ----

## Packages:
library(LEA)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/"
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## PART 2: Writing some functions we'll need ----

## Function: Estimating missing data per SNP position.
get_miss_by_SNP <- function(SNP) {

  ## Let's now estimate levels of missing data per sample:
  print(paste0("Now processing SNP ", SNP))

  ## How many sites have no data?
  miss <- table(gendata[, SNP] == 9)["TRUE"]
  
  ## What proportion this is from the total?
  if (is.na(miss) == TRUE) { SNP_miss <- 0 
  } else { SNP_miss <- round(miss/(dim(gendata)[1]), digits = 2) }
  
  ## Storing:
  SNP_miss <- data.frame(SNP = SNP, miss_prop = SNP_miss)
  SNP_miss$SNP <- as.character(SNP_miss$SNP)
  
  ## Return:
  return(SNP_miss)
  
} ## End of function.
    
## Function: Excluding SNPs with too much missing data.
prep_SNPs <- function(miss_SNP) {
  
  ## Testing:
  #miss_SNP <- 0.7
  
  ## Read missing information for SNPs:
  SNP_df <- read.csv(file = here(paste0("ipyrad/", dataset, "_outfiles/missing_data/missing_by_SNP.csv")), header = TRUE)
  
  ## Which ones have more than maximum allowed?
  SNPs_exclude <- SNP_df$SNP[SNP_df$miss_prop > miss_SNP]
  
  ## Exclude:
  gen_clean <- gendata[setdiff(names(gendata), SNPs_exclude)]
  #dim(gen_clean)
  
  ## Import sample IDs:
  sample_IDs <- read.table(here(paste0("ipyrad/", dataset, "_outfiles/usnps.012.indv")), header = FALSE)
  gen_clean$SAMPLE_ID <- sample_IDs$V1
  gen_clean <- gen_clean[c("SAMPLE_ID", setdiff(names(gen_clean), "SAMPLE_ID"))]
  
  ## Save:
  write.table(gen_clean, here(paste0("ipyrad/", dataset, "_outfiles/ms", miss_SNP, ".usnps")), row.names = FALSE)
  
} ## End of function.

## Function: Estimating levels of missing data per sample.
get_miss_by_ind <- function(miss_SNP) {

  ## Testing:
  #miss_SNP <- 0.7
  
  ## Genetic data after removing poorly sampled SNPs:
  gen_clean <- read.table(file = here(paste0("ipyrad/", dataset, "_outfiles/ms", miss_SNP, ".usnps")), header = TRUE)
  dim(gen_clean)
  
  ## List to save results:
  list_ind <- vector("list", length(gen_clean$SAMPLE_ID))
  names(list_ind) <- gen_clean$SAMPLE_ID
  
  ## Status:
  print("Now estimating missing data for samples!")
  
  ## Function:
  for (sample in names(list_ind)) {
      
    ## Testing:
    #sample <- "WAMR_166392_Ct_hele"
    
    ## Let's now estimate levels of missing data per sample:
    print(paste0("Now processing sample ", sample))
      
    ## Get SNPs for sample:
    sample_data <- gen_clean[gen_clean$SAMPLE_ID == sample, ]
    
    ## What proportion is this from the total number of loci? 
    miss <- table(sample_data == 9)["TRUE"]
    ind_miss <- round(miss/(dim(gen_clean)[2]), digits = 2)
    
    ## Storing:
    ind_miss <- data.frame(SAMPLE_ID = sample, miss_prop = ind_miss)
    ind_miss$SAMPLE_ID <- as.character(ind_miss$SAMPLE_ID)
    list_ind[[sample]] <- ind_miss
    
  } ## Close loop.
  
  ## List to dataframe:
  ind_df <- do.call("rbind", list_ind)
  
  ## Save this information:
  write.csv(ind_df, file = here(paste0("ipyrad/", dataset, "_outfiles/missing_data/missing_by_sample_ms", miss_SNP, ".csv")), row.names = FALSE, quote = FALSE)
  
} ## End of function.

## Function: Filtering individuals based on maximum missing data allowed.
prep_inds <- function(miss_ind, miss_SNP, max_n) {

  ## Import data:
  gen_clean <- read.table(file = here(paste0("ipyrad/", dataset, "_outfiles/ms", miss_SNP, ".usnps")), header = TRUE)
  dim(gen_clean)
  
  ## Missing data per individual:
  ind_df <- read.csv(file = here(paste0("ipyrad/", dataset, "_outfiles/missing_data/missing_by_sample_ms", miss_SNP, ".csv")), header = TRUE)
  
  ## What samples have less than max_miss missing data?
  inds_to_keep <- gen_clean$SAMPLE_ID[ind_df$miss_prop <= miss_ind]
  
  ## Keep samples with less than max_miss missing data:
  gen_clean <- gen_clean[gen_clean$SAMPLE_ID %in% inds_to_keep, ]
  dim(gen_clean)
  
  ## Now let's rarefy overly sampled sites.
  ## First, read sample information:
  sample_info <- read.csv(file = paste0(path, "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv"), header = TRUE)
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% gen_clean$SAMPLE_ID, ]
  
  ## How many samples per site?
  site_list <- group_by(sample_info, UPDATED_SP, LAT, LON) %>% group_split()
  n_per_site <- purrr::map(.x = site_list, .f = nrow)
   
  ## Sites that have less or equal the maximum number of samples:
  norm_sampled <- site_list[n_per_site <= max_n]
  norm_sampled <- dplyr::bind_rows(norm_sampled)
   
  ## Oversampled sites:
  over_sampled <- site_list[n_per_site > max_n]
  over_sampled_filtered <- purrr::map_dfr(over_sampled, sample_n, size = max_n)
   
  ## Combine "normal" and filtered oversampled sites:
  sample_info <- rbind(norm_sampled, over_sampled_filtered)
  
  ## Let's only keep genetic data for the individuals that passed all filters:
  gen_clean <- gen_clean[gen_clean$SAMPLE_ID %in% sample_info$SAMPLE_ID, ]
  dim(gen_clean)
  
  ## Function: Remove invariant sites:
  get_var_SNP <- function(SNP) {
  if (length(unique(gen_clean[[SNP]])) >= 2) { 
    SNP_df <- data.frame(SNP = gen_clean[[SNP]]) 
    names(SNP_df) <- SNP
  } else { SNP_df <- NULL }
  return(SNP_df) } ## End of function.
  gen_only <- gen_clean[setdiff(names(gen_clean), "SAMPLE_ID")] ; dim(gen_only)
  gen_only <- map_dfc(names(gen_only), get_var_SNP) ; dim(gen_only)
  gen_only$SAMPLE_ID <- gen_clean$SAMPLE_ID ; dim(gen_only)
  dim(gen_only)
  
  ## Reorder and saving data with Sample IDs:
  gen_clean <- gen_only[c("SAMPLE_ID", setdiff(names(gen_only), "SAMPLE_ID"))]
  write.table(gen_clean, here(paste0("ipyrad/", dataset, "_outfiles/mi", miss_ind, "_ms", miss_SNP, ".usnps")), row.names = FALSE)
  
} ## End of function.

## PART 3: Use functions we wrote to estimate missing data and generate filtered SNP datasets ----

## Genetic dataset we'll use (ipyrad output):
#dataset <- "spaldingi_complex_n72_m50"
#dataset <- "superciliaris_complex_n46_m50"
dataset <- "inornatus_complex_n103_m50"

## New folder:
dir.create(here(paste0("ipyrad/", dataset, "_outfiles/missing_data")))

## Reading genetic dataset:  
gendata <- read.table(here(paste0("ipyrad/", dataset, "_outfiles/usnps.012")), sep = "\t", row.names = 1)
dim(gendata)

## Estimate levels of missing data per SNP position and save results:
print("Now estimating missing data for SNPs!")
SNP_df <- map_df(names(gendata), get_miss_by_SNP)
write.csv(SNP_df, file = here(paste0("ipyrad/", dataset, "_outfiles/missing_data/missing_by_SNP.csv")), row.names = FALSE, quote = FALSE)

## Exclude SNPs with too much missing data:
for (miss_SNP in (c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8))) { prep_SNPs(miss_SNP) }

## Estimate levels of missing data per sample under each missing SNP level:
print("Now estimating missing data for samples!")
for (miss_SNP in (c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8))) { get_miss_by_ind(miss_SNP) }

## Finally, generate filtered genetic datasets under the maximum levels of missing data by SNP position and individual:
for (miss_ind in (c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8))) {
  for (miss_SNP in (c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8))) {
    prep_inds(miss_ind, miss_SNP, max_n = 5)
  }
}

## End of script.
