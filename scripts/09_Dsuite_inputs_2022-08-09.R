###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of the script are:
### To prepare input files for Dsuite analyses.
### Good resources:
### https://github.com/millanek/Dsuite
### https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data

## PART 1: Getting ready ----

## Packages:
library(ape)
library(phytools)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Folders:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group"
setwd(path)
library(here)

## Data:
dataset <- "inornatusgr_R1_c90_n281_m70"
#dataset <- "inornatusgr_R1_c90_n264_m65"
#dataset <- "inornatusgr_R1_c90_n268_m80"

## Taxonomic scheme:
#scheme <- "CANDIDATE_SP_I"
scheme <- "CANDIDATE_SP_III"

## Outgroup:
outgroup <- "essingtonii"

## Folder:
dir.create(path = paste0("Dsuite/", dataset))

## PART 2: Assignments file for Dsuite ----

## Sample info:
sample_info <- read.csv(header = TRUE, file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv")
sample_info$CANDIDATE_SP_I <- gsub(sample_info$CANDIDATE_SP_I, pattern = "cf. ", replacement = "")
sample_info$CANDIDATE_SP_III <- gsub(sample_info$CANDIDATE_SP_III, pattern = "cf. ", replacement = "")
sample_info$CANDIDATE_SP_III[sample_info$CANDIDATE_SP_III == "superciliaris-E (N)"] <- "superciliaris-EN"  
sample_info$CANDIDATE_SP_III[sample_info$CANDIDATE_SP_III == "superciliaris-E (S)"] <- "superciliaris-ES"
sample_info$CANDIDATE_SP_III <- gsub(sample_info$CANDIDATE_SP_III, pattern = "-", replacement = "_")

## Samples:
ddRAD_samples <- read.table(header = FALSE, file = paste0("ipyrad/", dataset, "_outfiles/snps.012.indv"))
names(ddRAD_samples) <- "SAMPLE_ID"

## Merge:
sample_info <- merge(ddRAD_samples, sample_info, by = "SAMPLE_ID")

## Taxa to include and exclude:
if (scheme == "CANDIDATE_SP_III") {
  include <- c(outgroup, "spaldingi_NW", "rimacola", "spaldingi_S", "spaldingi_CY", "spaldingi_NE", "spaldingi_TE", "superciliaris_K", 
             "superciliaris_EN", "superciliaris_ES", "superciliaris_W", "inornatus_E", "lateralis", "inornatus_N", "inornatus_S") 
  exclude <- setdiff(unique(sample_info$CANDIDATE_SP_III), include)
} 
if (scheme == "CANDIDATE_SP_I") {
  include <- c(outgroup, "spaldingi_NW", "rimacola", "spaldingi_S", "spaldingi_CY", "spaldingi_NE", "spaldingi_TE", "superciliaris_K", 
             "superciliaris_E", "superciliaris_W", "eutaenius", "lateralis", "inornatus_N", "severus", "fallens", "helenae", "brachyonyx") 
  exclude <- setdiff(unique(sample_info$CANDIDATE_SP_I), include)
}

## Duplicating assignments column to modify it:
sample_info$assignment <- as.character(sample_info[[scheme]])

## Rename outgroup as required by Dsuite:
sample_info$assignment[sample_info$assignment == outgroup] <- "Outgroup"

## Exclude remaining outgroups from the analyses:
sample_info$assignment[sample_info$assignment %in% exclude] <- "xxx"

## If using helenae complex taxa, remove some potentially misidentified ones:
if (scheme == "CANDIDATE_SP_I") {
  misID <- c("SAMR_42918_Ct_saxa", "WAMR_102423_Ct_saxa", "WAMR_139523_Ct_hele", "WAMR_131371_Ct_hele")
  sample_info$assignment[sample_info$SAMPLE_ID %in% misID] <- "xxx"
}

## If needed, reduce the number of outgroup samples to keep only those used in the phylogenetic tree:
# out_keep <- c("NTMR_22188_Ct_essi", "NTMR_22191_Ct_essi") 
# out_drop <- setdiff(sample_info$SAMPLE_ID[sample_info$CANDIDATE_SP_III == outgroup], out_keep)
# sample_info$assignment[sample_info$SAMPLE_ID %in% out_drop] <- "xxx"

## Save:
assignments <- sample_info[c("SAMPLE_ID", "assignment")]
write.table(assignments, file = paste0("Dsuite/", dataset, "/", scheme, "_assignments.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## PART 3: Guide tree for Dsuite ----

## Using SVD Quartets tree:
include[include == outgroup] <- "essingtonii"
if (scheme == "CANDIDATE_SP_III") {
  tree <- read.nexus("paup_SVDQuartets/inornatusgr_n242_m70_all_outgroups_two-supe-E.tre")
} 
if (scheme == "CANDIDATE_SP_I") {
  tree <- read.nexus("paup_SVDQuartets/inornatusgr_n242_m70_paup_brac_fall_hele_seve_no_misID.tre")
}
 
## Function for adding a tip to a tree where a single tip was before:
add_tip <- function(tree, tip, new.tips) {

    ## Find the edge leading to the tip:
    tip_id <- match(tip, tree$tip.label)

    ## Create the new cherry:
    tree_to_add <- ape::stree(length(c(tip, new.tips)))

    ## Naming the tips:
    tree_to_add$tip.label <- c(tip, new.tips)

    ## Add 0 branch length:
    tree_to_add$edge.length <- rep(0, Nedge(tree_to_add))

    ## Binding both trees:
    return(bind.tree(tree, tree_to_add, where = tip_id))
} ## End of function.

## If splitting deep sister clades into OTUs:
#tree <- add_tip(tree, tip = "superciliaris_E", new.tips = "superciliaris_ES")
#tree <- add_tip(tree, tip = "spaldingi_NW", new.tips = "spaldingi_TE")
plot(tree)

## Other edits and save tree:
tree$tip.label[tree$tip.label == "superciliaris_E"] <- "superciliaris_EN"
tree <- keep.tip(tree, include)
tree$tip.label[tree$tip.label == "essingtonii"] <- "Outgroup"
tree$node.label <- NULL ## Remove node support.
plot(tree) ; nodelabels()
tree <- rotate(phy = tree, node = 18)
tree <- rotate(phy = tree, node = 21)
tree <- ladderize(tree, right = FALSE)
write.tree(tree, file = paste0("Dsuite/", dataset, "/", scheme, "_tree.tre"), tree.names = FALSE)

## Save plot order:
include[include == outgroup] <- "Outgroup"
#include <- include[c(1:13, 15, 14)] ## Reorder.
write.table(include, file = paste0("Dsuite/", dataset, "/", scheme, "_plot_order.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## End of script.
