###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To process and plot mitochondrial cytb trees.

## PART 1: Setting up ----

## Packages:
library(ggtree)
library(patchwork)
library(phytools)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/" ## Suindara.
setwd(path)
#detach("package:here", unload = TRUE)
library(here)

## PART 2: Defining tips to be included in tree ----

## Read tree:
tree <- read.tree(file = "cytb/RAxML_bipartitions.n485.tre")
tree <- ladderize(tree, right = FALSE)
  
## Exclude outgroups:
exclude <- c("CUMV_14681_Ct_scho", "CUMV_14700_Ct_scho", "NTMR_22188_Ct_essi", "NTMR_22191_Ct_essi",
             "SAMAR_33571_Ct_taen", "SAMAR_33727_Ct_taen", #"SAMAR_45424_Le_bipe", "WAMR_111809_Le_bipe",
             "UMMZ_242633_Ct_pant", "UMMZ_242639_Ct_pant", #"WAMR_157958_Le_ips", "WAMR_135152_Le_ips",
             "WAMR_139414_Ct_nigr", "WAMR_139415_Ct_nigr", "SAMR_37942_Ct_atla", "SAMR_57397_Ct_atla",
             "UMMZ_242606_Ct_aust", "UMMZ_242607_Ct_aust", "SAMR_22246_Ct_leon", "SAMR_46944_Ct_leon",
             "SAMAR_33571_Ct_taen", "SAMAR_33727_Ct_taen")
tree <- drop.tip(tree, exclude)

## Add sample info:
tips <- data.frame(SAMPLE_ID = tree$tip.label)
sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
sample_info <- sample_info[sample_info$SAMPLE_ID %in% tree$tip.label, ]
sample_info <- sample_info[c("SAMPLE_ID", "VOUCHER", "REGO", "CANDIDATE_SP_III")]

## Combine:
#sample_info <- merge(tips, sample_info, all.x = TRUE) ; sample_info[is.na(sample_info)] <- "x" ## If using all cytb samples.
sample_info <- merge(tips, sample_info) ## If keeping only samples with ddRAD data.

## Adding burbidgei and the Cape York taxon back:
add_back <- data.frame(CANDIDATE_SP_III = c(rep("spaldingi-CY", 6), rep("burbidgei", 5)),
                       SAMPLE_ID = c("ANWC_R05242_Ct_spal", "ANWC_R05268_Ct_spal", "ANWC_R05240_Ct_spal", 
                                     "ANWC_R05241_Ct_spal", "QM_87533_Ct_spal", "QM_87556_Ct_spal",
                                     "WAMR_136872_Ct_burb", "WAMR_136874_Ct_burb", "WAMR_171066_Ct_burb",
                                     "WAMR_171061_Ct_burb", "WAMR_171060_Ct_burb"))
add_back$VOUCHER <- gsub(add_back$SAMPLE_ID, pattern = "(.+)_Ct_.+", replacement = "\\1")
add_back$REGO <- gsub(add_back$VOUCHER, pattern = "_", replacement = "")
sample_info <- rbind(sample_info, add_back)

## If keeping only samples with ddRAD data or some other subset, prune further:
keep <- tree$tip.label[tree$tip.label %in% sample_info$SAMPLE_ID]
tree <- keep.tip(tree, keep)

## PART 3: Editing tip labels ----

## Change labels:
sample_info$label <- sample_info$VOUCHER
unvouchered <- sample_info$SAMPLE_ID[sample_info$VOUCHER == "no_voucher"]
for (sample in unvouchered) { sample_info$label[sample_info$SAMPLE_ID == sample] <- sample_info$REGO[sample_info$SAMPLE_ID == sample] }
unvouchered <- sample_info$SAMPLE_ID[is.na(sample_info$label)]
for (sample in unvouchered) { sample_info$label[sample_info$SAMPLE_ID == sample] <- sample }
sample_info$label <- gsub(sample_info$label, pattern = "NA_", replacement = "")
sample_info$label <- gsub(sample_info$label, pattern = "_Ct_.+", replacement = "")
sample_info$label <- gsub(sample_info$label, pattern = "SAMR", replacement = "SAMAR")
sample_info$label <- gsub(sample_info$label, pattern = "([A-Z]+)([0-9]+)", replacement = "\\1_\\2")
sample_info$label <- gsub(sample_info$label, pattern = "(R)_([0-9]+)", replacement = "_\\1\\2")
sample_info$label <- gsub(sample_info$label, pattern = "__", replacement = "_")
sample_info$label <- gsub(sample_info$label, pattern = " ", replacement = "_")

## Order sample info by tree position:
sample_info$SAMPLE_ID <- factor(x = sample_info$SAMPLE_ID, levels = tree$tip.label)
sample_info <- arrange(sample_info, SAMPLE_ID)

## Change tree tips:
tree$tip.label <- as.character(sample_info$label)

## PART 4: Setting up clade bars ----

## Define clades to add bars:
list_span <- list()
list_span$"inornatus-N and S, superciliaris-E (S) and W" <- c("WAM_R153812", "CCM_3739")
list_span$"spaldingi-NW" <- c("WAM_R141379", "NTM_R22620")
list_span$"supeciliaris-E (N)" <- c("NTM_R26117", "CCM_2292")
list_span$"C. rimacola" <- c("WAM_R126015", "WAM_R126010")
list_span$"C. lateralis "<- c("SAMA_R55259", "SAMA_R42768")
list_span$"inornatus-E" <- c("SAMA_R55874", "SAMA_R55799")
list_span$"inornatus-S" <- c("WAM_R140720", "WAM_R146913")
list_span$"inornatus-S, superciliaris-E (S) and W" <- c("ABTC_60781", "SAMA_R44367")
list_span$"spaldingi-TE" <- c("NTM_R20378", "NTM_R22166")
list_span$"inornatus-N, superciliaris-K" <- c("CCM_0848", "WAM_R174676")
list_span$"C. burbidgei" <- c("WAM_R136872", "WAM_R171060")
list_span$"spaldingi-S" <- c("AM_R130167", "SAMA_R46204")
list_span$"spaldingi-NE, spaldingi S" <- c("SAMA_R55731", "CCM_0044")
list_span$"spaldingi-CY" <- c("ANWC_R05242", "ANWC_R05240")

## Find nodes to add bars:
node_list <- list_span
for (clade in names(node_list)) { 
  node_list[[clade]] <- findMRCA(tree = tree, type = "node", tips = node_list[[clade]]) 
}

names(node_list) <- paste0("italic(\"", names(node_list), "\")") ## To allow spaces in ggtree.

## PART 5: A few more tip label changes ----

## Remove spaces:
sample_info$label <- gsub(x = sample_info$label, pattern = "_", replacement = " ")
sample_info$label <- paste0("plain(\"", sample_info$label, "\")") ## To allow spaces in ggtree.

## Change tree tips:
tree$tip.label <- as.character(sample_info$label)

## PART 5: Setting up tip symbols ----

## List of samples in each taxon and group tips by taxon to color tips in tree:
t_list <- split(x = sample_info, f = sample_info$CANDIDATE_SP_III)
t_list <- purrr::map(t_list, .f = function(x) as.character(x$label))
tree <- groupOTU(tree, t_list) ## Grouping.

## Colors and symbols in tree:
palette_df <- data.frame(CANDIDATE_SP = c("burbidgei", "inornatus-E", "inornatus-N", "inornatus-S", "lateralis", 
                                          "rimacola", "spaldingi-CY", "spaldingi-NE", "spaldingi-NW", "spaldingi-S", "spaldingi-TE",
                                          "superciliaris-E (N)", "superciliaris-E (S)", "superciliaris-K", "superciliaris-W", "x"),
                         symbol = c(23, 25, 21, 21, 21, 
                                    22, 24, 21, 21, 21, 23, 
                                    23, 22, 21, 21, NA),
                         color = c(rep("black", 15), NA), 
                         fill = c("gray50", "gray50", "#e7969c", "#7b4173", "#843c39", ## burbidgei, eutaenius, helenae, inornatus, lateralis.
                                  "white", "gray50", "#e39c39", "#3e2331", "#c9622f", "white", ## rimacola, robustus, spaldingi.
                                  "#66bba3", "#66bba3", "#162c3b", "#b5cf6b", ## superciliaris.
                                  NA)) ## x.
palette_df <- palette_df[palette_df$CANDIDATE_SP %in% names(t_list), ]

## PART 6: Plot tree ----

## Plot tree:
plot1 <- ggtree(tree, size = 0.5, ladderize = FALSE, color = "black") + ## Size = branch line thickness.
  
  ## Editing tree tips:
  geom_tiplab(color = "black", size = 1.55, offset = 0.005, parse = TRUE) +
  geom_tippoint(aes(fill = group, color = group, shape = group), size = 1.75, stroke = 0.5, alpha = 1) +
  
  ## Symbols:
  scale_color_manual(values = palette_df$color, name = "OTU or taxon") +
  scale_fill_manual(values = palette_df$fill, name = "OTU or taxon") +
  scale_shape_manual(values = palette_df$symbol, name = "OTU or taxon") +
  
  ## Add clade bars:
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[1]], label = names(node_list)[[1]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[2]], label = names(node_list)[[2]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[3]], label = names(node_list)[[3]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[4]], label = names(node_list)[[4]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[5]], label = names(node_list)[[5]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[6]], label = names(node_list)[[6]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[7]], label = names(node_list)[[7]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[8]], label = names(node_list)[[8]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[9]], label = names(node_list)[[9]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[10]], label = names(node_list)[[10]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[11]], label = names(node_list)[[11]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[12]], label = names(node_list)[[12]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[13]], label = names(node_list)[[13]], parse = TRUE, color = "black") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 3.5, offset.text = 0.01, offset = 0.07, node = node_list[[14]], label = names(node_list)[[14]], parse = TRUE, color = "black") +
  
  ## Legend:
  #theme(legend.position = "left")
  theme(legend.position = "none")
  #theme(legend.position = c(0.078, 0.83))
  
## Plot:
plot2 <- plot1 + xlim(0, 0.95)

## Saving plot:
ggsave(plot = plot2, filename = "cytb/Fig7_cytb.pdf", width = 8.5, height = 11, units = "in", limitsize = FALSE, dpi = 500)
ggsave(plot = plot2, filename = "cytb/Fig7_cytb.jpg", width = 8.5, height = 11, units = "in", limitsize = FALSE, dpi = 500)

## PART 7: Plot complete mitochondrial tree for supplementary ----

## If needed, clearing working space:
rm(list = ls())

## Read tree:
tree <- read.tree(file = "cytb/RAxML_bipartitions.n485.tre")
tree <- ladderize(tree, right = FALSE)

## Exclude outgroups:
exclude <- c("SAMR_37942_Ct_atla", "SAMR_57397_Ct_atla", "UMMZ_242607_Ct_aust", "UMMZ_242606_Ct_aust", "NTMR_22191_Ct_essi",
             "NTMR_22188_Ct_essi", "SAMR_46944_Ct_leon", "WAMR_139415_Ct_nigr", "WAMR_139414_Ct_nigr", "UMMZ_242633_Ct_pant",
             "UMMZ_242639_Ct_pant", "CUMV_14681_Ct_scho", "CUMV_14700_Ct_scho", "SAMAR_33571_Ct_taen", "SAMAR_33727_Ct_taen")
tree <- drop.tip(tree, exclude)

## Add sample info:
tips <- data.frame(SAMPLE_ID = tree$tip.label)
sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
sample_info <- sample_info[sample_info$SAMPLE_ID %in% tree$tip.label, ]
sample_info <- sample_info[c("SAMPLE_ID", "VOUCHER", "REGO", "CANDIDATE_SP_III")]
sample_info <- merge(tips, sample_info, all.x = TRUE)

## Table S1:
Table_S2 <- read.csv(file = "supplementary_files/Table_S2_cytb.csv", header = TRUE)
names(Table_S2)[names(Table_S2) == "Voucher.or.Tissue.Sample"] <- "label"
names(Table_S2)[names(Table_S2) == "Sample.ID"] <- "SAMPLE_ID"
sample_info <- merge(sample_info, Table_S2, by = "SAMPLE_ID", all.x = TRUE)

## ADD burbidgei:
sample_info$CANDIDATE_SP_III[sample_info$Original.Taxon == "burbidgei"] <- "burbidgei"

## Order sample info by tree position:
sample_info$SAMPLE_ID <- factor(x = sample_info$SAMPLE_ID, levels = tree$tip.label)
sample_info <- arrange(sample_info, SAMPLE_ID)

## Change tree tips:
tree$tip.label <- as.character(sample_info$label)

## PART 8: Setting up tip symbols ----

## List of samples in each taxon and group tips by taxon to color tips in tree:
sample_info$CANDIDATE_SP_III[is.na(sample_info$CANDIDATE_SP_III)] <- "No nuclear data"
t_list <- split(x = sample_info, f = sample_info$CANDIDATE_SP_III)
t_list <- purrr::map(t_list, .f = function(x) as.character(x$label))
tree <- groupOTU(tree, t_list) ## Grouping.

## Colors and symbols in tree:
palette_df <- data.frame(CANDIDATE_SP = c("burbidgei", "inornatus-E", "inornatus-N", "inornatus-S", "lateralis", 
                                          "No nuclear data",
                                          "rimacola", "spaldingi-CY", "spaldingi-NE", "spaldingi-NW", "spaldingi-S", "spaldingi-TE",
                                          "superciliaris-E (N)", "superciliaris-E (S)", "superciliaris-K", "superciliaris-W"),
                         symbol = c(23, 25, 21, 21, 21, 
                                    NA,
                                    22, 24, 21, 21, 21, 23, 
                                    23, 22, 21, 21),
                         fill = c("gray50", "gray50", "#e7969c", "#7b4173", "#843c39", ## burbidgei, eutaenius, helenae, inornatus, lateralis.
                                  NA,
                                  "white", "gray50", "#e39c39", "#3e2331", "#c9622f", "white", ## rimacola, robustus, spaldingi.
                                  "#66bba3", "#66bba3", "#162c3b", "#b5cf6b") ## superciliaris.
                                  )

## PART 9: Plot tree ----

## Plot tree:
plot1 <- ggtree(tree, size = 0.5, ladderize = FALSE, color = "black") + ## Size = branch line thickness.
  
  ## Editing tree tips:
  geom_tiplab(color = "black", size = 1.58, offset = 0.005, parse = TRUE) +
  geom_tippoint(aes(fill = group, shape = group), color = "black", size = 1.75, stroke = 0.5, alpha = 1) +
  
  ## Symbols:
  scale_fill_manual(values = palette_df$fill, name = "Nuclear unit:") +
  scale_shape_manual(values = palette_df$symbol, name = "Nuclear unit:") +
  
  ## Setting some other elements:
  ggtitle(label = "Figure S3. Phylogenetic relationships based on all mitochondrial samples.") +
  theme(plot.title = element_text(size = 16, margin = margin(t = 0, r = 0, b = 20, l = 0)),
        plot.margin = margin(t = 30, r = 30, b = 30, l = 30))
  
## Saving plot:
ggsave(plot = plot1, filename = "supplementary_files/Fig_S3_cytb.pdf", width = 8.5, height = 25, units = "in", limitsize = FALSE, dpi = 500)

## End of script.
