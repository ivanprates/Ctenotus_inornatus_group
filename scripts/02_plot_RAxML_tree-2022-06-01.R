###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To plot a RAxML tree based on the ddRAD data.

## PART 1: Setting up ----

## Packages:
library(ape)
library(ggtree)
library(pals)
library(patchwork)
library(phytools)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/" ## Suindara.
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## PART 2: Trees and tip labels ----

## Read tree:
tree <- read.tree(here("RAxML/n242_m70_bs500/RAxML_bipartitions.n242_m70_bs500_7"))
tree <- ladderize(tree, right = FALSE)
  
## Sample info:
sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
sample_info <- sample_info[sample_info$SAMPLE_ID %in% tree$tip.label, ]
sample_info$CANDIDATE_SP_I <- gsub(sample_info$CANDIDATE_SP_I, pattern = "cf. ", replacement = "")

## Exclude outgroups and mixed-up sample:
exclude <- c("atlas", "australis", "bipes", "essingtonii", "ips", "leonhardii", "nigrilineatus",
             "pantherinus", "schomburgkii", "taeniolatus") ## Outgroups.
sample_info <- sample_info[sample_info$CANDIDATE_SP_I %in% setdiff(sample_info$CANDIDATE_SP_I, exclude), ]
sample_info <- sample_info[sample_info$SAMPLE_ID %in% setdiff(sample_info$SAMPLE_ID, "WAMR_163003_Ct_duri"), ]
tree <- keep.tip(tree, tree$tip.label[match(sample_info$SAMPLE_ID, tree$tip.label)])

## Tree labels:
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

## Order by tree position:
sample_info$SAMPLE_ID <- factor(x = sample_info$SAMPLE_ID, levels = tree$tip.label)
sample_info <- arrange(sample_info, SAMPLE_ID)

## Change tree tips:
tree$tip.label <- as.character(sample_info$label)

## Define clades to add bars:
list_span <- list()
list_span$robustus <- c("WAM_R174671", "NTM_R17738")
list_span$rimacola <- c("WAM_R126010", "WAM_R126015")
list_span$spaldingi <- c("QM_82086", "SAMA_R52303")
list_span$"cf. mastigura" <- c("WAM_R174676", "CCM_0850")
list_span$superciliaris <- c("NTM_R22432", "WAM_R158376")
list_span$"cf. eutaenius" <- c("SAMA_R55874", "SAMA_R55799")
list_span$lateralis <- c("SAMA_R54463", "SAMA_R42768")
list_span$inornatus <- c("CCM_3739", "WAM_R141131") 

## Find nodes to add bars:
node_list <- list_span
for (clade in names(node_list)) { node_list[[clade]] <- findMRCA(tree = tree, type = "node", tips = node_list[[clade]]) }

## Create group labels:
names(node_list) <- paste0("italic(\"C. ", names(node_list), "\")")
names(node_list) <- paste0("paste(", names(node_list), ", )")

## Remove spaces from tip labels:
sample_info$label <- gsub(x = sample_info$label, pattern = "_", replacement = " ")
sample_info$label <- paste0("plain(\"", sample_info$label, "\")") ## To allow spaces in ggtree.
tree$tip.label <- as.character(sample_info$label)

## Plot tree:
plot1 <- ggtree(tree, size = 0.5, color = "black", ladderize = FALSE) + # aes(color = group)

  ## Editing tree tips:
  #geom_tiplab(color = "black", size = 1.5, offset = 0.0001, parse = TRUE) +
  
  ## Node symbols and text:
  geom_point2(aes(subset = (node == 221)), shape = 23, size = 4, fill = "red") +
  geom_point2(aes(subset = (node == 292)), shape = 23, size = 4, fill = "red") +
  geom_point2(aes(subset = (node == 337)), shape = 23, size = 4, fill = "red") +
  geom_text(aes(x = 0.0035, y = 14), label = "spaldingi", check_overlap = TRUE, color = "red", size = 5, fontface = "italic") +
  geom_text(aes(x = 0.0035, y = 9), label = "complex", check_overlap = TRUE, color = "red", size = 5) +
  geom_text(aes(x = 0.0048, y = 80), label = "superciliaris", check_overlap = TRUE, color = "red", size = 5, fontface = "italic") +
  geom_text(aes(x = 0.0048, y = 75), label = "complex", check_overlap = TRUE, color = "red", size = 5) +
  geom_text(aes(x = 0.0048, y = 143), label = "inornatus", check_overlap = TRUE, color = "red", size = 5, fontface = "italic") +
  geom_text(aes(x = 0.0048, y = 138), label = "complex", check_overlap = TRUE, color = "red", size = 5) +
  
  ## Other edits:
  theme(legend.position = "none",
        plot.title = element_blank(),
        plot.caption = element_text(size = 16, hjust = 0.5, margin = margin(t = 0, b = 0, l = 0, r = 0))) +
  
  ## Add clade bars:
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[1]], label = names(node_list)[[1]], parse = TRUE, color = "gray20") + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[2]], label = names(node_list)[[2]], parse = TRUE, color = "gray20") + 
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[3]], label = names(node_list)[[3]], parse = TRUE, color = "gray20") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[4]], label = names(node_list)[[4]], parse = TRUE, color = "gray20") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[5]], label = names(node_list)[[5]], parse = TRUE, color = "gray20") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[6]], label = names(node_list)[[6]], parse = TRUE, color = "gray20") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[7]], label = names(node_list)[[7]], parse = TRUE, color = "gray20") +
  geom_cladelabel(barsize = 1.5, align = TRUE, fontsize = 4.5, offset.text = 0.00025, offset = 0.0002, node = node_list[[8]], label = names(node_list)[[8]], parse = TRUE, color = "gray20")

## Plot:
plot2 <- plot1 + xlim(0, 0.0195)

## Flip nodes:
#plot2 + geom_text(aes(label = node)) ## Check nodes.
plot3 <- ggtree::flip(plot2, 399, 417) %>% ggtree::flip(401, 409) ## Flip nodes.
plot3

## Saving plot:
ggsave(plot = plot3, width = 7, height = 10, units = "in", limitsize = FALSE, dpi = 300, filename = here("RAxML/Fig1_RAxML_ddRAD.pdf"))
ggsave(plot = plot3, width = 7, height = 10, units = "in", limitsize = FALSE, dpi = 300, filename = here("RAxML/Fig1_RAxML_ddRAD.jpg"))

## End of script.
