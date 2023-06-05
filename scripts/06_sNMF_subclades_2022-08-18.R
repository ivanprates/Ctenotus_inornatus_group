###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To run genetic clustering analyses using sNMF;
### To make bar plots based on the ancestry coefficients.
### To make maps based on cluster assignments.
### To extract and organize cluster assignments for downstream analyses.

### PART 1: Getting ready ----

## Packages:
library(ape)
library(ggrepel)
library(ggtree)
library(LEA)
library(pals)
library(patchwork)
library(phytools)
library(raster)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyverse)

## If needed, clearing working space:
#rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/" ## Suindara.
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## If needed, deleting old files:
unlink(here("sNMF"), recursive = TRUE)

## Creating directories to save results:
dir.create(path = here("sNMF"))
dir.create(path = here("sNMF/data"))
dir.create(path = here("sNMF/plots"))
dir.create(path = here("sNMF/qmatrices"))

## PART 2: Writing some functions we'll need ----

## Function: Run genotypic clustering analyses using sNMF.  
run_sNMF <- function(a, cpu, maxK, minK, miss_ind, miss_SNP, rep) {
  
  ## Testing:
  #miss_ind <- 0.5 ; miss_SNP <- 0.4 ; a <- 100 ; minK <- 1 ; maxK <- 5 ; rep <- 3 ; cpu = 8
  #dataset <- "inornatus_complex_n103_m50"

  ## Genetic data:
  gendata <- read.table(file = here(paste0("ipyrad/", dataset, "_outfiles/mi", miss_ind, "_ms", miss_SNP, ".usnps")), header = TRUE)
  
  ## Remove samples:
  remove_IDs <- c("WAMR_126010_Ct_rima", "WAMR_126015_Ct_rima", ## Too few samples for rimacola.
                  "SAMR_55874_Ct_euta", "SAMAR_55799_Ct_late", "NA_CCM0090_Ct_late", ## Too few samples for eutaenius.
                  "NTMR_17738_Ct_cogg", "NTMR_20378_Ct_robu", "NTMR_22166_Ct_robu", ## C. coggeri?
                  "AMSR_111493_Ct_spal", "AMSR_111494_Ct_spal") ## A species from Cape York?
                  #"SAMR_55881_Ct_robu") ## Probably different from spaldingi. But not in this assembly anyways.
                  #"QM_82086_Ct_spal") ## A species from Queensland?
  gendata <- gendata[gendata$SAMPLE_ID %in% setdiff(gendata$SAMPLE_ID, remove_IDs), ]
  
  ## Writing down data in geno format outside of R as requested by the LEA package:
  geno_df <- gendata[setdiff(names(gendata), "SAMPLE_ID")]
  write.geno(geno_df, output.file = here(paste0("sNMF/data/", dataset, "_mi", miss_ind, "_ms", miss_SNP, ".geno")))
  
  ## Running sNMF:
  project.snmf <- snmf(input.file = here(paste0("sNMF/data/", dataset, "_mi", miss_ind, "_ms", miss_SNP, ".geno")), 
                       entropy = TRUE, ploidy = 2, project = "new", seed = 2022,
                       CPU = cpu, K = minK:maxK, alpha = a, repetitions = rep)
    
  ## Showing summary of project results:
  snmf_summary <- summary(project.snmf)
  snmf_summary
  
  ## Selecting criterion to determine the best K:
  #crossEntropy <- snmf_summary$crossEntropy[1,] # Using min cross entropy across runs.
  crossEntropy <- snmf_summary$crossEntropy[2,] # Using mean cross entropy across runs.
  names(crossEntropy) <- gsub(x = names(crossEntropy), pattern = "K = ", replacement = "")
  
  ## Formatting:
  crossEntropy_df <- data.frame(K = names(crossEntropy), ce = crossEntropy)
  crossEntropy_df <- arrange(crossEntropy_df, by = ce)
  
  ## Selecting best K based on minimum (or mean) cross entropy among runs:
  K <- as.numeric(as.character(crossEntropy_df$K[which.min(crossEntropy_df$ce)]))
 
  ## Plotting mean entropy scores:
  png(filename = here(paste0("sNMF/plots/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, "_K", K, ".png")))
  plot(y = crossEntropy, x = names(crossEntropy), cex = 1.2, col = "blue", pch = 19, xlab = paste0("K values"), ylab = "Cross-entropy")
  lines(y = crossEntropy, x = names(crossEntropy), col = "blue")
  dev.off() ## Closing plot.
  
  ## Getting the best run under the best K based on cross-entropy scores:
  bestrun <- which.min(cross.entropy(project.snmf, K = K))
  
  ## Getting the Q matrix for the best run:
  qmatrix <- as.data.frame(Q(project.snmf, K = K, run = bestrun))
      
  ## Replace column names:
  names(qmatrix) <- gsub(names(qmatrix), pattern = "V", replacement = "cluster_\\1")
      
  ## Adding individual IDs:
  qmatrix$SAMPLE_ID <- gendata$SAMPLE_ID
    
  ## "Melt" dataframe using to assign samples to clusters based on max qscores:
  qmatrix_melt <- gather(qmatrix, key = cluster_assigned, value = coeff, 1:all_of(K))
    
  ## Assign specimens to cluster based on the highest qscore (coeff) values:
  cluster_assigned <- qmatrix_melt %>% group_by(SAMPLE_ID) %>% top_n(n = 1, wt = coeff)
  qmatrix <- merge(qmatrix, cluster_assigned, by = "SAMPLE_ID")
  
  ## Now, let's add sample information (including lat-longs) to the qmatrix:
  sample_info <- read.csv(file = here("sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv"), header = TRUE)
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% qmatrix$SAMPLE_ID, ]
  qmatrix <- merge(qmatrix, sample_info, by = "SAMPLE_ID")
  
  ## Using a phylogeny to order samples in structure plots:
  tree <- read.tree("RAxML/n242_m70_bs500/RAxML_bipartitions.n242_m70_bs500_7")
  tree <- ladderize(tree, right = FALSE)
  
  ## Obtaining the correct (i.e., rotated, ladderized) order of tips:
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  phylo_order <- tree$tip.label[ordered_tips]

  ## Using tip labels to guide sample order in plots:
  qmatrix$SAMPLE_ID <- factor(qmatrix$SAMPLE_ID, levels = phylo_order)
  qmatrix <- arrange(qmatrix, SAMPLE_ID)
  qmatrix$phylo_order <- 1:nrow(qmatrix)

  ## Save qmatrix:
  qmatrix <- qmatrix[setdiff(names(qmatrix), c("REGION", "VoucherNotes"))] ## Remove some columns.
  for (column in paste0("cluster_", 1:K)) { qmatrix[column] <- round(qmatrix[column], 2) } ## Round qscores.
  write.csv(qmatrix, file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), row.names = FALSE)

  ## Remove sNMF project:
  remove.snmfProject(file = here(paste0("sNMF/data/", dataset, "_mi", miss_ind, "_ms", miss_SNP, ".snmfProject")))
  
  ## Return:
  return(length(qmatrix$SAMPLE_ID))
  
} ## End of function.

## Function: Create color palette for plots.
get_palette <- function(dataset) {
  if (dataset == "spaldingi_complex_n72_m50") { palette <- c("#3e2331", "#e39c39", "#c9622f") }
  if (dataset == "superciliaris_complex_n46_m50") { palette <- c("#162c3b", "#66bba3", "#B5CF6B") }
  if (dataset == "inornatus_complex_n103_m50") { palette <- c("#843C39", "#E7969C", "#7B4173") }
  return(palette)
}

## Function: Plotting ancestry coefficients.
plot_sNMF_bars <- function(a, miss_ind, miss_SNP) {

  ## Testing:
  #miss_ind <- 0.7 ; miss_SNP <- 0.5 ; a <- 100
  
  ## Print status:
  print(paste0("Now making barplot!"))
   
  ## Read qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
  qmatrix <- qmatrix[!is.na(qmatrix$SAMPLE_ID), ]
  
  ## Arrange:
  qmatrix <- arrange(qmatrix, phylo_order)

  ## Exclude admixed samples to set factor levels:
  qlevels <- qmatrix[qmatrix$coeff > 0.95, ]

  ## Set factors levels to preserve order in plots:
  qmatrix$SAMPLE_ID <- factor(qmatrix$SAMPLE_ID, levels = qmatrix$SAMPLE_ID)
    
  ## Number of clusters:
  K <- length(unique(qmatrix$cluster_assigned))
  
  ## "Melt":
  qmatrix_m <- gather(qmatrix, key = sNMF_cluster, value = qscores, 2:(K+1))
        
  ## Set factors levels to preserve order in plots:
  qmatrix_m$sNMF_cluster <- factor(qmatrix_m$sNMF_cluster, levels = unique(qlevels$cluster_assigned))
  
  ## Color palette:
  palette <- get_palette(dataset)

  ## Creating stacked bar plot of ancestry coefficients:
  plot <- ggplot(data = qmatrix_m, aes(y = SAMPLE_ID)) +
      
    ## Adding bars that represent ancestry coefficients:
    geom_bar(aes(x = qscores, fill = sNMF_cluster), 
             color = "white", ## Color of bar lines.
             size = 0.01, ## Thickness of bar lines.
             stat = "identity", position = "fill",
	           show.legend = FALSE) +
        
    ## Filling bars by cluster assigment:
    scale_fill_manual(values = palette) + 
          
    ## Adjusting labels:
    labs(x = "Ancestry proportions", y = "") +
          
    ## Adjusting limits of the x axis:
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = round(c(0, 1), 0)) +
            
    ## Changing theme:
    theme_minimal() +
  
    ## Theme parameters:
    theme(#axis.text.y = element_text(color = "gray30", angle = 0, vjust = 0.5, hjust = 1, size = 3),
          axis.text.y = element_blank(), ## Removing IDs from plot.
          axis.title.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid = element_blank(), 
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 0, margin = margin(0, 0, 0, 0)),
          plot.margin = margin(r = 15, l = 10, t = 0, b = 0),
          plot.title = element_blank())
	
  ## Return:
  plot
  return(plot)
  
} ## End of function.

## Function: Plotting maps with the resulting cluster assignments.
plot_maps <- function(a, miss_ind, miss_SNP) {
  
  ## Testing:
  #miss_ind <- 0.5 ; miss_SNP <- 0.4 ; 
  #dataset <- "rima-robu-spal_n72_m50" ; a <- 3500
  #dataset <- "euta-inor-hele-late_n103_m50" ; a <- 1500
  #dataset <- "mast-supe_n46_m50" ; a <- 1500
  
  ## Print status:
  print(paste0("Now making maps!"))
  
  ## Read qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
  qmatrix <- qmatrix[!is.na(qmatrix$SAMPLE_ID), ]
  
  ## Edit labels to plot:
  qmatrix$UPDATED_SP <- gsub(qmatrix$UPDATED_SP, pattern = "cf. ", replacement = "")
  qmatrix$label <- gsub(qmatrix$UPDATED_SP, pattern = "(^[a-z]{4}).+", replacement = "\\1")

  ## Arrange:
  qmatrix <- arrange(qmatrix, desc(phylo_order))

  ## Exclude admixed samples to set factor levels:
  qlevels <- qmatrix[qmatrix$coeff > 0.95, ]

  ## Set order of factors to preserve sample order in plots:
  qmatrix$cluster_assigned <- factor(qmatrix$cluster_assigned, levels = unique(qlevels$cluster_assigned))
      
  ## Number of clusters:
  K <- length(unique(qmatrix$cluster_assigned))
  
  ## Change order of columns to keep cluster color order in map plots:
  qmatrix[c(paste0("cluster_", 1:K))] <- qmatrix[c(paste0("cluster_", 1:K))][levels(qmatrix$cluster_assigned)]
  
  ## Make qscores numeric:
  qmatrix[c(paste0("cluster_", 1:K))] <- apply(qmatrix[c(paste0("cluster_", 1:K))], 2, as.numeric)
  
  ## Creating labels for map facets:
  label_df <- data.frame(cluster = qmatrix$cluster_assigned, taxon = qmatrix$UPDATED_SP)
  label_df$taxon <- paste0("C. ",  label_df$taxon) 
  table_lb <- data.frame(table(label_df))
  table_lb <- arrange(table_lb, cluster)
  label_df <- table_lb %>% group_by(cluster) %>% top_n(Freq, n = 1)
  label_df <- label_df %>% group_by(cluster) %>% sample_n(size = 1)
  
  ## Create labels for duplicated taxon name across clusters:
  label_df$label <- as.character(label_df$taxon)
  duplicated <- names(table(label_df$taxon))[table(label_df$taxon) >= 2]
  for (taxon in duplicated) {
    n_taxon <- table(label_df$taxon)[taxon]
    new_labels <- paste0(taxon, " (", 1:n_taxon, ")")
    for (plabel in 1:length(new_labels)) { 
      label_df$label[label_df$taxon == taxon][plabel] <- as.character(new_labels[plabel]) 
  }} ## Close loops.
  
  ## Add labels to qmatrix:
  qmatrix$label <- factor(qmatrix$cluster_assigned, levels = label_df$cluster, labels = label_df$label)
  
  ## Adding samples not included in sNMF:
  sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
  sample_info <- sample_info[sample_info$CANDIDATE_SP_I %in% qmatrix$CANDIDATE_SP_I, ]
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% setdiff(sample_info$SAMPLE_ID, qmatrix$SAMPLE_ID), ]
  not_in_assembly <- c("SAMR_55881_Ct_robu", "SAMR_55891_Ct_robu", "SAMR_36875_Ct_regi", "SAMR_51089_Ct_saxa", "WAMR_163003_Ct_duri")
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% setdiff(sample_info$SAMPLE_ID, not_in_assembly), ]
  for (sample in sample_info$SAMPLE_ID) {
    #sample <- sample_info$SAMPLE_ID[1]
    sp <- sample_info$CANDIDATE_SP_I[sample_info$SAMPLE_ID == sample]
    label <- as.character(unique(qmatrix$label[qmatrix$CANDIDATE_SP_I == sp]))
    cluster <- as.character(unique(qmatrix$cluster_assigned[qmatrix$CANDIDATE_SP_I == sp]))
    sample_info$label[sample_info$SAMPLE_ID == sample] <- label
    sample_info$cluster_assigned[sample_info$SAMPLE_ID == sample] <- cluster }
    
  ## Combine:
  map_df <- rbind(qmatrix[c("SAMPLE_ID", "cluster_assigned", "LAT", "LON", "label")], sample_info[c("SAMPLE_ID", "cluster_assigned", "LAT", "LON", "label")])
  map_df$cluster_assigned <- factor(map_df$cluster_assigned, levels = unique(qlevels$cluster_assigned))
  
  ## Color palette:
  palette <- get_palette(dataset)
  
  ## Australia map:
  AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
  AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.
  
  ## Plot map with ggplot:
  sNMF_map <- ggplot() +
        
    ## Adding baseline map:
    geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray70", size = 0.15) +
        
    ## Adding ID labels on top of map:
    #geom_text_repel(data = qmatrix, aes(x = LON, y = LAT, label = label), color = "black", segment.color = "black", size = 3, segment.size = 0.25) +
        
    ## Adding lat-longs for sampled individuals:
    geom_point(data = map_df, aes(x = LON, y = LAT, fill = cluster_assigned), 
              shape = 21, alpha = 0.7, size = 2, color = "black", stroke = 0.2) +
     
    ## Setting a bunch of aesthetic parameters:
    theme_void() +
    guides(fill = "none", size = "none", alpha = "none", shape = "none") +
    scale_fill_manual(values = rev(palette)) +
        
    ## Setting up facets by cluster:
    facet_wrap(~label, ncol = 1) +
                   
    ## Setting some other elements:
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
          strip.text = element_text(size = 9, face = "italic", margin = margin(t = 0, r = 0, b = 2, l = 0)),
          #strip.text = element_blank(),
          strip.background = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank())
          
  ## Return:
  return(sNMF_map)
    
} ## End of function.

## Function: Plot trees.
plot_tree <- function(a, miss_ind, miss_SNP) {

  ## Testing:
  #miss_ind <- 0.5 ; miss_SNP <- 0.4 ; a <- 500
  #dataset <- "inornatus_complex_n103_m50"
  
  ## Print status:
  print("Now making the trees!")
  
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
   
  ## Read qmatrix and keep only samples from taxa used in sNMF analyses:
  qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
  MRCA_node <- findMRCA(tree, tips = as.character(qmatrix$SAMPLE_ID))
  keep_tips <- getDescendants(tree, node = MRCA_node)
  keep_tips <- keep_tips[keep_tips <= length(tree$tip.label)]
  keep_tips <- tree$tip.label[keep_tips]
  remove_tips <- tree$tip.label[!tree$tip.label %in% keep_tips]
  tree <- drop.tip(tree, tree$tip.label[match(remove_tips, tree$tip.label)])
  
  ## Tree labels:
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% tree$tip.label, ]
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
  
  ## Fixing missIDs:
  sample_info$COLLECTOR_fixed <- sample_info$COLLECTOR_SP 
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NMVD_67793_Ct_leae"] <- "saxatilis" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "SAMAR_34180_Ct_bore"] <- "saxatilis" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "SAMR_34123_Ct_hill"] <- "inornatus" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM4920_Ct_essi"] <- "inornatus" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM2518_Ct_vert"] <- "inornatus" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM2565_Ct_vert"] <- "inornatus" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM1061_Ct_sp"] <- "inornatus" ## In tree (as WAMR_174673).
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM0044_Ct_euta"] <- "robustus" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NTMR_17738_Ct_cogg"] <- "robustus" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM5477_Ct_vert"] <- "spaldingi" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM0090_Ct_late"] <- "eutaenius" ## In tree (as QMJ_94601).
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "SAMAR_55799_Ct_late"] <- "eutaenius" ## In tree.
  sample_info$COLLECTOR_fixed[sample_info$SAMPLE_ID == "NA_CCM1530_Ct_robu"] <- "inornatus" ## In tree (as WAMR_174689).

  ## Order by tree position and change tip labels:
  sample_info$SAMPLE_ID <- factor(x = sample_info$SAMPLE_ID, levels = tree$tip.label)
  sample_info <- arrange(sample_info, SAMPLE_ID)
  
  ## Sample labels:
  tree$tip.label <- as.character(sample_info$label)
  
  ## List of samples in each taxon:
  t_list <- split(x = sample_info, f = sample_info$COLLECTOR_fixed)
  t_list <- purrr::map(t_list, .f = function(x) as.character(x$label))
  
  ## Group tips by taxon to edit tips in ggtree:
  tree <- groupOTU(tree, t_list)
  
  ## Color palette and shapes:
  symbols_df <- data.frame(taxa = c("borealis", "brachyonyx", "eutaenius", "fallens", "helenae", "inornatus",
                                    "lateralis", "mastigura", "rimacola", "robustus", "saxatilis", "severus", "spaldingi"),
                           shapes = c(21, 25, 24, 25, 24, 23, 24, 23, 21, 21, 23, 25, 21),
                           palette = c("gray50", "gray50", "black", "white", "gray50", "black", 
                                       "white", "gray50", "gray70", "white", "white", "black", "black")) 
  symbols_df <- symbols_df[symbols_df$taxa %in% unique(sample_info$COLLECTOR_fixed), ]

  ## Plot tree:
  plot <- ggtree(tree, size = 0.5, ladderize = FALSE, color = "black") + # aes(color = group)) + ## Size = branch line thickness.
    
    ## Editing tree tips:
    geom_tiplab(color = "black", size = 1.5, offset = 0.0002, parse = TRUE) +
    
    ## Tip points:
    geom_tippoint(aes(shape = group, fill = group), color = "black", size = 1.2, stroke = 0.3) +
    
    ## Editing tree tips:
    scale_fill_manual(values =  symbols_df$palette, name = "Original taxon ID") +
    scale_shape_manual(values =  symbols_df$shapes, name = "Original taxon ID") +
    
    ## Other edits:
    theme(legend.position = "left",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_blank(),
          plot.caption = element_text(size = 16, hjust = 0.5, margin = margin(t = 0, b = 0, l = 0, r = 0))) 
  
  ## Plot:
  if (dataset == "spaldingi_complex_n72_m50") { phyx <- 0.0105 }
  if (dataset == "superciliaris_complex_n46_m50") { phyx <- 0.0099 }
  if (dataset == "inornatus_complex_n103_m50") { phyx <- 0.0105 }
  plot2 <- plot + xlim(0, phyx)
  
  ## Flip nodes:
  if (dataset == "inornatus_complex_n103_m50") { 
  #plot2 + geom_text(aes(label = node)) ## Check nodes.
  plot3 <- ggtree::flip(plot2, 165, 183) %>% ggtree::flip(167, 175) 
  plot2 <- plot3 }
  
  ## Return:
  return(plot2)

} ## End of function.

## Function: Run sNMF and make plot under combinations of alpha and missing data.
run_all <- function(miss_ind, miss_SNP) {

  ## Set alpha (based on preliminary runs) and run sNMF:
  a <- 500
  run_sNMF(cpu = 8, minK = 1, maxK = 8, rep = 20, a = a, miss_ind = miss_ind, miss_SNP = miss_SNP) 
  
  ## Create plots post-analysis:
  plot_t <- plot_tree(a = a, miss_ind = miss_ind, miss_SNP = miss_SNP) 
  plot_b <- plot_sNMF_bars(a = a, miss_ind = miss_ind, miss_SNP = miss_SNP)
  plot_m <- plot_maps(a = a, miss_ind = miss_ind, miss_SNP = miss_SNP)
  
  ## Combining plots:
  plot_f <- ( plot_t + plot_b + plot_m ) + plot_layout(widths = c(1.25, 0.75, 1))
  
  ## To get the number of samples per run, read qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
  nsamples <- nrow(qmatrix)
  
  ## Saving plot:
  for (format in c("jpg", "pdf")) {
  ggsave(plot = plot_f, width = 20, height = 15, units = "cm", limitsize = FALSE, device = format, dpi = 300,
         filename = here(paste0("sNMF/plots/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, "_n", nsamples, ".", format))) }

} ## End of function.

## PART 3: Use functions we wrote to run clustering analyses and make plots ----

## Genetic dataset we'll use (ipyrad output):
dataset1 <- "spaldingi_complex_n72_m50"
dataset2 <- "superciliaris_complex_n46_m50"
dataset3 <- "inornatus_complex_n103_m50"
datasets <- c(dataset1, dataset2, dataset3)
#datasets <- dataset3

## Run sNMF and get plots:
for (dataset in datasets) { 
  run_all(miss_ind = 0.5, miss_SNP = 0.4) 
}

# ## PART 4: Checking samples not used in sNMF analyses ----
# 
# ## sNMF parameters:
# miss_ind <- 0.5 ; miss_SNP <- 0.4
# 
# dataset <- "rima-robu-spal_n72_m50" ; a <- 3500
# qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
# 
# dataset <- "euta-inor-hele-late_n103_m50" ; a <- 1500 
# qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
# 
# dataset <- "mast-supe_n46_m50" ; a <- 1500
# qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
# 
# tree <- read.tree(here("RAxML/n242_m70_bs500/RAxML_bipartitions.n242_m70_bs500_7"))
# MRCA_node <- findMRCA(tree, tips = as.character(qmatrix$SAMPLE_ID))
# keep_tips <- getDescendants(tree, node = MRCA_node)
# keep_tips <- keep_tips[keep_tips <= length(tree$tip.label)]
# keep_tips <- tree$tip.label[keep_tips]
# remove_tips <- tree$tip.label[!tree$tip.label %in% keep_tips]
# tree <- drop.tip(tree, tree$tip.label[match(remove_tips, tree$tip.label)])
# excluded <- sort(setdiff(tree$tip.label, tree$tip.label[tree$tip.label %in% qmatrix$SAMPLE_ID]))
# sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
# sample_info <- sample_info[sample_info$SAMPLE_ID %in% excluded, ]
# sample_info <- sample_info[sample_info$UPDATED_SP %in% setdiff(unique(sample_info$UPDATED_SP), c("rimacola", "eutaenius")), ]
# 
# ## Plot map:
# AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
# AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.
#   
# ## Plot map with ggplot:
# excluded_map <- ggplot() +
#   geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray70", size = 1) +
#   geom_point(data = sample_info , aes(x = LON, y = LAT, fill = CANDIDATE_SP_I), shape = 21, size = 10, color = "black", stroke = 0.2) +
#   geom_label(data = sample_info , aes(x = LON, y = LAT, label = SAMPLE_ID), size = 3) +
#   theme_void()  
# excluded_map

## End of script.
