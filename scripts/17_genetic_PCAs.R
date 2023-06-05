###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### April 2023.
### The goals of this script are:
### To perform genetic PCA based on the ddRAD data.
### To plot the results.

## PART 1: Getting ready ----

## Packages:
library(LEA)
library(patchwork)
library(reshape2)
library(scales)
library(tidyverse)

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/" ## Suindara.
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## Function: Create color palette for plots.
get_palette <- function(dataset) {
  if (dataset == "spaldingi_complex_n72_m50") { palette <- c("#3e2331ff", "#e39c39", "#c9622f") }
  if (dataset == "superciliaris_complex_n46_m50") { palette <- c("#162c3b", "#66bba3", "#B5CF6B") }
  if (dataset == "inornatus_complex_n103_m50") { palette <- c("#E7969C", "#7B4173", "#cd5c5cff") }
  return(palette)
}

## PART 2: Genetic PCA ----

## Function:
run_PCA <- function(dataset, title) {

  ## Testing:
  #dataset <- "spaldingi_complex_n72_m50"
  #dataset <- "inornatus_complex_n103_m50"
  #dataset <- "superciliaris_complex_n46_m50"
  
  ## Color palette:
  palette <- get_palette(dataset)
  
  ## Run PCA analysis. Command includes location of geno file:
  PCA <- pca(paste0("sNMF/data/", dataset, "_mi0.5_ms0.4.geno"))
  
  ## Displaying information on the PCA analysis:
  summary(PCA)
  
  ## Plotting eigenvalues:
  plot(PCA, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")
  
  ## Storing the proportion of variation represented by selected PC axes"
  PC_perc <- vector("list", 4) # We'll use the four first PC axes and save them into this list.
  for (p in c(1:4)) { PC_perc[[p]] <- round((summary(PCA)[2, p]*100), digits = 0) 
  names(PC_perc) <- c(paste0("PC", rep(1:4))) }
  
  ## Saving genetic PCA axes to plot later:
  pcadata <- as.data.frame(PCA$projections)
  pcadata <- pcadata[, 1:4]
  names(pcadata) <- c(paste0("PC", rep(1:4)))
  
  ## Add sample IDs:
  qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi0.5_ms0.4_a500.csv")), header = TRUE)
  qmatrix <- arrange(qmatrix, SAMPLE_ID)
  pcadata$SAMPLE_ID <- qmatrix$SAMPLE_ID
  pcadata$cluster_assigned <- qmatrix$cluster_assigned
  
  ## Add labels:
  si <- read.csv(file = "sample_information/ddRAD_sample_information_n1315_IP_2022-08-18.csv", header = TRUE)
  si <- si[si$SAMPLE_ID %in% pcadata$SAMPLE_ID, ]
  si$WORKING_TAXON <- gsub(si$WORKING_TAXON, pattern = "superciliaris_TE", replacement = "superciliaris_E")
  si$WORKING_TAXON <- gsub(si$WORKING_TAXON, pattern = "spaldingi_NW", replacement = "robustus_NW")
  for (sample in pcadata$SAMPLE_ID) {
    pcadata$WORKING_TAXON[pcadata$SAMPLE_ID == sample] <- si$WORKING_TAXON[si$SAMPLE_ID == sample]
  }
  
  ## Plotting:
  plot_PCA <- ggplot(data = pcadata, aes(x = PC1, y = PC2)) +
    
    ## Configuring point size and shape:
    geom_point(aes(fill = WORKING_TAXON), shape = 21,
               size = 3, color = "black", alpha = 1, show.legend = TRUE) +
    
    ## Adding ID labels on top of map:
    #geom_text_repel(aes(label = cluster_assigned), color = "black", segment.color = "black", size = 3, segment.size = 0.25) +
    
    ## Expand plot limits:
    expand_limits(y = c(min(pcadata$PC2)-2, max(pcadata$PC2)+2), x = c(min(pcadata$PC1)-2, max(pcadata$PC1)+2)) +
    
    ## Labels:
    labs(y = paste0("PC2 (", PC_perc$PC2, "%)"), 
         x = paste0("PC1 (", PC_perc$PC1, "%)"),
         title = title) +
    
    ## Setting color scheme:
    scale_fill_manual(values = palette) + 
    
    ## Also setting values on axes:
    scale_x_continuous(breaks = pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5)) +
    
    ## Other edits:
    theme_bw() +
    theme(plot.margin = margin(b = 10, r = 10, l = 10, t = 10),
          panel.border = element_rect(size = 1, colour = "gray20"),
          axis.text = element_text(color = "gray20", size = 12),
          axis.title = element_text(color = "gray20", size = 14),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
          
  ## Return:
  plot_PCA
  return(plot_PCA)
  
  ## Remove PCA-related files:
  unlink("*.pcaProject")
  unlink("sNMF/data/*.lfmm")
  unlink("sNMF/data/*.pca", recursive = TRUE)

}

## Run:
inor <- run_PCA(dataset = "inornatus_complex_n103_m50", title = "inornatus complex") 
spal <- run_PCA(dataset = "spaldingi_complex_n72_m50", title = "robustus complex")
supe <- run_PCA(dataset = "superciliaris_complex_n46_m50", title = "superciliaris complex")

## Combine:
pf <- spal / supe / inor
pf

## Saving plot:
ggsave(plot = pf, width = 15, height = 25, units = "cm", limitsize = FALSE, device = pdf, 
       dpi = 300, filename = "sNMF/plots/PCAs.pdf")

## End of script.
