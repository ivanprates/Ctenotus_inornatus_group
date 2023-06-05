###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### April 2023.
### The goals of this script are:
### To make maps based on sNMF cluster assignments using pie charts.

### PART 1: Getting ready ----

## Packages:
library(pals)
library(patchwork)
library(raster)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(sf)
library(tidyverse)

## If needed, clearing working space:
#rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/" ## Suindara.
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## Function: Create color palette for plots.
get_palette <- function(dataset) {
  if (dataset == "spaldingi_complex_n72_m50") { palette <- c("#3e2331ff", "#e39c39", "#c9622f") }
  if (dataset == "superciliaris_complex_n46_m50") { palette <- c("#162c3b", "#66bba3", "#B5CF6B") }
  if (dataset == "inornatus_complex_n103_m50") { palette <- c("#cd5c5cff", "#E7969C", "#7B4173") }
  return(palette)
}

## Function: Plotting maps with the resulting cluster assignments.
plot_maps <- function(dataset, miss_ind = 0.5, miss_SNP = 0.4, a = 500) {
  
  ## Testing:
  #dataset <- "inornatus_complex_n103_m50" 
  #dataset <- "spaldingi_complex_n72_m50"
  #dataset <- "superciliaris_complex_n46_m50"
  
  ## Print status:
  print(paste0("Now making maps!"))
  
  ## Read qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF/qmatrices/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, ".csv")), header = TRUE)
  
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
  
  ## Calculating mean qscores per locality:
  qmatrix_mean <- qmatrix %>% group_by(cluster_assigned, LAT, LON) %>% summarise_at(c(paste0("cluster_", 1:K)), mean) 
  
  ## Round qscores:
  for (column in paste0("cluster_", 1:K)) { qmatrix_mean[column] <- round(qmatrix_mean[column], 2) }
  
  ## Color palette:
  palette <- get_palette(dataset)
  
  ## Australia map:
  AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
  AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.
  
  ## Plot map with ggplot:
  sNMF_map <- ggplot() +
        
    ## Adding baseline map:
    geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray70", size = 0.15) +
        
    ## Pie chart of the ancestry coefficients: 
    geom_scatterpie(data = qmatrix_mean, aes(x = LON, y = LAT), 
                    pie_scale = 2, alpha = 1, color = "transparent",
                    cols = paste0("cluster_", 1:K)) +
    
    ## Setting a bunch of aesthetic parameters:
    theme_void() +
    guides(fill = "none", size = "none", alpha = "none", shape = "none") +
    scale_fill_manual(values = rev(palette)) +
        
    ## Setting up facets by cluster:
    facet_wrap(~cluster_assigned, ncol = 1) +
                   
    ## Setting some other elements:
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray30", linewidth = 0.5),
          strip.text = element_text(size = 9, face = "italic", margin = margin(t = 0, r = 0, b = 2, l = 0)),
          #strip.text = element_blank(),
          strip.background = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank())
          
  ## Saving plot:
  ggsave(plot = sNMF_map, width = 6.66, height = 15, units = "cm", limitsize = FALSE, device = pdf, dpi = 300,
         filename = here(paste0("sNMF/plots/", dataset, "_mi", miss_ind, "_ms", miss_SNP, "_a", a, "_sNMF_pie_maps.pdf")))

  ## Return:
  return(sNMF_map)
    
} ## End of function.

## Run:
plot_maps(dataset = "inornatus_complex_n103_m50") 
plot_maps(dataset = "spaldingi_complex_n72_m50")
plot_maps(dataset = "superciliaris_complex_n46_m50")

## End of script.
