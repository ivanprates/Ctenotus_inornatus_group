###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To estimate patterns of geographic isolation by distance (IBD) in focal subclades of the Ctenotus inornatus species group.
### To make scatter plots of the IBD relationships.
### To plot samples by candidate species on a map.
### Useful resources:
### https://popgen.nescent.org/2015-05-18-Dist-SNP.html
### https://www.rdocumentation.org/packages/hierfstat/versions/0.5-7/topics/genet.dist
### https://www.rdocumentation.org/packages/poppr/versions/2.8.5/topics/nei.dist
### https://pubmed.ncbi.nlm.nih.gov/20067366/ ## Study comparing metrics of genetic distance.

## PART 1: Getting ready ----

## Packages:
library(BEDASSLE)
library(fossil)
library(ggrepel)
library(here)
library(LEA)
library(patchwork)
library(reshape2)
library(raster)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/"
setwd(path)
detach("package:here", unload = TRUE)
library(here)

## Directories:
dir.create(here("IBD"))
dir.create(here("IBD/data"))
dir.create(here("IBD/plots"))

## Genetic dataset:
miss_SNP <- 0.6 ## Maximum missing data per SNP.
miss_ind <- 0.8 ## Maximum missing data per sample.
dataset <- "inornatusgr_R1_c90_n223_m95"
gendata <- read.table(paste0("ipyrad/", dataset, "_outfiles/mi", miss_ind, "_ms", miss_SNP, ".usnps"), sep = " ", header = TRUE)
snps <- dim(gendata)[2]

## PART 2: Writing some functions we'll need ----

## Function: Estimating genetic and geographic distances:
get_dist <- function(subclade) {

  ## Testing:
  #subclade <- "spaldingi_suppl"
  #subclade <- "inornatus_suppl"
  #subclade <- "superciliaris_suppl"

  ## Print status:
  print(paste0("Starting to process ", subclade, "!"))
  
  ## Candidate species:
  if (subclade == "inornatus_suppl") { taxa <- c("inornatus-N", "inornatus-S", "lateralis") }
  if (subclade == "inornatus") { taxa <- c("inornatus-N", "inornatus-S") }
  if (subclade == "superciliaris") { taxa <- c("superciliaris-E (S)", "superciliaris-E (N)", "superciliaris-W") }
  if (subclade == "superciliaris_suppl") { taxa <- c("superciliaris-E (S)", "superciliaris-E (N)", "superciliaris-W") }
  if (subclade == "spaldingi_suppl") { taxa <- c("spaldingi-CY", "spaldingi-NE", "spaldingi-S") } 
  if (subclade == "spaldingi") { taxa <- c("spaldingi-CY", "spaldingi-NE", "spaldingi-S") } 
  if (subclade == "robustus") { taxa <- c("spaldingi-NW", "spaldingi-TE") } 
  
  ## Loading sample information:
  sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
  sample_info <- sample_info[sample_info$CANDIDATE_SP_III %in% taxa, ]
  if (subclade == "spaldingi_suppl") { 
    sample_info$CANDIDATE_SP_III[sample_info$CANDIDATE_SP_III %in% c("spaldingi-NW", "spaldingi-TE")] <- "spaldingi-NW-TE" }
  if (subclade == "superciliaris_suppl") { 
    sample_info$CANDIDATE_SP_III[sample_info$CANDIDATE_SP_III %in% c("superciliaris-E (S)", "superciliaris-E (N)")] <- "superciliaris-EN-ES" }
  sample_info$LAT <- round(sample_info$LAT, 4)
  sample_info$LON <- round(sample_info$LON, 4)
  sample_info <- arrange(sample_info, by = SAMPLE_ID)
  
  ## Removing seemingly mixed-up samples and those corresponding to distinct, undersampled species:
  remove <- c("WAMR_163003_Ct_duri", "AMSR_111493_Ct_spal", "AMSR_111494_Ct_spal")
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% setdiff(sample_info$SAMPLE_ID, remove), ]
  
  ## Restricting genetic samples and sample info to focal subclade:
  gen_df <- gendata[gendata$SAMPLE_ID %in% sample_info$SAMPLE_ID, ]
  gen_df$SAMPLE_ID <- factor(gen_df$SAMPLE_ID, levels = gen_df$SAMPLE_ID)
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% gen_df$SAMPLE_ID, ]
  
  ## Estimating a geographic distance matrix:
  geo_dist <- earth.dist(sample_info[, c("LON", "LAT")], dist = FALSE)
  colnames(geo_dist) = sample_info$SAMPLE_ID
  rownames(geo_dist) = sample_info$SAMPLE_ID
  diag(geo_dist) <- NA ## Replacing diagonals (individual compared to itself) with NA.
  geo_dist[geo_dist == 0] <- NA ## Replacing zero distances (i.e., samples from same site) with NA.
  rownames(gen_df) <- gen_df$SAMPLE_ID ## Row names.
      
  ## Remove SAMPLE_ID column:
  SNP_df <- gen_df[setdiff(names(gen_df), "SAMPLE_ID")]
  dim(SNP_df)
  
  ## Now estimating Fst using BEDASSLE.
  ## Let's first prepare the data to the required format.
  ## Convert to matrix with numeric column values:
  allele_counts <- apply(SNP_df, 2, as.numeric) 
  
  ## Add row names again:
  row.names(allele_counts) <- row.names(SNP_df)
    
  ## Create a sample size object:
  sample_sizes <- allele_counts
    
  ## Change 9 (i.e., missing data) for zeros (i.e., locus not sampled) in the allele counts matrix:
  allele_counts[allele_counts == 9] <- 0     
    
  ## Removing invariable sites:
  SNPs_to_keep <- apply(allele_counts, 2, function(x) length(unique(x)) > 1)
  allele_counts <- allele_counts[, SNPs_to_keep]
    
  ## Now, in the sample size object, keep only SNPs present in the final allele count matrix:
  sample_sizes <- sample_sizes[, colnames(allele_counts)]
    
  ## Since organism is diploid, replace with total number of alleles per sample (= 2):
  sample_sizes[sample_sizes == 1 ] <- 2
  sample_sizes[sample_sizes == 0 ] <- 2
    
  ## Change NAs for zeros in the sample size matrix:
  sample_sizes[sample_sizes == 9] <- 0
    
  ## Estimate Fst using BEDASSLE:
  gen_dist <- calculate.all.pairwise.Fst(allele.counts = allele_counts, sample.sizes = sample_sizes)
  gen_dist <- as.matrix(gen_dist)
  
  ## Checkpoint: Is matrix symetric? Must be "TRUE":
  isSymmetric(gen_dist)
  
  ## Change column and row names, set diagonal as NA:
  colnames(gen_dist) <- row.names(SNP_df)
  rownames(gen_dist) <- row.names(SNP_df)
  diag(gen_dist) <- NA
  
  ## Melting distance matrices to a list of pairwise distances:
  geo_vector <- melt(geo_dist)
  gen_vector <- melt(gen_dist)
  
  ## Checkpoint: Are data paired across matrices? Must all be "TRUE".
  table(geo_vector$Var1 == gen_vector$Var1 & geo_vector$Var2 == gen_vector$Var2)
  
  ## Combine:
  dist_df <- data.frame(subclade = subclade, SAMPLE1 = gen_vector$Var1, SAMPLE2 = gen_vector$Var2, 
                        Genetic_distance = gen_vector$value, Geographic_distance = geo_vector$value)
  
  ## Adding candidate species information: 
  dist_df$CANDIDATE_SP1 <- sample_info$CANDIDATE_SP_III[dist_df$SAMPLE1]
  dist_df$CANDIDATE_SP2 <- sample_info$CANDIDATE_SP_III[dist_df$SAMPLE2]
  
  ## Establishing whether comparison is intra-group or inter-group:
  dist_df$comparison <- "intragroup"
  for (i in 1:nrow(dist_df)) {
    if (dist_df$CANDIDATE_SP1[i] != dist_df$CANDIDATE_SP2[i]) { dist_df$comparison[i] <- "intergroup" }
  } ## Close loop.
  
  ## Saving:
  write.csv(x = dist_df, row.names = FALSE, file = paste0("IBD/data/", subclade, ".csv"), quote = FALSE)
  
  ## Print status:
  print(paste0("Done processing ", subclade, "!"))
  
} ## End of function.
  
## Function: Plotting IBD ----

## Function:
plot_IBD <- function(subclade) {
  
  ## Testing:
  #subclade <- "spaldingi_suppl"
  #subclade <- "superciliaris_suppl"
  
  ## Import distances estimated previously:
  dist_df <- read.csv(file = paste0("IBD/data/", subclade, ".csv"), header = TRUE)
  dist_df <- dist_df[!is.na(dist_df$Genetic_distance), ]
  
  ## Plot:
  plot <- ggplot(data = dist_df) + 
      
    ## Add points and regression lines:
    geom_point(aes(x = Geographic_distance, y = Genetic_distance, fill = comparison), color = "black", size = 1.5, shape = 21,
               alpha = 0.5, stroke = 0.2) +
    #geom_smooth(data = filter(dist_df, comparison == "intragroup"), aes(x = Geographic_distance, y = Genetic_distance), method = "lm", se = FALSE, size = 0.5, color = "gray30") +
    #geom_smooth(data = filter(dist_df, comparison == "intergroup"), aes(x = Geographic_distance, y = Genetic_distance), method = "lm", se = FALSE, size = 0.5, color = "gray30") +
    
    ## Titles and axes:
    
    xlab("Geographic distances") +
    ylab("Genetic distances (Fst)") +
    ylim(0, 1) +
    
    ## Colors and legends:
    guides(fill = "none") +
    scale_fill_manual(values = c("black", "gray70")) + 
    #labs(title = titlelabel) +
    
    ## Other params:
    theme_bw() +
    theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 5),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          #plot.title = element_text(size = 12, hjust = 0, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          panel.border = element_rect(size = 0.75, colour = "gray20"),
          axis.ticks = element_line(size = 0.75, color = "gray20"),
          panel.grid = element_blank())
    
  ## Return:
  return(plot)
  
} ## End of function.

## Function: Mapping groups ----
map_group <- function(subclade) {
  
  #subclade <- "inornatus"
  
  ## Candidate species:
  if (subclade == "inornatus") { taxa <- c("inornatus-N", "inornatus-S") }
  if (subclade == "superciliaris") { taxa <- c("superciliaris-E (S)", "superciliaris-E (N)", "superciliaris-W") }
  if (subclade == "superciliaris_suppl") { taxa <- c("superciliaris-E (S)", "superciliaris-E (N)", "superciliaris-W") }
  if (subclade == "spaldingi_suppl") { taxa <- c("rimacola", "spaldingi-CY", "spaldingi-NE", "spaldingi-NW", "spaldingi-S", "spaldingi-TE") } 
  if (subclade == "spaldingi") { taxa <- c("spaldingi-CY", "spaldingi-NE", "spaldingi-S") } 
  if (subclade == "robustus") { taxa <- c("spaldingi-NW", "spaldingi-TE") } 
  
  ## Import distances estimated previously:
  dist_df <- read.csv(file = here(paste0("IBD/data/", subclade, ".csv")), header = TRUE)
  dist_df <- dist_df[!is.na(dist_df$Genetic_distance), ]
  
  ## Loading sample information:
  sample_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
  sample_info <- sample_info[sample_info$CANDIDATE_SP_III %in% taxa, ]
  if (subclade == "spaldingi_suppl") { 
    sample_info$CANDIDATE_SP_III[sample_info$CANDIDATE_SP_III %in% c("spaldingi-NW", "spaldingi-TE")] <- "spaldingi-NW-TE" }
  if (subclade == "superciliaris_suppl") { 
    sample_info$CANDIDATE_SP_III[sample_info$CANDIDATE_SP_III %in% c("superciliaris-E (S)", "superciliaris-E (N)")] <- "superciliaris-EN-ES" }
  if (subclade %in% c("robustus", "spaldingi_suppl")) { 
    sample_info$CANDIDATE_SP_III <- gsub(sample_info$CANDIDATE_SP_III, pattern = "spaldingi-NW", replacement = "robustus-NW")
    sample_info$CANDIDATE_SP_III <- gsub(sample_info$CANDIDATE_SP_III, pattern = "spaldingi-TE", replacement = "robustus-TE")
  }
  
  ## Get rid of apparent sample mix-ups:
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% setdiff(sample_info$SAMPLE_ID, c("WAMR_163003_Ct_duri")), ]
  
  ## Keep one sample per site for a cleaner map:
  sample_info <- sample_info %>% group_by(LAT, LON) %>% sample_n(1)
  
  ## Australia map:
  AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
  AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.
  
  ## Legend title:
  if (subclade == "robustus") { titlelabel <- expression(italic("C. robustus")) }
  if (subclade %in% c("spaldingi", "spaldingi_suppl")) { titlelabel <- expression(italic("C. spaldingi")) }
  if (subclade %in% c("superciliaris", "superciliaris_suppl")) { titlelabel <- expression(italic("C. superciliaris")) }
  if (subclade %in% c("inornatus", "inornatus_suppl")) { titlelabel <- expression(italic("C. inornatus")) }
  
  ## Color palette:
  if (subclade == "robustus") { 
    palette <- c("#3e2331", "white")
    symbols <- c(21, 23)
   }
  if (subclade == "spaldingi") { 
    palette <- c("gray50", "#e39c39", "#c9622f")
    symbols <- c(24, 21, 21)
  }
   if (subclade == "spaldingi_suppl") { 
    palette <- c("white", "gray50", "#e39c39", "#3e2331", "#c9622f")
    symbols <- c(22, 24, 21, 21, 21)
  }
  if (subclade == "superciliaris") { 
    palette <- c("#66bba3", "#66bba3", "#B5CF6B")
    symbols <- c(23, 22, 21)
  }
  if (subclade == "superciliaris_suppl") { 
    palette <- c("#66bba3","#162c3b", "#B5CF6B")
    symbols <- c(21, 21, 21)
  }
  if (subclade == "inornatus") { 
    palette <- c("#E7969C", "#7B4173")
    symbols <- c(21, 21)
  }
  if (subclade == "inornatus_suppl") { 
    palette <- c("gray50", "#E7969C", "#7B4173", "#843C39")
    symbols <- c(25, 21, 21, 21)
  }

  ## Plot map with ggplot:
  map <- ggplot() +
        
    ## Adding baseline map:
    geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray50", size = 0.3) +
        
    ## Setting a bunch of aesthetic parameters:
    theme_void() +
    scale_fill_manual(values = palette, name = titlelabel) +
    scale_shape_manual(values = symbols, name = titlelabel) +
       
    ## Setting some other elements:
    theme(panel.background = element_rect(fill = "aliceblue"),
          plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
          panel.border = element_rect(fill = NA, color = "gray20", size = 0.75),
          strip.background = element_blank(),
          legend.title = element_text(size = 10))
          
  ## Points on top:
  if (subclade == "robustus") {
  map <- map + geom_point(data = filter(sample_info, CANDIDATE_SP_III %in% setdiff(sample_info$CANDIDATE_SP_III, "robustus-TE")), aes(x = LON, y = LAT, fill = CANDIDATE_SP_III, shape = CANDIDATE_SP_III), 
               size = 2.5, color = "black", stroke = 0.4, alpha = 1)
  map <- map + geom_point(data = filter(sample_info, CANDIDATE_SP_III == "robustus-TE"), aes(x = LON, y = LAT, fill = CANDIDATE_SP_III, shape = CANDIDATE_SP_III), 
               size = 2.5, color = "black", stroke = 0.4, alpha = 1)
  } else {
  map <- map + geom_point(data = sample_info, aes(x = LON, y = LAT, fill = CANDIDATE_SP_III, shape = CANDIDATE_SP_III), 
               size = 2.5, color = "black", stroke = 0.4, alpha = 1)
  }
  
  ## Return:
  return(map)
    
} ## End of function.

## PART 3: Estimate IBD for each subclade ----

## Taxa and missing data:
subclades <- c("inornatus",
               "robustus",
               "spaldingi",
               "superciliaris",
               "inornatus_suppl",
               "spaldingi_suppl",
               "superciliaris_suppl")

## Estimate genetic distances:
for (subclade in subclades) { get_dist(subclade) }

## PART 4: Figure for manuscript ----

## Plot:
plot_robu1 <- plot_IBD("robustus")
map_robu1 <- map_group("robustus")
plot_spal1 <- plot_IBD("spaldingi")
map_spal1 <- map_group("spaldingi")
plot_supe1 <- plot_IBD("superciliaris")
map_supe1 <- map_group("superciliaris")
plot_inor1 <- plot_IBD("inornatus")
map_inor1 <- map_group("inornatus")

## Combine and save:
robu <- ( plot_robu1 + map_robu1 ) + plot_layout(widths = c(0.95, 1))
spal <- ( plot_spal1 + map_spal1 ) + plot_layout(widths = c(0.95, 1))
supe <- ( plot_supe1 + map_supe1 ) + plot_layout(widths = c(0.95, 1))
inor <- ( plot_inor1 + map_inor1 ) + plot_layout(widths = c(0.95, 1))
comb1 <- ( robu / spal / supe / inor )
ggsave(plot = comb1, height = 11, width = 8, units = "in", dpi = 600, 
       filename = paste0("IBD/plots/ms", miss_SNP, "_mi", miss_ind, "_snps", snps, ".pdf"))
ggsave(plot = comb1, height = 11, width = 8, units = "in", dpi = 600, 
       filename = paste0("IBD/plots/ms", miss_SNP, "_mi", miss_ind, "_snps", snps, ".jpg"))

## PART 5: Figure for supplementary material ----

## Plot:
plot_spal2 <- plot_IBD("spaldingi_suppl")
map_spal2 <- map_group("spaldingi_suppl")
plot_supe2 <- plot_IBD("superciliaris_suppl")
map_supe2 <- map_group("superciliaris_suppl")
plot_inor2 <- plot_IBD("inornatus_suppl") + theme(plot.margin = margin(t = 0, r = 0, b = 25, l = 25))
map_inor2 <- map_group("inornatus_suppl")

## For supplementary material:
comb2 <- ( plot_spal2 | map_spal2 ) / ( plot_supe2 | map_supe2 ) / ( plot_inor2 | map_inor2 ) +
   plot_annotation(title = 
   "Fig. S6. Pairwise Fst between individuals from same (gray) or different (black) genetic 
   groups, as a function of the geographic distances between them. Groups are the 
   operational taxonomic units (OTUs), whose ranges are shown on maps. Plots for the 
   spaldingi complex combine OTUs spaldingi-NW and NE. Plots for the superciliaris
   complex combine OTUs superciliaris-E (N) and superciliaris-E (S). Plots for the
   inornatus complex show all OTUs in this complex. For details, see text.",
                   theme = theme(plot.title = element_text(size = 14, hjust = 0, margin = margin(25, 25, 0, 25))))
ggsave(plot = comb2, height = 11, width = 8.5, units = "in", dpi = 300, 
       filename = "supplementary_files/Fig_S6_IBD.pdf")

## End of script.
