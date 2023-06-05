###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of the script are:
### To process outputs of Dsuite analyses.
### To make a heat map plot from the resulting D statistics values.
### Good resources:
### https://github.com/millanek/Dsuite
### https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data

## PART 1: Getting ready ----

## Packages:
library(reshape2)

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

## Samples:
#Dsuite_out <- read.table(header = TRUE, file = paste0("Dsuite/", dataset, "/", scheme, "_BBAA.txt")) ## Assuming no tree.
Dsuite_out <- read.table(header = TRUE, file = paste0("Dsuite/", dataset, "/", scheme, "_tree.txt")) ## Assuming tree provided.

## Adjusting p-values for multiple testing correction using Benjamini-Hochberg to control for the false discovery rate:
Dsuite_out$adj_pvalue <- p.adjust(Dsuite_out$p.value, method = "BH")

## Plot D-stats and p-values adjusted or not:
#plot(Dsuite_out$Dstatistic, ylab = "D" , xlab = "trio number")
#plot(Dsuite_out$p.value, ylab = "p value", xlab = "trio number", ylim = c(0, 1))
#plot(Dsuite_out$adj_pvalue, ylab = "p value", xlab = "trio number", ylim = c(0, 1))

## Remove spaces from OTU names:
Dsuite_out$P2 <- gsub(Dsuite_out$P2, pattern = "_", replacement = " ")
Dsuite_out$P3 <- gsub(Dsuite_out$P3, pattern = "_", replacement = " ")

## Extracting maximum D values and associated p:
Dsuite_out$pair <- paste0(Dsuite_out$P2, "_", Dsuite_out$P3)
get_max_D <- function(pair) {
  pair_df <- Dsuite_out[Dsuite_out$pair == pair, ]
  pair_stat <- pair_df[which.max(pair_df$Dstatistic), ]
  return(pair_stat)
}
Dsuite_max <- map_df(.x = unique(Dsuite_out$pair), .f = get_max_D)

## Plot order:
p_order <- read.table(file = paste0("Dsuite/", dataset, "/", scheme, "_plot_order.txt"), header = FALSE)
p_order <- setdiff(p_order$V1, "Outgroup")
p_order <- gsub(p_order, pattern = "_", replacement = " ")
Dsuite_max$P2 <- factor(Dsuite_max$P2, levels = rev(p_order))
Dsuite_max$P3 <- factor(Dsuite_max$P3, levels = p_order)

## Background tiles:
tile_df <- matrix(nrow = length(unique(Dsuite_max$P2)), ncol = length(unique(Dsuite_max$P2)))
rownames(tile_df) <- unique(Dsuite_max$P2)
colnames(tile_df) <- unique(Dsuite_max$P2)
tile_df[is.na(tile_df)] <- 1
tile_df <- melt(tile_df)
names(tile_df) <- c("P2", "P3", "value")
tile_df$P2 <- gsub(tile_df$P2, pattern = "_", replacement = " ")
tile_df$P3 <- gsub(tile_df$P3, pattern = "_", replacement = " ")
tile_df$P2 <- factor(tile_df$P2, levels = rev(p_order))
tile_df$P3 <- factor(tile_df$P3, levels = p_order)

## Heatmap:
plot_h <- ggplot() +
  geom_tile(data = tile_df, aes(x = P3, y = P2), fill = "white", colour = "gray80", size = 0.75) +
  geom_tile(data = Dsuite_max, aes(x = P3, y = P2, fill = Dstatistic), colour = "gray80", size = 0.75) +
  geom_tile(data = filter(Dsuite_max, Dsuite_max$adj_pvalue < 0.05), aes(x = P3, y = P2), 
            fill = "transparent", colour = "black", size = 1) +
  ggtitle(label = "Figure S6. Excess allele sharing among operational taxonomic units\nbased on the D-statistic. Black outlines indicate significant values.") +
  scale_fill_gradientn(colors = c("white", "#fdcdb8ff", "#fc8d6dff", "#f24734ff", "#bd151aff", "#6c010eff")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, margin = margin(t = 0, r = 0, b = 20, l = 0)),
        plot.margin = margin(t = 160, r = 20, b = 160, l = 20))
plot_h

## Save:
#ggsave(plot = plot_h, filename = "supplementary_files/Fig_S7_Dstat_heatmap.jpg", width = 8.5, height = 1, units = "in", dpi = 300, bg = "white")
ggsave(plot = plot_h, filename = "supplementary_files/Fig_S7_Dstat_heatmap.pdf", width = 8.5, height = 11, units = "in", dpi = 300, bg = "white")

## End of script.

## AFTERTHOUGHT: Gene flow versus time ----

## Pairs:
pairs <- c("inor-H_supe-W", ## Running.
           "inor-N_supe-W",
           "inor-N_supe-K",
           "inor-N_supe-E",
           "inor-N_spal-Rx",
           "inor-H_inor-N",
           "supe-K_supe-W",
           "supe-E_supe-K",
           "supe-E_supe-W",
           "spal-C_spal-N",
           "spal-C_spal-S",
           "spal-N_spal-S",
           "spal-N_spal-Rx",
           # "inor-B_inor-D", ## Simpson.
           # "inor-B_inor-M", ## Gibson.
           # "inor-F_inor-P",
           # "inor-F_inor-M",
           # "inor-P_inor-S",
           # "inor-M_inor-S",
           # "spal-B_spal-Rrr",
           #"inor-P1_inor-P2",
           #"inor-M1_inor-M2", 
           #"inor-N1_inor-N2",
           "inor-M_inor-P",
           "inor-N_inor-P",
           "inor-M_inor-N")

## PART 4: Function: Parameter summaries ----
get_stats <- function(pair) {

  ## Testing:
  #pair <- "inor-N_spal-Rx"
  #pair <- "inor-N_inor-P"
  
  ## Status:
  print(paste0("Now doing pair ", pair, "!"))
  
  ## data:
  out <- read.csv(file = paste0("G-PhoCS/data/", pair, ".csv"), header = TRUE)
  stats <- out[c("N_a", "N_b", "T_ab", "M_ab", "M_ba", "NM2_ab", "NM2_ba")]
  
  ## Median values:
  stats <- data.frame(Median = apply(stats, 2, median, na.rm = TRUE))
  stats <- as.data.frame(t(stats))
  stats$NM2_avg <- (stats$NM2_ab + stats$NM2_ba)/2
  stats$M_avg <- (stats$M_ab + stats$M_ba)/2
  stats$pair <- pair
  
  ## Return:
  return(stats)
  
} ## End of function.

## Apply:
stat_df <- map_df(pairs, get_stats)
stat_df$pair <- gsub(stat_df$pair, pattern = "spal-C", replacement = "spal-CY")
stat_df$pair <- gsub(stat_df$pair, pattern = "spal-Rx", replacement = "spal-NW")
stat_df$pair <- gsub(stat_df$pair, pattern = "inor-H", replacement = "inor-S")
stat_df$pair <- gsub(stat_df$pair, pattern = "supe", replacement = "superciliaris")
stat_df$pair <- gsub(stat_df$pair, pattern = "inor", replacement = "inornatus")
stat_df$pair <- gsub(stat_df$pair, pattern = "spal", replacement = "spaldingi")

## Plots:
#plot(y = stat_df$T_ab, x = stat_df$NM2_ab)
#plot(y = stat_df$T_ab, x = stat_df$NM2_ba)
#plot(y = stat_df$T_ab, x = stat_df$M_ba)
#plot(y = stat_df$T_ab, x = stat_df$M_ab)
plot(y = stat_df$T_ab, x = stat_df$NM2_avg)
plot(y = stat_df$T_ab, x = stat_df$M_avg)
ggplot(data = stat_df, mapping = aes(x = NM2_avg, y = log10(T_ab))) + theme_classic() +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

## ABBA-BABA stats:
Dmax_df <- Dsuite_max 
Dmax_df$pair <- gsub(Dmax_df$pair, pattern = " ", replacement = "-")
Dmax_df$pairinv <- gsub(Dmax_df$pair, pattern = "(.+)_(.+)", replacement = "\\2_\\1")

## Merge:
merge_1 <- merge(stat_df, Dmax_df, by = "pair")
merge_2 <- merge(stat_df, Dmax_df, by.x = "pair", by.y = "pairinv")
merge_1 <- merge_1[setdiff(names(merge_1), "pairinv")]
merge_2 <- merge_2[setdiff(names(merge_2), "pair")]
names(merge_2)[names(merge_2) == "pair.y"] <- "pair"
merged <- rbind(merge_1, merge_2)

## Plots:
plot(y = merged$T_ab, x = merged$Dstatistic)
plot(y = merged$T_ab, x = merged$f4.ratio)

## End of Afterthought. 
