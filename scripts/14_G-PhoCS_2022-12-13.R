###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To process outputs of G-PhoCS analyses.
### To make plots.

## PART 1: Getting ready ----

## Packages:
library(HDInterval)
library(patchwork)
library(raster)
library(reshape2)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Paths:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/" ## Suindara.
setwd(path)
#detach("package:here", unload = TRUE)
library(here)

## Groups considered: OTUs, taxa, and subpopulations within each complex:
inornatus <- c("Murchison", "Murchison1", "Murchison2", "Pilbara", "Pilbara1", "Pilbara2", "Simpson", 
               "brachyonyx", "fallens", "severus", "inornatus-N", "inornatus-N1", "inornatus-N2", 
               "inornatus-S", "inornatus-E", "lateralis")
spaldingi <- c("borealis", "spaldingi-CY", "spaldingi-NE", "spaldingi-NW", "spaldingi-S", "spaldingi-TE")
superciliaris <- c("superciliaris-E", "superciliaris-EN", "superciliaris-ES", "superciliaris-K", "superciliaris-W")

## Burn-in:
burnin <- 75000
steps <- 300000

## Mutation rate:
u = 7.6*10^-9 ## From Gottscho et al. (2017 MPE).

## Generation time:
g <- 2  ## 22 months in Ctenotus "helenae" from James 1991 Copeia and James 1991 J. Herpetology.

## PART 2: Function: Estimate pop. gen. stats and change population labels  ----
get_stats <- function(pair) {

  ## Testing:
  #pair <- "inornatus-N_superciliaris-K"
  #pair <- "Murchison_severus"
  
  ## Status:
  print(paste0("Now doing pair ", pair, "!"))
  
  ## Outfiles for pair:
  files <- list.files(path = "G-PhoCS/great_lakes/", pattern = paste0(pair, "_run"))
  
  ## For pairs with more than three runs:
  pair_12 <- c("Pilbara1_Pilbara2") #, "brachyonyx_Simpson")
  if (pair %in% pair_12) { files <- files[c(1, 2)] }
  pair_13 <- c("spaldingi-NW_superciliaris-W", "spaldingi-CY_spaldingi-NE") #, "brachyonyx_Simpson")
  if (pair %in% pair_13) { files <- files[c(1, 3)] }
  pair_23 <- c("borealis_spaldingi-NW", "brachyonyx_Murchison", "Murchison_severus", #"brachyonyx_Simpson",
               "Pilbara_severus", "spaldingi-NW_superciliaris-E", "Murchison1_Murchison2")
  if (pair %in% pair_23) { files <- files[c(2, 3)] }
  pair_14 <- c("inornatus-E_lateralis", "brachyonyx_Simpson")
  if (pair %in% pair_14) { files <- files[c(1, 4)] }
  
  ## Data:
  get_out <- function(file) {
    out <- read.table(file = paste0("G-PhoCS/great_lakes/", file), header = TRUE)
    out <- out[1:steps, ] ## use the same maximum number of steps across all runs (some went longer).
    out <- out[burnin:nrow(out), ] ## Remove burn-in (as number of steps).
    out <- out[seq(1, nrow(out), 100), ] ## Keep every 100th step.
    return(out)
  }
  out <- map_df(files, get_out)
  out <- out[rowSums(is.na(out)) != ncol(out), ] ## Remove NAs.
  
  ## Estimate M and average:
  out$M_ab <- out$m_pop_a..pop_b*u*g
  out$M_ba <- out$m_pop_b..pop_a*u*g
  out$M_avg <- (out$M_ab + out$M_ba)/2
  
  ## Estimate N and T:
  out$N_a <- out$theta_pop_a/(4*u*g) 
  out$N_b <- out$theta_pop_b/(4*u*g) 
  out$T_ab <- out$tau_root/u 
  
  ## Estimate 2NM: 
  #out$NM2_ab <- 2*out$N_b*out$M_ab
  #out$NM2_ba <- 2*out$N_a*out$M_ba
  ## Or:
  out$NM2_ab <- out$m_pop_a..pop_b*out$theta_pop_b/2
  out$NM2_ba <- out$m_pop_b..pop_a*out$theta_pop_a/2
  
  ## Highest 2NM band:
  if (mean(out$NM2_ab, na.rm = TRUE) >  mean(out$NM2_ba, na.rm = TRUE)) { out$NM2_high <- out$NM2_ab 
  } else { out$NM2_high <- out$NM2_ba }
  
  ## Average:
  out$NM2_avg <- (out$NM2_ab + out$NM2_ba)/2
  
  ## Add complex info:
  complex1 <- gsub(x = pair, pattern = "(.+)_.+", replacement = "\\1")
  complex2 <- gsub(x = pair, pattern = ".+_(.+)", replacement = "\\1")
  complex1 <- gsub(x = complex1, pattern = "(.+)-.+", replacement = "\\1")
  complex2 <- gsub(x = complex2, pattern = "(.+)-.+", replacement = "\\1")
  if (complex1 %in% inornatus | complex1 == "inornatus") { complex1 <- "inornatus" }
  if (complex2 %in% inornatus | complex2 == "inornatus") { complex2 <- "inornatus" }
  if (complex1 %in% spaldingi | complex1 == "spaldingi") { complex1 <- "spaldingi" }
  if (complex2 %in% spaldingi | complex2 == "spaldingi") { complex2 <- "spaldingi" }
  if (complex1 %in% superciliaris | complex1 == "superciliaris") { complex1 <- "superciliaris" }
  if (complex2 %in% superciliaris | complex2 == "superciliaris") { complex2 <- "superciliaris" }
  if (complex1 == complex2) { out$complex <- complex1 }
  if (complex1 != complex2) { out$complex <- "cross-complex" }
  
  ## Add pair info:
  out$pair <- pair
  
  ## Changing labels:
  out$pair <- gsub(out$pair, pattern = "brachyonyx", replacement = "C. \"brachyonyx\"")
  out$pair <- gsub(out$pair, pattern = "fallens", replacement = "C. \"fallens\"")
  out$pair <- gsub(out$pair, pattern = "severus", replacement = "C. \"severus\"")
  out$pair <- gsub(out$pair, pattern = "borealis", replacement = "C. \"borealis\"")
  
  ## Split pair:
  out <- separate(data = out, col = pair, into = c("pop_a", "pop_b"), sep = "_")
  
  ## Save complete dataframe:
  write.csv(out, file = paste0("G-PhoCS/data/", pair, ".csv"), row.names = FALSE)
  
  ## Return average 2NM estimates for plots:
  return(out[c("complex", "pop_a", "pop_b", "NM2_avg", "M_avg", "NM2_high")])
  
} ## End of function.

## PART 3: Function: Assessing mixing ----
assess_mix <- function(pair) {

  ## Testing:
  #pair <- "Murchison_inornatus-N"
  #pair <- "lateralis_superciliaris-E"
  #pair <- "Murchison1_Murchison2"
  
  ## Status:
  print(paste0("Now doing pair ", pair, "!"))
  
  ## Outfiles for pair:
  files <- list.files(path = "G-PhoCS/great_lakes/", pattern = paste0(pair, "_run"))
  
  ## Data:
  out_ls <- list()
  for (i in 1:length(files)) {
    out <- read.table(file = paste0("G-PhoCS/great_lakes/", files[i]), header = TRUE)
    out <- out[seq(1, nrow(out), 100), ] ## Keep every 100th step.
    out$NM2_ab <- out$m_pop_a..pop_b*out$theta_pop_b/2
    out$NM2_ba <- out$m_pop_b..pop_a*out$theta_pop_a/2
    out$NM2_avg <- (out$NM2_ab + out$NM2_ba)/2
    out_ls[[i]] <- out 
  }
  
  ## Plot limits (y-axis):
  get_lim <- function(i) {
    return(data.frame(max_NM2_ab = max(out_ls[[i]]$NM2_ab),
                      max_NM2_ba = max(out_ls[[i]]$NM2_ba),
                      max_NM2_avg = max(out_ls[[i]]$NM2_avg),
                      min_NM2_ab = min(out_ls[[i]]$NM2_ab),
                      min_NM2_ba = min(out_ls[[i]]$NM2_ba),
                      min_NM2_avg = min(out_ls[[i]]$NM2_avg)))
  } ## End of function.
  lim_df <- map_df(1:length(files), get_lim)  
  
  ## Plot trace:
  jpeg(filename = paste0("G-PhoCS/traces/", pair, ".jpg"), res = 200, width = 15, height = 8, units = "in")
    par(mfrow = c(length(files), 3))
    for (i in 1:length(files)) {
      plot(out_ls[[i]]$NM2_ab, main = "a --> b", ylab = "2NM", cex = 0.75, ylim = c(min(lim_df$min_NM2_ab), max(lim_df$max_NM2_ab)))
      plot(out_ls[[i]]$NM2_ba, main = "b --> a", ylab = "2NM", cex = 0.75, ylim = c(min(lim_df$min_NM2_ba), max(lim_df$max_NM2_ba)))
      plot(out_ls[[i]]$NM2_avg, main = "b <--> a", ylab = "2NM", cex = 0.75, ylim = c(min(lim_df$min_NM2_avg), max(lim_df$max_NM2_avg)))
    }
  dev.off()
  
} ## End of function.
    
## PART 4: Function: Parameter summaries ----
summ_stats <- function(pair) {

  ## Testing:
  #pair <- "inornatus-N_Pilbara"
  
  ## Status:
  print(paste0("Now doing pair ", pair, "!"))
  
  ## data:
  out <- read.csv(file = paste0("G-PhoCS/data/", pair, ".csv"), header = TRUE)
  out_r <- out[c("N_a", "N_b", "T_ab", "M_ab", "M_ba", "NM2_ab", "NM2_ba")]
  
  ## Median values:
  medians <- apply(out_r, 2, median, na.rm = TRUE)
  
  ## Highest posterior densities:
  hpds <- as.data.frame(hdi(out_r))
   
  ## Combine and name rows:
  stats <- rbind(medians, hpds)
  row.names(stats) <- c("Median", "Lower 95% HPD", "Upper 95% HPD")
  
  ## Round:
  stats$N_a <- round(stats$N_a, 0)
  stats$N_b <- round(stats$N_b, 0)
  stats$T_ab <- round(stats$T_ab, 0)
  stats$M_ab <- sprintf(stats$M_ab, fmt = "%#.10f")
  stats$M_ba <- sprintf(stats$M_ba, fmt = "%#.10f")
  stats$NM2_ab <- round(stats$NM2_ab, 2)
  stats$NM2_ba <- round(stats$NM2_ba, 2)
  
  ## Transpose:
  stats <- as.data.frame(t(stats))
  
  ## Name rows:
  pop_a <- unique(out$pop_a)
  pop_b <- unique(out$pop_b)
  stats$Parameter <- c(paste0("N ", pop_a),
                       paste0("N ", pop_b),
                       paste0("T ", pop_a, " + ", pop_b),
                       paste0("M ", pop_a, " --> ", pop_b),
                       paste0("M ", pop_b, " --> ", pop_a),
                       paste0("2NM ", pop_a, " --> ", pop_b),
                       paste0("2NM ", pop_b, " --> ", pop_a))
  return(stats[c("Parameter", "Median", "Lower 95% HPD", "Upper 95% HPD")])
             
} ## End of function.

## PART 5: Function: Plots ----
get_plot <- function(dataset) {
  
  ## Testing:
  #dataset <- "OTUs"
  #dataset <- "subpops"
  
  ## Title:
  if (dataset == "cross") { 
    labtitle <- "Cross-complex pairs" 
    palette <- "#152868" }
  if (dataset == "inor") { 
    labtitle <- "inornatus complex" 
    palette <- "#15b6ad" }
  if (dataset == "supe") { 
    labtitle <- "superciliaris complex" 
    palette <- "#f3e2c6" }
  if (dataset == "spal") { 
    labtitle <- "spaldingi complex" 
    palette <- "orange" }
  
  ## Data:
  df <- df_list[[dataset]]
  df$pair <- paste0(df$pop_a, " vs. ", df$pop_b)
  
  ## Remove randomly split pairs from plot:
  remove <- c("inornatus-N1 vs. inornatus-N2", "Murchison1 vs. Murchison2", "Pilbara1 vs. Pilbara2", "Simpson1 vs. Simpson2")
  df <- df[df$pair %in% setdiff(df$pair, remove), ]
  
  ## Switch order and renaming units:
  df$pair <- gsub(df$pair, pattern = "Murchison vs. C. \"severus\"", replacement = "C. \"severus\" vs. Murchison")
  df$pair <- gsub(df$pair, pattern = "Murchison vs. inornatus-N", replacement = "inornatus-N vs. Murchison")
  df$pair <- gsub(df$pair, pattern = "Murchison vs. Pilbara", replacement = "Pilbara vs. Murchison")
  df$pair <- gsub(df$pair, pattern = "Pilbara vs. C. \"severus\"", replacement = "C. \"severus\" vs. Pilbara")
  df$pair <- gsub(df$pair, pattern = "inornatus-N vs. Pilbara", replacement = "inornatus-N vs. Pilbara")
  df$pair <- gsub(df$pair, pattern = "Murchison", replacement = "western shrublands")
  df$pair <- gsub(df$pair, pattern = "Simpson", replacement = "central deserts")
  df$pair <- gsub(df$pair, pattern = "inornatus-E", replacement = "eutaenius")
  df$pair <- gsub(df$pair, pattern = "superciliaris-K", replacement = "mastigura")
  df$pair <- gsub(df$pair, pattern = "spaldingi-NW", replacement = "robustus-NW")
  df$pair <- gsub(df$pair, pattern = "spaldingi-TE", replacement = "robustus-TE")
  df$pair <- gsub(df$pair, pattern = "C. \"borealis\"", replacement = "\"borealis\"")
  df$pair <- gsub(df$pair, pattern = "C. \"brachyonyx\"", replacement = "\"brachyonyx\"")
  df$pair <- gsub(df$pair, pattern = "C. \"fallens\"", replacement = "\"fallens\"")
  df$pair <- gsub(df$pair, pattern = "C. \"severus\"", replacement = "\"severus\"")
  
  ## Reorder:
  df$complex <- factor(df$complex, levels = c("spaldingi", "superciliaris", "inornatus", "cross-complex"))
  df$pair <- factor(df$pair, levels = unique(df$pair))
  df <- arrange(df, pair)
  
  ## Plot:
  plot_d <- ggplot() + 
    
    ## Line, plot:
    geom_violin(data = filter(df, !is.na(NM2_avg)), aes(x = NM2_avg, y = pair, fill = complex), 
                trim = FALSE, scale = "width", alpha = 0.90, color = "black", size = 0.3) +
    
    ## Titles:
    labs(title = labtitle, fill = "Complex", x = "2NM") +
    
    ## Axes:
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0, 1),
                       breaks = waiver(), n.breaks = 10) +
    
    ## Aesthetics:
    scale_fill_manual(values = palette) +
    theme_classic() +
    theme(#legend.position = c(0.85, 0.170),
          legend.position = "none",
          legend.background = element_rect(fill = "gray90"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "gray90"),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks = element_line(size = 0.5, color = "black"),
          plot.title = element_text(size = 12, margin = margin(t = 10, l = 0, b = 10, r = 0)),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 14, margin = margin(t = 10, l = 0, b = 0, r = 0)),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 10))
  plot_d
  
  ## Return:
  return(plot_d)

} ## End of function.

## PART 6: Define target pairs ----

## Pairs: 
p_cross <- c("inornatus-S_superciliaris-K", 
             "inornatus-S_superciliaris-W", 
             "inornatus-S_superciliaris-ES", 
             "inornatus-S_superciliaris-EN",
             "inornatus-N_superciliaris-K", 
             "inornatus-N_superciliaris-W", 
             "inornatus-N_superciliaris-ES", 
             "inornatus-N_superciliaris-EN",            
             "inornatus-N_spaldingi-NW", 
             "inornatus-N_spaldingi-TE", 
             "lateralis_spaldingi-S", 
             "spaldingi-NW_superciliaris-W") 
            
p_inor <- c("Pilbara1_Pilbara2",
            "Murchison1_Murchison2", 
            "inornatus-N1_inornatus-N2",
            "Murchison_Pilbara",
            "Murchison_inornatus-N",
            "inornatus-N_Pilbara",
            "brachyonyx_Simpson",
            "brachyonyx_Murchison", 
            "fallens_Pilbara",
            "fallens_Murchison",
            "Pilbara_severus",
            "Murchison_severus",
            "inornatus-N_inornatus-S", 
            "inornatus-S_lateralis", 
            "inornatus-E_lateralis")
      
p_supe <- c("superciliaris-K_superciliaris-W", 
            "superciliaris-ES_superciliaris-K", 
            "superciliaris-ES_superciliaris-W", 
            "superciliaris-EN_superciliaris-K", 
            "superciliaris-EN_superciliaris-W", 
            "superciliaris-EN_superciliaris-ES")

p_spal <- c("borealis_spaldingi-NW",
            "spaldingi-NE_spaldingi-S", 
            "spaldingi-CY_spaldingi-NE", 
            "spaldingi-NE_spaldingi-NW", 
            "spaldingi-CY_spaldingi-S", 
            "spaldingi-NW_spaldingi-TE")

## Failed:
            #"inornatus-E_inornatus-S",
            #"lateralis_superciliaris-EN", 
            #"inornatus-S_spaldingi-NW", 
            #"inornatus-S_spaldingi-TE",
            #"inornatus-N_spaldingi-NE",
            #"inornatus-S_superciliaris-E",
            #"inornatus-N_superciliaris-E",
            #"lateralis_superciliaris-W", 
            #"lateralis_superciliaris-E",
            #"spaldingi-NW_superciliaris-E",
            #"superciliaris-E_superciliaris-K",
            #"superciliaris-E_superciliaris-W",
            #"spaldingi-TE_superciliaris-W", ##
            #"spaldingi-TE_superciliaris-EN", ## Failed.            
            #"inornatus-E_spaldingi-NE", ## Failed.
            #"inornatus-E_spaldingi-S", ## Failed.
            #"lateralis_spaldingi-NE", ## Failed.
            #"lateralis_superciliaris-ES", ## Failed.
            #"spaldingi-NW_superciliaris-K", ## Failed.
            #"spaldingi-TE_superciliaris-ES", ## Failed.
            #"spaldingi-TE_superciliaris-K") ## Failed.

## PART 7: Assess mixing and convergence ----

## Traces:
map_df(.x = p_cross, .f = assess_mix)
map_df(.x = p_inor, .f = assess_mix)
map_df(.x = p_supe, .f = assess_mix)
map_df(.x = p_spal, .f = assess_mix)

## PART 8: Estimate migration ----

## Data:
df_list <- vector("list")
df_list$cross <- map_df(.x = p_cross, .f = get_stats)
df_list$inor <- map_df(.x = p_inor, .f = get_stats)
df_list$supe <- map_df(.x = p_supe, .f = get_stats)
df_list$spal <- map_df(.x = p_spal, .f = get_stats)

## PART 9: Summary stats ----
s_cross <- map_df(.x = p_cross, .f = summ_stats)
s_inor <- map_df(.x = p_inor, .f = summ_stats)
s_supe <- map_df(.x = p_supe, .f = summ_stats)
s_spal <- map_df(.x = p_spal, .f = summ_stats)
params <- rbind(s_cross, s_inor, s_supe, s_spal)
write.csv(params, file = "G-PhoCS/data/G-PhoCS_parameters.csv", row.names = FALSE)

## PART 10: Plot average 2NM across pairs ----

## Plots:
plot_cross <- get_plot(dataset = "cross")
plot_inor <- get_plot(dataset = "inor") 
plot_supe <- get_plot(dataset = "supe") 
plot_spal <- get_plot(dataset = "spal") 

## Combine and add title:
plot_a <- ( plot_spal + theme(axis.title.x = element_blank()) ) / 
          ( plot_supe + theme(axis.title.x = element_blank())  ) /
          ( plot_inor + theme(axis.title.x = element_blank())  ) /
            plot_cross 
plot_a <- plot_a + plot_layout(heights = c(1, 1, 2, 2))
plot_a <- plot_a +  
                    plot_annotation(title = "Population migration rate (2NM)",
                    subtitle = "Effective number of genes received by a population per generation (both directions averaged)",
                    theme = theme(plot.title = element_text(size = 16),
                                  plot.subtitle = element_text(size = 12),
                                  plot.margin = margin(t = 20, r = 20, b = 10, l = 20)))
## Check:
plot_a

## Save:
ggsave(plot = plot_a, filename = paste0("G-PhoCS/plot_G-PhoCS_2NM_average_bi_bycomplex", format(burnin, scientific = FALSE), ".jpg"), dpi = 300, units = "in", width = 8.5, height = 10)
ggsave(plot = plot_a, filename = paste0("G-PhoCS/plot_G-PhoCS_2NM_average_bi_bycomplex", format(burnin, scientific = FALSE), ".pdf"), dpi = 300, units = "in", width = 8.5, height = 10)

## PART 11: Maps of units ----

## Units:
units_ls <- vector("list")

## OTUs:
units_ls$"inornatus-E" <- c("SAMR_55874_Ct_euta", "NA_CCM0090_Ct_late", "SAMAR_55799_Ct_late")
units_ls$"inornatus-N" <- c("NA_CCM1710_Ct_inor", "NA_CCM0848_Ct_inor", "NA_CCM1249_Ct_inor", "NA_CCM1444_Ct_inor", "NA_CCM1445_Ct_inor", "NA_CCM1530_Ct_robu", "NA_CCM1580_Ct_inor", "NA_CCM3739_Ct_inor", "NA_CCM4920_Ct_essi", "NTMR_20687_Ct_robu", "SAMR_34123_Ct_hill", "WAMR_126017_Ct_inor")
units_ls$"inornatus-S" <- c("CUMV_14611_Ct_hele", "NA_ABTC31797_Ct_saxa", "SAMR_36252_Ct_brac", "SAMR_42918_Ct_saxa", "SAMR_50208_Ct_saxa", "SAMR_54008_Ct_saxa", "UMMZ_242614_Ct_hele", "UMMZ_244288_ct_inor", "WAMR_094929_Ct_hele", "WAMR_129923_Ct_hele", "WAMR_135396_Ct_hele", "WAMR_139523_Ct_hele")
units_ls$"lateralis" <- c("SAMAR_34261_Ct_late", "SAMAR_42687_Ct_late", "SAMAR_42768_Ct_late", "SAMAR_42837_Ct_late", "SAMAR_54433_Ct_late", "SAMAR_54463_Ct_late", "SAMR_42758_Ct_late", "SAMR_54340_Ct_late", "SAMR_54434_Ct_late", "SAMR_65339_Ct_late", "SAMR_65407_Ct_late", "SAMR_65416_Ct_late")
units_ls$"spaldingi-CY" <- c("AMSR_111494_Ct_spal", "AMSR_111493_Ct_spal")
units_ls$"spaldingi-NE" <- c("NTMR_21670_Ct_spal", "SAMAR_54485_Ct_spal", "SAMAR_51093_Ct_spal", "SAMR_55725_Ct_spal", "SAMAR_55675_Ct_spal", "NTMR_22299_Ct_spal", "NA_CCM2497_Ct_spal", "NTMR_22298_Ct_spal", "QM_82086_Ct_spal", "NTMR_22622_Ct_spal", "NTMR_22325_Ct_spal", "SAMAR_34202_Ct_spal")
units_ls$"spaldingi-NW" <- c("NA_CCM1389_Ct_robu", "NA_CCM1412_Ct_robu", "NA_CCM1415_Ct_robu", "NA_CCM2866_Ct_robu", "NTMR_13838_Ct_bore", "NTMR_17739_Ct_bore", "NTMR_22175_Ct_robu", "NTMR_22620_Ct_robu", "NTMR_22934_Ct_robu", "WAMR_137950_Ct_robu", "WAMR_141379_Ct_robu", "WAMR_145686_Ct_robu")
units_ls$"spaldingi-TE" <- c("NTMR_17738_Ct_cogg", "NTMR_20378_Ct_robu", "NTMR_22166_Ct_robu")
units_ls$"spaldingi-S" <- c("SAMAR_55731_Ct_robu", "SAMR_28536_Ct_robu", "SAMR_33525_Ct_robu", "SAMR_33560_Ct_robu", "SAMR_33695_Ct_robu", "SAMR_33876_Ct_robu", "SAMR_36313_Ct_robu", "SAMR_39480_Ct_robu", "SAMR_40718_Ct_robu", "SAMR_47028_Ct_robu" , "SAMR_48962_Ct_robu", "SAMR_55746_Ct_robu")
#units_ls$"superciliaris-E" <- c("SAMR_53975_Ct_saxa", "NA_CCM2565_Ct_vert", "NA_CCM2518_Ct_vert", "NA_CCM2292_Ct_inor", "NA_CCM3869_Ct_inor", "SAMR_38776_Ct_saxa", "SAMAR_34180_Ct_bore", "NTMR_26117_Ct_inor")
units_ls$"superciliaris-EN" <- c("NA_CCM2565_Ct_vert", "NA_CCM3869_Ct_inor", "NTMR_26117_Ct_inor", "NTMR_20846_Ct_saxa", "NA_CCM2518_Ct_vert", "NA_CCM2292_Ct_inor")
units_ls$"superciliaris-ES" <- c("SAMR_38776_Ct_saxa", "SAMR_53975_Ct_saxa", "NTMR_16340_Ct_inor", "SAMAR_34180_Ct_bore", "NTMR_22432_Ct_saxa")
units_ls$"superciliaris-K" <- c("NA_CCM0980_Ct_mast", "NA_CCM1579_Ct_inor", "NA_CCM0830_Ct_inor", "NA_CCM0823_Ct_robu", "NA_CCM1061_Ct_sp", "NA_CCM1390_Ct_inor", "NA_CCM1529_Ct_inor", "NA_CCM0946_Ct_inor", "NA_CCM1206_Ct_robu", "NA_CCM0847_Ct_robu", "NA_CCM0753_Ct_robu", "NA_CCM0850_Ct_inor")
units_ls$"superciliaris-W" <- c("NA_CCM0535_Ct_inor", "NTMR_25983_Ct_spal", "UMMZ_242625_Ct_inor", "UMMZ_242630_Ct_inor", "WAMR_108698_Ct_saxa", "WAMR_108766_Ct_saxa", "WAMR_108816_Ct_saxa", "WAMR_132522_Ct_saxa", "WAMR_132682_Ct_saxa", "WAMR_153812_Ct_saxa", "WAMR_158204_Ct_saxa", "WAMR_158376_Ct_saxa")

## Subgroups, including taxa:
units_ls$Murchison <- c("CUMV_14600_Ct_hele", "CUMV_14601_Ct_hele", "CUMV_14602_Ct_hele", "CUMV_14604_Ct_hele", "CUMV_14606_Ct_hele", "CUMV_14611_Ct_hele", "CUMV_14612_Ct_hele", "CUMV_14655_Ct_hele", "CUMV_14656_Ct_hele", "WAMR_102671_Ct_hele", "WAMR_141131_Ct_hele", "WAMR_145926_Ct_hele")
units_ls$Murchison1 <- c("CUMV_14602_Ct_hele", "CUMV_14604_Ct_hele", "CUMV_14611_Ct_hele", "CUMV_14612_Ct_hele", "CUMV_14655_Ct_hele", "WAMR_102671_Ct_hele")
units_ls$Murchison2 <- c("CUMV_14600_Ct_hele", "CUMV_14601_Ct_hele", "CUMV_14606_Ct_hele", "CUMV_14656_Ct_hele", "WAMR_141131_Ct_hele", "WAMR_145926_Ct_hele")
units_ls$"inornatus-N1" <- c("NA_CCM4920_Ct_essi", "NA_CCM1580_Ct_inor", "NA_CCM1445_Ct_inor", "NA_CCM1444_Ct_inor", "NA_CCM1710_Ct_inor", "SAMR_34123_Ct_hill")
units_ls$"inornatus-N2" <- c("NA_CCM0848_Ct_inor", "NA_CCM1249_Ct_inor", "NA_CCM1530_Ct_robu", "NA_CCM3739_Ct_inor", "NTMR_20687_Ct_robu", "WAMR_126017_Ct_inor")
units_ls$Pilbara <- c("WAMR_102423_Ct_saxa", "WAMR_129923_Ct_hele", "WAMR_131009_Ct_hele", "WAMR_131371_Ct_hele", "WAMR_135396_Ct_hele", "WAMR_139296_Ct_hele", "WAMR_139523_Ct_hele", "WAMR_140720_Ct_hele", "WAMR_141301_Ct_hele", "WAMR_145567_Ct_hele", "WAMR_145698_Ct_hele", "WAMR_157646_Ct_hele")
units_ls$Pilbara1 <- c("WAMR_102423_Ct_saxa", "WAMR_131009_Ct_hele", "WAMR_139523_Ct_hele", "WAMR_140720_Ct_hele", "WAMR_141301_Ct_hele", "WAMR_145698_Ct_hele")
units_ls$Pilbara2 <- c("WAMR_129923_Ct_hele", "WAMR_131371_Ct_hele", "WAMR_135396_Ct_hele", "WAMR_139296_Ct_hele", "WAMR_145567_Ct_hele", "WAMR_157646_Ct_hele")
units_ls$Simpson <- c("NMVD_67682_Ct_saxa", "NMVD_67793_Ct_leae", "NA_ABTC31797_Ct_saxa", "NTMR_20628_Ct_saxa", "SAMR_46880_Ct_saxa", "SAMR_50208_Ct_saxa", "SAMR_50209_Ct_saxa")
units_ls$"C. \"brachyonyx\"" <- c("SAMR_29707_Ct_brac", "SAMR_36252_Ct_brac", "SAMR_36349_Ct_brac", "SAMR_41330_Ct_brac", "SAMR_45215_Ct_brac")
units_ls$"C. \"fallens\"" <- c("WAMR_131018_Ct_fall", "WAMR_154016_Ct_fall", "UMMZ_244288_ct_inor")
units_ls$"C. \"severus\"" <- c("WAMR_146913_Ct_seve", "WAMR_152991_Ct_seve", "WAMR_156159_Ct_seve")
units_ls$"C. \"borealis\"" <- c("NTMR_13838_Ct_bore", "NTMR_17739_Ct_bore", "NTMR_22167_Ct_bore", "NTMR_23004_Ct_bore")

## Samples:
sample_info <- read.csv("sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)

## For units:
get_coords <- function(unit) {
  samples <- units_ls[[unit]]
  unit_df <- sample_info[sample_info$SAMPLE_ID %in% samples, ]
  unit_df$unit <- unit
  unit_df <- unit_df[c("unit", "LAT", "LON")]
  return(unit_df)
} ## End of function.
map_df <- map_df(names(units_ls), get_coords)
map_df$unit <- factor(map_df$unit, levels = names(units_ls))

## AUS map:
AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.

## Plot:
map <- ggplot(data = map_df) + theme_void() +
 geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray70", size = 0.25) +
 geom_point(data = map_df, aes(x = LON, y = LAT), fill = "blue", shape = 21, size = 2, color = "black", stroke = 0.2) +
 facet_wrap(~unit, ncol = 4) +
 labs(title = "Fig. S1. The geographic distribution of groups used in G-PhoCS analyses.") +
 theme(plot.margin = margin(20, 20, 20, 20),
       plot.title = element_text(margin = margin(20, 0, 20, 0)),
       strip.text = element_text(face = "italic", margin = margin(0, 0, 5, 0)))

## Save:
ggsave(plot = map, filename = "supplementary_files/Fig_S1_G-PhoCS_sampling_maps.pdf", dpi = 300, units = "in", width = 8.5, height = 11)

## End of script.
