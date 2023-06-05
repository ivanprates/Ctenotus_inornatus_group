###############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### May 2023.
### The goals of this script are:
### To organize sample information for Table S1 and Table S2.
### To organize data of newly generated mitochondrial sequences to upload to GenBank.
### To plot maps for Fig. S1 and Fig. S2.

## PART 1: Setting up ----

## Packages:
library(ggrepel)
library(grid)
library(gridExtra)
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

## PART 2: Table S1: ddRAD samples ----

## Read tree:
ddRAD_tree <- read.nexus("RAxML/n242_m70_bs500/RAxML_bipartitions.n242_m70_bs500_7")

## Sample info:
ddRAD_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)
ddRAD_info <- ddRAD_info[ddRAD_info$SAMPLE_ID %in% ddRAD_tree$tip.label, ]

## Exclude mixed-up sample:
ddRAD_info <- ddRAD_info[ddRAD_info$SAMPLE_ID %in% setdiff(ddRAD_info$SAMPLE_ID, "WAMR_163003_Ct_duri"), ]

## Labels:
ddRAD_info$label <- ddRAD_info$VOUCHER
unvouchered <- ddRAD_info$SAMPLE_ID[ddRAD_info$VOUCHER == "no_voucher"]
for (sample in unvouchered) { ddRAD_info$label[ddRAD_info$SAMPLE_ID == sample] <- ddRAD_info$REGO[ddRAD_info$SAMPLE_ID == sample] }
unvouchered <- ddRAD_info$SAMPLE_ID[is.na(ddRAD_info$label)]
for (sample in unvouchered) { ddRAD_info$label[ddRAD_info$SAMPLE_ID == sample] <- sample }
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "NA_", replacement = "")
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "_Ct_.+", replacement = "")
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "SAMR", replacement = "SAMAR")
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "([A-Z]+)([0-9]+)", replacement = "\\1_\\2")
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "(R)_([0-9]+)", replacement = "_\\1\\2")
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "__", replacement = "_")
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = " ", replacement = "_")

## Subset columns:
ddRAD_info <- ddRAD_info[c("label", "SAMPLE_ID", "TISSUE_ID", "ALTERNATIVE_ID", "GENUS", "COLLECTOR_SP",
                          "CANDIDATE_SP_III", "STATE", "LOCATION", "LAT", "LON")]

## Minor changes:
ddRAD_info$TISSUE_ID[ddRAD_info$SAMPLE_ID == "NTMR_22167_Ct_bore" & !is.na(ddRAD_info$SAMPLE_ID)] <- "ABTC29698"
ddRAD_info$ALTERNATIVE_ID[ddRAD_info$SAMPLE_ID == "NTMR_22167_Ct_bore" & !is.na(ddRAD_info$SAMPLE_ID)] <- "BOREA29698"
ddRAD_info$label <- gsub(ddRAD_info$label, pattern = "_", replacement = " ")
ddRAD_info$TISSUE_ID <- gsub(ddRAD_info$TISSUE_ID, pattern = "_", replacement = "")
ddRAD_info$LOCATION <- gsub(ddRAD_info$LOCATION, pattern = "- ", replacement = "")
ddRAD_info$LOCATION <- gsub(ddRAD_info$LOCATION, pattern = " Western Australia", replacement = "")
ddRAD_info$LOCATION <- gsub(ddRAD_info$LOCATION, pattern = " NSW\\.", replacement = "")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "([0-9]+)km ", replacement = "\\1 km ")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "([0-9]+)k ", replacement = "\\1 km ")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "([0-9]+)K ", replacement = "\\1 km ")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "H/S\\.", replacement = "HS")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "R\\.", replacement = "R")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "SA\\.", replacement = "SA")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "NP\\.", replacement = "NP")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, "KIMBERLEY RESEARCH STATIO", "KIMBERLEY RESEARCH STATION")
ddRAD_info$LOCATION <- gsub(ddRAD_info$LOCATION, pattern = "  ", replacement = " ")

## Capitalization:
replace <- ddRAD_info$LOCATION[grep(x = ddRAD_info$LOCATION, pattern = "[A-Z][A-Z][A-Z]")]
ignore <- ddRAD_info$LOCATION[grep(x = ddRAD_info$LOCATION, pattern = " [A-Z][A-Z][A-Z] ")]
replace <- setdiff(replace, ignore)
for (label in replace) {
  ddRAD_info$LOCATION[ddRAD_info$LOCATION == label] <- str_to_title(label)
}

## Fixing some problems from capitalization:
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "([0-9]+) KM ", replacement = "\\1 km ")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = "Bdrs11", replacement = "BDRS11")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = " Hs", replacement = " HS")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = " Np", replacement = " NP")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, pattern = " Pp", replacement = " PP")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, "GaTE", "Gate")
ddRAD_info$LOCATION <- gsub(x = ddRAD_info$LOCATION, "Of Windorah", "of Windorah")

## NAs:
ddRAD_info$TISSUE_ID[is.na(ddRAD_info$TISSUE_ID)] <- ""
ddRAD_info$ALTERNATIVE_ID[is.na(ddRAD_info$ALTERNATIVE_ID)] <- ""

## Arrange:
ddRAD_info <- arrange(ddRAD_info, by = label)

## SRA accessions:
ddRAD_info$Accession <- ""
#in_lib_names <- ddRAD_info$SAMPLE_ID[ddRAD_info$SAMPLE_ID %in% SRA$Library.Name]
#ni_lib_names <- setdiff(ddRAD_info$SAMPLE_ID, in_lib_names)
SRA <- read.csv(file = "sample_information/SRA_Accessions_Ctenotus_2023-05.csv", header = TRUE)
for (sample in SRA$Library.Name) {
  ddRAD_info$Accession[ddRAD_info$SAMPLE_ID == sample] <- SRA$Experiment.Accession[SRA$Library.Name == sample]
}
ddRAD_info$Accession[ddRAD_info$SAMPLE_ID == "WAMR_157958_Le_ips"] <- "SRX11811099"
ddRAD_info$Accession[ddRAD_info$SAMPLE_ID == "WAMR_135152_Le_ips"] <- "SRX11811083"
ddRAD_info$Accession[ddRAD_info$SAMPLE_ID == "WAMR_111809_Le_bipe"] <- "SRX11811067"
ddRAD_info$Accession[ddRAD_info$SAMPLE_ID == "SAMAR_45424_Le_bipe"] <- "SRX11811124"
  
## Updating OTU labels:
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "inornatus-E", replacement = "eutaenius")
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "spaldingi-TE", replacement = "robustus-TE")
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "spaldingi-NW", replacement = "robustus-NW")
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "superciliaris-K", replacement = "mastigura")

## Putative taxa:
ddRAD_info$Putative_taxon <- ddRAD_info$CANDIDATE_SP_III
ddRAD_info$Putative_taxon[ddRAD_info$Putative_taxon %in% c("inornatus-N", "inornatus-S")] <- "inornatus"
ddRAD_info$Putative_taxon[ddRAD_info$Putative_taxon %in% c("robustus-NW", "robustus-TE")] <- "robustus"
ddRAD_info$Putative_taxon[ddRAD_info$Putative_taxon %in% c("spaldingi-CY", "spaldingi-NE", "spaldingi-S")] <- "spaldingi"
ddRAD_info$Putative_taxon[ddRAD_info$Putative_taxon %in% c("superciliaris-E (N)", "superciliaris-E (S)", "superciliaris-W")] <- "superciliaris"
ddRAD_info$Putative_taxon[ddRAD_info$Putative_taxon == "eutaenius"] <- "cf. eutaenius"
ddRAD_info$Putative_taxon[ddRAD_info$Putative_taxon == "mastigura"] <- "cf. mastigura"

## Checks:
sort(unique(ddRAD_info$Putative_taxon))
sort(unique(ddRAD_info$CANDIDATE_SP_III))

## Change column names:
names(ddRAD_info)[names(ddRAD_info) == "label"] <- "Voucher or Tissue Sample"
names(ddRAD_info)[names(ddRAD_info) == "SAMPLE_ID"] <- "Sample ID"
names(ddRAD_info)[names(ddRAD_info) == "TISSUE_ID"] <- "Tissue ID"
names(ddRAD_info)[names(ddRAD_info) == "ALTERNATIVE_ID"] <- "Alternative ID"
names(ddRAD_info)[names(ddRAD_info) == "GENUS"] <- "Genus"
names(ddRAD_info)[names(ddRAD_info) == "COLLECTOR_SP"] <- "Original Taxon"
names(ddRAD_info)[names(ddRAD_info) == "Putative_taxon"] <- "Putative Taxon"
names(ddRAD_info)[names(ddRAD_info) == "CANDIDATE_SP_III"] <- "Putative Unit"
names(ddRAD_info)[names(ddRAD_info) == "LAT"] <- "Latitude"
names(ddRAD_info)[names(ddRAD_info) == "LON"] <- "Longitude"
names(ddRAD_info)[names(ddRAD_info) == "LOCATION"] <- "Location"
names(ddRAD_info)[names(ddRAD_info) == "STATE"] <- "State or Territory"
names(ddRAD_info)[names(ddRAD_info) == "Accession"] <- "SRA Exp. Accession"

## Reorder columns:
ddRAD_info <- ddRAD_info[c(1:7, 13, 8:12)]

## SRA NAs:
ddRAD_info$"SRA Exp. Accession"[ddRAD_info$"SRA Exp. Accession" == ""] <- "NA"

## Save:
write.csv(x = ddRAD_info, file = "2nd_submission_MolEcol/supplementary_files_2nd_submission/Table_S1_ddRAD.csv", row.names = FALSE)

## And as pdf:
df <- ddRAD_info
maxrow <- 31
bs <- 11
npages <- ceiling(nrow(df)/maxrow)      
pdf(file = "2nd_submission_MolEcol/supplementary_files_2nd_submission/Table_S1_ddRAD.pdf", width = 20, height = 10)
idx <- seq(1, maxrow)  
grid.table(df[idx, ], rows = NULL, theme = ttheme_default(base_size = bs))  
for (i in 2:npages) {
  grid.newpage()
  if (i*maxrow <= nrow(df)) { idx = seq(1+((i-1)*maxrow), i*maxrow) }
  else { idx <- seq(1+((i-1)*maxrow), nrow(df)) }
  grid.table(df[idx, ], rows = NULL, theme = ttheme_default(base_size = bs))
}
dev.off()

## PART 3: cytb sample information ----

## Read tree:
mt_tree <- read.tree("cytb/RAxML_bipartitions.n485.tre")

## Sample info:
mt_info <- read.csv(file = "sample_information/raboskylab_SQL_individuals_all.csv", header = TRUE)
mt_info <- mt_info[mt_info$SAMPLE_ID %in% mt_tree$tip.label, ]

## Fixing some samples:
for (sample in c("NA_CCM3918_Ct_inor", "NA_CCM4044_Ct_robu", "NA_CCM3923_Ct_inor", "NA_CCM4016_Ct_quir")) {
  mt_info$LAT[mt_info$SAMPLE_ID == sample] <- -1*mt_info$LAT[mt_info$SAMPLE_ID == sample]
}
mt_info$LAT[mt_info$SAMPLE_ID == "NA_CCM3568_Ct_inor" & !is.na(mt_info$SAMPLE_ID)] <- -16.79282	
mt_info$LON[mt_info$SAMPLE_ID == "NA_CCM3568_Ct_inor" & !is.na(mt_info$SAMPLE_ID)] <- 137.24896
mt_info$LAT[mt_info$SAMPLE_ID == "NA_CCM3569_Ct_inor" & !is.na(mt_info$SAMPLE_ID)] <- -16.79801
mt_info$LON[mt_info$SAMPLE_ID == "NA_CCM3569_Ct_inor" & !is.na(mt_info$SAMPLE_ID)] <- 137.24995
mt_info$LOCATION[mt_info$SAMPLE_ID == "NA_CCM0044_Ct_euta" & !is.na(mt_info$SAMPLE_ID)] <- "The Vines Magnetic Island"
mt_info$LOCATION[mt_info$SAMPLE_ID == "NA_CCM4006_Ct_bore" & !is.na(mt_info$SAMPLE_ID)] <- "Garig Gunak Barlu National Park"
mt_info$REGO <- gsub(mt_info$REGO, pattern = "QM", replacement = "QMJ")
mt_info$VOUCHER <- gsub(mt_info$VOUCHER, pattern = "QM", replacement = "QMJ")
mt_info$SPECIES <- gsub(mt_info$SPECIES, pattern = "\\?", replacement = "")
mt_info$SPECIES <- gsub(mt_info$SPECIES, pattern = "[a-z]+_([a-z]+)", replacement = "\\1")
mt_info$TISSUE_ID[mt_info$SAMPLE_ID == "NTMR_22167_Ct_bore" & !is.na(mt_info$SAMPLE_ID)] <- "ABTC29698"
mt_info$ALTERNATIVE_ID[mt_info$SAMPLE_ID == "NTMR_22167_Ct_bore" & !is.na(mt_info$SAMPLE_ID)] <- "BOREA29698"
mt_info$REGO[mt_info$REGO == "wbj1027" & !is.na(mt_info$REGO)] <- "WBJ1027"
smiths <- mt_info$REGO[mt_info$LOCATION == "Smith's Lake" & !is.na(mt_info$LOCATION)]
for (sample in smiths) {
  mt_info$REGO[mt_info$REGO == sample] <- paste0("WAM", sample)
}

## Decimals:
mt_info$LAT <- sprintf("%.4f", mt_info$LAT)
mt_info$LON <- sprintf("%.4f", mt_info$LON)

## Change labels:
mt_info$label <- mt_info$VOUCHER
unvouchered <- mt_info$SAMPLE_ID[mt_info$VOUCHER == "no_voucher" & !is.na(mt_info$VOUCHER)]
for (sample in unvouchered) { mt_info$label[mt_info$SAMPLE_ID == sample] <- mt_info$REGO[mt_info$SAMPLE_ID == sample] }
unvouchered <- mt_info$SAMPLE_ID[is.na(mt_info$label)]
for (sample in unvouchered) { mt_info$label[mt_info$SAMPLE_ID == sample] <- sample }
mt_info$label <- gsub(mt_info$label, pattern = "NA_", replacement = "")
mt_info$label <- gsub(mt_info$label, pattern = "_Ct_.+", replacement = "")
mt_info$label <- gsub(mt_info$label, pattern = "SAMR", replacement = "SAMAR")
mt_info$label <- gsub(mt_info$label, pattern = "([A-Z]+)([0-9]+)", replacement = "\\1_\\2")
mt_info$label <- gsub(mt_info$label, pattern = "(R)_([0-9]+)", replacement = "_\\1\\2")
mt_info$label <- gsub(mt_info$label, pattern = "__", replacement = "_")
mt_info$label <- gsub(mt_info$label, pattern = " ", replacement = "_")

## Subset columns:
mt_info <- mt_info[c("label", "SAMPLE_ID", "TISSUE_ID", "ALTERNATIVE_ID", "GENUS", "SPECIES",
                     "STATE", "LOCATION", "LAT", "LON")]

## Minor changes:
mt_info$label <- gsub(mt_info$label, pattern = "_", replacement = " ")
mt_info$TISSUE_ID <- gsub(x = mt_info$TISSUE_ID, pattern = "_", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "~", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = ",", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = ";", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "- ", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = " Western Australia", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = " NSW\\.", replacement = "")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "([0-9]+)km ", replacement = "\\1 km ")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "([0-9]+)k ", replacement = "\\1 km ")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "([0-9]+)K ", replacement = "\\1 km ")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "([0-9]+)KM ", replacement = "\\1 km ")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "H/S\\.", replacement = "HS")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "R\\.", replacement = "R")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "SA\\.", replacement = "SA")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "NP\\.", replacement = "NP")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "SF\\.", replacement = "SF")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "Pk\\.", replacement = "Pk")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "Stn\\.", replacement = "Stn")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, "BARLEE RANGE NATURE RESER", "BARLEE RANGE NATURE RESERVE")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, "KIMBERLEY RESEARCH STATIO", "KIMBERLEY RESEARCH STATION")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = "  ", replacement = " ")

## Capitalization:
replace <- mt_info$LOCATION[grep(x = mt_info$LOCATION, pattern = "[A-Z][A-Z][A-Z]")]
ignore <- mt_info$LOCATION[grep(x = mt_info$LOCATION, pattern = " [A-Z][A-Z][A-Z] ")]
replace <- setdiff(replace, ignore)
for (label in replace) {
  mt_info$LOCATION[mt_info$LOCATION == label] <- str_to_title(label)
}

## Fixing some problems from capitalization:
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = " Hs", replacement = " HS")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, pattern = " Np", replacement = " NP")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, "GaTE", "Gate")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, "Of Windorah", "of Windorah")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, "Of Winton", "of Winton")
mt_info$LOCATION <- gsub(x = mt_info$LOCATION, "MOUNT TOM PRICE MINE", "Mount Tom Price Mine")

## NAs:
mt_info$TISSUE_ID[is.na(mt_info$TISSUE_ID)] <- ""
mt_info$ALTERNATIVE_ID[is.na(mt_info$ALTERNATIVE_ID)] <- ""

## Arrange:
mt_info <- arrange(mt_info, by = label)

## GenBank accession numbers of all Ctenotus cytb sequences in GenBank as per 2022-10-04:
accession1 <- read.table(file = "sample_information/Ctenotus_cytb_NCBI_2022-10-04_samples.txt", header = TRUE)

## Filter to the inornatus species group:
accession2 <- merge(mt_info, accession1, by.x = "ALTERNATIVE_ID", by.y = "Isolate", all.x = TRUE)
accession2$label <- gsub(accession2$label, pattern = " ", replacement = "")

## Arrange:
accession2 <- arrange(accession2, label)

## For those not identified by isolate, try voucher or other label:
noGB <- accession2[is.na(accession2$Accession), ]
noGB <- noGB[c("label")]
noGB <- merge(noGB, accession1, by.x = "label", by.y = "Voucher", all.x = TRUE)
noGB <- noGB[!is.na(noGB$Accession), ]

## Add Accessions:
for (sample in noGB$label) {
  #sample <- noGB$label[1]
  accession2$Accession[accession2$label == sample & !is.na(accession2$label)] <- noGB$Accession[noGB$label == sample]
}

## Adding existing accessions for outgroups:
accession2$Accession[accession2$label == "NTMR22191"] <- "KJ505005"
accession2$Accession[accession2$label == "UMMZ242633"] <- "ON036022"
accession2$Accession[accession2$label == "UMMZ242639"] <- "ON036025"

## Change column names:
accession3 <- accession2
names(accession3)[names(accession3) == "label"] <- "Voucher or Tissue Sample"
names(accession3)[names(accession3) == "SAMPLE_ID"] <- "Sample ID"
names(accession3)[names(accession3) == "TISSUE_ID"] <- "Tissue ID"
names(accession3)[names(accession3) == "ALTERNATIVE_ID"] <- "Alternative ID"
names(accession3)[names(accession3) == "GENUS"] <- "Genus"
names(accession3)[names(accession3) == "SPECIES"] <- "Original Taxon"
names(accession3)[names(accession3) == "CANDIDATE_SP_III"] <- "Putative Unit"
names(accession3)[names(accession3) == "LAT"] <- "Latitude"
names(accession3)[names(accession3) == "LON"] <- "Longitude"
names(accession3)[names(accession3) == "LOCATION"] <- "Location"
names(accession3)[names(accession3) == "STATE"] <- "State or Territory"

## Columns:
accession3 <- accession3[setdiff(names(accession3), c("Genus", "GenBank_sp", "Voucher"))]
accession3 <- accession3[c(2, 3, 4, 1, 5:10)]

## PART 4: Tentative OTU and taxon assignment of mtDNA samples ----

## Samples with ddRAD data:
ddRAD_info <- read.csv(file = "sample_information/ddRAD_sample_information_organized_n1315_2022-04-13.csv", header = TRUE)

## Updating OTU labels:
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "inornatus-E", replacement = "eutaenius")
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "spaldingi-TE", replacement = "robustus-TE")
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "spaldingi-NW", replacement = "robustus-NW")
ddRAD_info$CANDIDATE_SP_III <- gsub(ddRAD_info$CANDIDATE_SP_III, pattern = "superciliaris-K", replacement = "mastigura")

## Populate with taxon info:
ddRAD_info <- ddRAD_info[ddRAD_info$SAMPLE_ID %in% accession3$"Sample ID", ]
accession4 <- accession3
accession4$OTU <- "No nuclear data"
for (sample in ddRAD_info$SAMPLE_ID) {
  accession4$OTU[accession4$"Sample ID" ==  sample] <- ddRAD_info$CANDIDATE_SP_III[ddRAD_info$SAMPLE_ID == sample]
}

## Rename terminals in tree:
mt_tree2 <- mt_tree
for (sample in accession4$"Sample ID") {
  mt_tree2$tip.label[mt_tree2$tip.label ==  sample] <- accession4$"Voucher or Tissue Sample"[accession4$"Sample ID" == sample]
}

## Now, to add putative OTU for samples with no ddRAD data, let's use mitochondrial tree:
## Define clades:
list_span <- list()
list_span$"inornatus-S" <- c("WAMR140720", "SAMAR44367")
list_span$"lateralis" <- c("SAMAR55259", "SAMAR42768")
list_span$"robustus-NW" <- c("WAMR146354", "CCM4044")
list_span$"robustus-TE" <- c("NTMR20378", "NTMR23946")
list_span$"superciliaris-E (N)" <- c("CCM2365", "CCM2292")
list_span$"superciliaris-W" <- c("WAMR153904", "WAMR135692")
list_span$"spaldingi-CY" <- c("ANWCR05242", "ANWCR05240")
list_span$"spaldingi-NE" <- c("ANWCR05369", "ERPQ31340")
list_span$"spaldingi-S" <- c("WAMR130175", "SAMAR46204")

## Find the samples within each clade:
no_ddRAD <- accession4$"Voucher or Tissue Sample"[accession4$OTU == "No nuclear data"]
for (clade in names(list_span)) { 
  node <- findMRCA(tree = mt_tree2, type = "node", tips = list_span[[clade]])
  descendants <- getDescendants(mt_tree2, node = node)
  descendants <- descendants[descendants <= length(mt_tree2$tip.label)]
  descendants <- sort(mt_tree2$tip.label[descendants])
  descendants <- descendants[descendants %in% no_ddRAD]
  for (sample in descendants) {
    accession4$OTU[accession4$"Voucher or Tissue Sample" ==  sample] <- clade
  }
}

## Some manual additions for samples not automated:
accession4$OTU[accession4$"Original Taxon" ==  "burbidgei"] <- "burbidgei"
#accession4$OTU[accession4$"Voucher or Tissue Sample" ==  "ERPQ31348"] <- "nullum"
accession4$OTU[accession4$"Voucher or Tissue Sample" ==  "ERPQ31348"] <- "inornatus-S"
accession4$OTU[accession4$"Voucher or Tissue Sample" %in% c("ABTC21767", "CCM2831", "CCM2855", "CCM3568", "CCM3569", "CCM3717", "CCM4018", "CCM4039", "WAMR126034")] <- "inornatus-N"
accession4$OTU[accession4$"Voucher or Tissue Sample" %in% c("CCM4016", "ERPQ31349", "QMJ87459", "ANWCR05463")] <- "spaldingi-NE"
accession4$OTU[accession4$"Voucher or Tissue Sample" %in% c("SAMAR55745", "SAMAR55880")] <- "spaldingi-S"
accession4$OTU[accession4$"Voucher or Tissue Sample" %in% c("CCM3004", "WAMR108775", "WAMR108699", "WAMR108746", "WAMR108776")] <- "superciliaris-W"
accession4$OTU[accession4$"Voucher or Tissue Sample" %in% c("CCM2023", "CCM2169", "CCM2175", "CCM2196", "CCM2197", 
                                                            "NTMR20243", "NTMR22185", "SAMAR36246", "SAMAR38760", "SAMAR38761", "SAMAR38762", "SAMAR38766",
                                                            "WAMR129245", "WAMR129291")] <- "superciliaris-E (S)"

## Taxon assignments:
accession4$Putative_taxon <- accession4$OTU
accession4$Putative_taxon[accession4$OTU == "eutaenius"] <- "cf. eutaenius"
accession4$Putative_taxon[accession4$OTU %in% c("inornatus-N", "inornatus-S")] <- "inornatus"
accession4$Putative_taxon[accession4$OTU %in% c("spaldingi-CY", "spaldingi-NE", "spaldingi-S")] <- "spaldingi"
accession4$Putative_taxon[accession4$OTU %in% c("robustus-NW", "robustus-TE")] <- "robustus"
accession4$Putative_taxon[accession4$OTU == "mastigura"] <- "cf. mastigura"
accession4$Putative_taxon[accession4$OTU %in% c("superciliaris-E (N)", "superciliaris-E (S)", "superciliaris-W")] <- "superciliaris"

## Reorder columns:
accession4 <- accession4[c(1:5, 11, 12, 6:10)]

## PART 5: Testing assignments based on maps ----

## Data to map:
accession5 <- accession4
names(accession5)[1] <- "label"
accession5$Latitude <- as.numeric(accession5$Latitude)
accession5$Longitude <- as.numeric(accession5$Longitude)
accession5 <- accession5[accession5$OTU %in% setdiff(accession5$OTU, c("atlas", "australis", "essingtonii", "leonhardii",
                                                                       "nigrilineatus", "pantherinus", "schomburgkii", "taeniolatus")), ]

## Australia map:
AUS_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), admin == "Australia")
AUS_map <- crop(as_Spatial(AUS_map), extent(matrix(c(113, -10, 154, -39.5), nrow = 2))) ## Extent.

## Plot map with ggplot:
map <- ggplot(data = accession5) +
        
  ## Adding baseline map:
  geom_sf(data = st_as_sf(AUS_map), fill = "white", color = "gray70", size = 0.15) +
        
  ## Adding ID labels on top of map:
  #geom_text_repel(aes(x = Longitude, y = Latitude, label = label), color = "black", segment.color = "black", size = 3, segment.size = 0.25) +
        
  ## Adding lat-longs for sampled individuals:
  geom_point(aes(x = Longitude, y = Latitude), fill = "blue", shape = 21, alpha = 0.5, size = 1.25, color = "black", stroke = 0.2) +
     
  ## Setting a bunch of aesthetic parameters:
  theme_void() +
        
  ## Setting some other elements:
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        strip.text = element_text(size = 5, face = "italic", margin = margin(t = 0, r = 0, b = 2, l = 0)),
        strip.background = element_blank(),
        plot.title = element_text(size = 6, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        plot.subtitle = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10))

## Check:
map_OTUs <- map + facet_wrap(~OTU, ncol = 3) + ggtitle(label = "Figure S1. Geographic distribution of mitochondrial samples partitioned\nby the corresponding nuclear operational taxonomic units.")
map_taxa <- map + facet_wrap(~Putative_taxon, ncol = 3) + ggtitle(label = "Figure S2. Geographic distribution of mitochondrial samples partitioned\nby the corresponding putative taxa.")

## Save:
ggsave(plot = map_OTUs, filename = "2nd_submission_MolEcol/supplementary_files_2nd_submission/Fig_S1_maps_cytb_samples_OTUs.pdf", width = 8.5, height = 11, units = "cm", limitsize = FALSE, dpi = 500)
ggsave(plot = map_taxa, filename = "2nd_submission_MolEcol/supplementary_files_2nd_submission/Fig_S2_maps_cytb_samples_taxa.pdf", width = 8.5, height = 11, units = "cm", limitsize = FALSE, dpi = 500)

## PART 6: Preparing GenBank submission ----

## Let's now create a source modifier file for the samples to be submitted to GenBank.
noGB <- accession5[is.na(accession5$Accession), ]

## Edits to match NCBI requirements:
noGB$label <- gsub(noGB$label, pattern = "R", replacement = ":R:")
noGB$label <- gsub(noGB$label, pattern = "CCM", replacement = "CCM:")
noGB$label <- gsub(noGB$label, pattern = "CUMV", replacement = "CUMV:")
noGB$label <- gsub(noGB$label, pattern = "PMO", replacement = "PMO:")
noGB$label <- gsub(noGB$label, pattern = "UMMZ", replacement = "UMMZ:")
noGB$label <- gsub(noGB$label, pattern = "DL:R", replacement = "DLR")
noGB$Latitude <- gsub(noGB$Latitude, pattern = "-", replacement = "")
noGB$Latitude <- paste0(noGB$Latitude, " S")
noGB$Longitude <- paste0(noGB$Longitude, " E")
noGB$lat_lon <- paste0(noGB$Latitude, " ", noGB$Longitude)
noGB$country <- paste0("Australia: ", noGB$"State or Territory", ", ", noGB$Location)

## Some taxonomic notes:
noGB$notes <- ""
noGB$notes[noGB$Putative_taxon == "cf. eutaenius"] <- "Appears to correspond to Ctenotus eutaenius."
noGB$notes[noGB$Putative_taxon == "cf. mastigura"] <- "Might correspond to Ctenotus mastigura."

## Columns to keep and their names:
noGB <- noGB[c("Sample ID", "country", "lat_lon", "label", "notes")]
names(noGB) <- c("sequence_ID", "country", "lat_lon", "specimen_voucher", "Note")

## Save:
#write.table(x = noGB, file = "GenBank/source_modifiers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## Organism file for GenBank:
org_file <- accession6[c("Sample ID", "Putative Taxon")]
names(org_file) <- c("sequence_ID", "organism")
org_file <- org_file[org_file$sequence_ID %in% noGB$sequence_ID, ]
org_file$organism[org_file$organism == "cf. eutaenius"] <- "sp. 1"
org_file$organism[org_file$organism == "cf. mastigura"] <- "sp. 2"
org_file$organism <- paste0("Ctenotus ", org_file$organism)

## Save:
#write.table(x = org_file, file = "GenBank/organism.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## Alignment:
alignment <- read.csv("cytb/alignments/inornatusgr_cytb_n485.csv", header = FALSE)
names(alignment) <- c("SAMPLE_ID", "SEQ")
alignment <- alignment[alignment$SAMPLE_ID %in% noGB$"sequence_ID", ] ; dim(alignment)

## Label column:
alignment$label <- paste0(">", alignment$SAMPLE_ID)

## For samples with no site info, label equal ID:
for (sample in alignment$SAMPLE_ID[!is.na(alignment$LOC_fixed)]) {
  alignment$seqlabel[alignment$SAMPLE_ID == sample] <- paste0(">", alignment$SAMPLE_ID[alignment$SAMPLE_ID == sample], "_", alignment$LOC_fixed[alignment$SAMPLE_ID == sample])
}

## Simplify:
alignment <- alignment[c("label", "SEQ")]

## Export alignment:
#write.table(alignment, file = "GenBank/alignment.fasta", 
#            row.names = FALSE, quote = FALSE, sep = "\n", eol = " \n", col.names = FALSE)

## PART 7: Final Supplementary Table S2 ----

## Final edits for Table S1:
accession6 <- accession5
names(accession6)[1] <- "Voucher or Tissue Sample"
names(accession6)[7] <- "Putative Taxon"
names(accession6)[12] <- "GenBank Accession"

## Add accessions of new sequences sent by NCBI:
GB_new <- read.table(file = "GenBank/seqids.txt", header = FALSE)
samples_GB_new <- accession6$"Sample ID"[accession6$"Sample ID" %in% GB_new$V2]
for (sample in samples_GB_new) {
  accession6$"GenBank Accession"[accession6$"Sample ID" == sample] <- GB_new$V3[GB_new$V2 == sample]
}

## Checks:
sort(unique(accession6$"Putative Taxon"))
sort(unique(accession6$OTU))

## Finally, save as csv:
write.csv(x = accession6, file = "2nd_submission_MolEcol/supplementary_files_2nd_submission/Table_S2_cytb.csv", row.names = FALSE)

## And as pdf:
df <- accession6
maxrow <- 30
bs <- 11
npages <- ceiling(nrow(df)/maxrow)      
pdf(file = "2nd_submission_MolEcol/supplementary_files_2nd_submission/Table_S2_cytb.pdf", width = 20, height = 10)
idx <- seq(1, maxrow)  
grid.table(df[idx, ], rows = NULL, theme = ttheme_default(base_size = bs))  
for (i in 2:npages) {
  grid.newpage()
  if (i*maxrow <= nrow(df)) { idx = seq(1+((i-1)*maxrow), i*maxrow) }
  else { idx <- seq(1+((i-1)*maxrow), nrow(df)) }
  grid.table(df[idx, ], rows = NULL, theme = ttheme_default(base_size = bs))
}
dev.off()

## End of script.
