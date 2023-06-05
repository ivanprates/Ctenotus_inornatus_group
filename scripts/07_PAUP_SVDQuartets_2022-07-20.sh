#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To run coalescent-based phylogenetic analyses using SVD Quartets as implemented in PAUP.
### Resources:
### http://www.phylosolutions.com/tutorials/ssb2018/svdquartets-tutorial.html
### https://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/svdquartets-astral-activity/
### https://phylosolutions.com/tutorials/svdq-qage/svdq-qage-tutorial.html

## Folder:
cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/paup_SVDQuartets/

## Read file and run analysis:
## Species assignments added by hand into file.
paup inornatusgr_R1_c90_n242_m70.nex -L paup.log
#paup inornatusgr_R1_c90_n242_m70_notaxa.nex -L paup_notaxa.log
#paup inornatusgr_R1_c90_n242_m70_brac_fall_hele_seve.nex -L paup_brac_fall_hele_seve.log
#paup inornatusgr_R1_c90_n242_m70_brac_fall_hele_seve_no_misID.nex -L paup_brac_fall_hele_seve_no_misID.log

## Check tree:
#figtree inornatusgr_n242_m70_paup.tre

## End of script.
