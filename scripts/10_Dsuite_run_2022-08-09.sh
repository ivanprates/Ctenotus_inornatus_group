#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of the script are:
### To run Dsuite to estimate genetic introgression in the Ctenotus inornatus group.
### Good resources:
### https://github.com/millanek/Dsuite
### https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data

## Project:
project=Ctenotus_inornatus_group

## Data:
dataset=inornatusgr_R1_c90_n281_m70
#dataset=inornatusgr_R1_c90_n264_m65
#dataset=inornatusgr_R1_c90_n268_m80

## Taxonomic scheme:
#scheme=CANDIDATE_SP_I
scheme=CANDIDATE_SP_III

## Create and change directory:
cd ~/Dropbox/Science/MYPAPERS_ongoing/${project}/Dsuite/${dataset}/

## Run program Dtrios:
#echo "Now running Dtrios!"
~/Dsuite/Build/Dsuite Dtrios -c -t ${scheme}_tree.tre -o ${scheme} ~/Dropbox/Science/MYPAPERS_ongoing/${project}/ipyrad/${dataset}_outfiles/snps.vcf ${scheme}_assignments.txt

## Plot heatmaps:
echo "Now plotting heatmaps!"
ruby ~/Dsuite/plot_d.rb ${scheme}_BBAA.txt ${scheme}_plot_order.txt 1 ${scheme}_heatmap_Dstats_BBAA.svg
ruby ~/Dsuite/plot_f4ratio.rb ${scheme}_BBAA.txt ${scheme}_plot_order.txt 1 ${scheme}_heatmap_f4-ratio_BBAA.svg

## Assess gene flow events with program Fbranch:
echo "Now running Fbranch!"
#~/Dsuite/Build/Dsuite Fbranch ${scheme}_tree.tre ${scheme}_tree.txt -p 0.01 > ${scheme}_Fbranch_p0.01.txt
~/Dsuite/Build/Dsuite Fbranch ${scheme}_tree.tre ${scheme}_tree.txt -p 0.05 > ${scheme}_Fbranch_p0.05.txt

## Plot Fbranch results:
echo "Now plotting Fbranch results!"
#python ~/Dsuite/utils/dtools.py ${scheme}_Fbranch_p0.01.txt ${scheme}_tree.tre
#mv fbranch.png ${scheme}_fbranch_p0.01.png
#mv fbranch.svg ${scheme}_fbranch_p0.01.svg
python ~/Dsuite/utils/dtools.py ${scheme}_Fbranch_p0.05.txt ${scheme}_tree.tre
mv fbranch.png ${scheme}_fbranch_p0.05.png
mv fbranch.svg ${scheme}_fbranch_p0.05.svg

echo "All done."

## End of script.
