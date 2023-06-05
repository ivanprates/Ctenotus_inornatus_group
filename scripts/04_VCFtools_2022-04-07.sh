#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To filter the SNP data using VCFtools.
### To change data format from vcf to geno (012) for sNMF analyses.
### Good resource: https://www.ddocent.com/filtering/

## Change folder:
#cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/ipyrad/mast-supe_n46_m50_outfiles
#cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/ipyrad/euta-inor-hele-late_n103_m50_outfiles
#cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/ipyrad/rima-robu-spal_n73_m50_outfiles
cd ~/Dropbox/Science/MYPAPERS_ongoing/Ctenotus_inornatus_group/ipyrad/inornatusgr_R1_c90_n264_m65_outfiles

## Apply MAF filter, residual non-biallelic sites, and if needed, remove samples:
#vcftools --vcf *.vcf --remove ~/Dropbox/Science/MYPAPERS_ongoing/2021_inornatus_group/ipyrad/remove_from_vcf.txt --recode --out snps --maf 0.05 --max-alleles 2
vcftools --vcf *.vcf --recode --out snps --maf 0.05 --max-alleles 2

## Renaming resulting vcf file:
mv snps.recode.vcf snps.vcf

## Now extracting one SNP per locus:
cat snps.vcf | grep "#" > usnps.vcf ## Extract the headers (whose line starts with #), save in different file.
cat snps.vcf | grep -v "#" | sort -u -k1,1 >> usnps.vcf ## Sort by RAD locus name and extract unique SNP per locus; concatenate with that file created.

## Saving in 012 format:
vcftools --vcf snps.vcf --out snps --012
vcftools --vcf usnps.vcf --out usnps --012

## Replacing -1 with 9 to indicate missing sites:
cat snps.012 | sed -e 's/-1/9/g' snps.012 > snps_tmp.012 ## Replace and save in temporary file.
cat usnps.012 | sed -e 's/-1/9/g' usnps.012 > usnps_tmp.012 ## Replace and save in temporary file.

## Lastly, renaming temporary file as a final 012 file:
mv snps_tmp.012 snps.012
mv usnps_tmp.012 usnps.012

## End of script.
