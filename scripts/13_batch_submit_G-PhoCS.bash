#!/bin/bash
###############################################
### Bash script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### August 2022.
### The goals of this script are:
### To create and submit multiple G-PhoCS analyses.

## Changing and submiting multiple jobs:
for control in *.ctl
do

## Set up second run:
#sed -i 's/_run1/_run2/' *.ctl 
#sed -i 's/12345/22345/' *.ctl 

## Status:
echo "Now doing ${control}!"

## Create temporary job file:
cp G-PhoCS.job ${control}.tjob

## Replace target control file in temporary job file:
sed -i "s/control.ctl/${control}/" ${control}.tjob

## Submit job:
sbatch ${control}.tjob

## Remove temporary job file:
rm ${control}.tjob
done

## Check:
squeue -u pratesi

## End of script.
~                            
