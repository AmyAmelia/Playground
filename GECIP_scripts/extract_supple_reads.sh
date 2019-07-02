#! bin/bash
# Curated 14/June/2019
#############################################
# This script loops through a list of bams to identify supplementary alignments (-f 2048)
# removing PCR and optical duplicated (-F 1024)
#############################################
#############################################
# a list of variables to be updated.
# input text file listing the full path the bam files to be analysised
input=Desktop/bamlist.txt
# path and name of file to output the results to.
output=Desktop/test_roi3.txt
# Region of intrest to be investigated. If not required remove $ROI from script below.
ROI="20:15599975-16599974"

# load GeL modules 
# script
while read bam;
do echo $bam >> $output;
samtools view -u -F 1024 $bam $ROI| samtools view -f 2048 >> $output;
done < $input