#!/bin/bash -e

# script modified from original
#13-10-2016 #protocol by Catherine Collins #scripting by Sophia Cameron-Christie
#29-05-2020 # edited by Anna Gosling
#22-06-2020 # edited by Catherine to run on NeSI. All programs and parameters remain the same.
#31-07-2020 #edited by Catherine to map fastq reads from shotgun sequencing libraries from Wairau Bar lesion samples to human genome build hg37
# 31-05-2021 # edited by Hugh Cross to import variables from file, and allow for multiple adapters
# 23-07-2021 # split into two scripts by Hugh Cross: one for trimming and one for mapping


#variables

#directory with the raw fastq 
datadir=$(grep -v '#' $1 | grep 'datadir' | cut -f 2 -d '=')
#suffix pattern of the first fastq file
fq1=$(grep -v '#' $1 | grep 'fq1' | cut -f 2 -d '=')
#suffix pattern of the second fastq file  
fq2=$(grep -v '#' $1 | grep 'fq2' | cut -f 2 -d '=')
#fastq files zipped? "yes" or "no"
zipped_fastqs=$(grep -v '#' $1 | grep 'zipped_fastqs' | cut -f 2 -d '=')
#list of samples (separated by spaces, with one pair of quotes around all)
samplist=$(grep -v '#' $1 | grep 'samplist' | cut -f 2 -d '=')
# [Optional] in case using different adapters, otherwise will use default
adapt1=$(grep -v '#' $1 | grep 'adapter1' | cut -f 2 -d '=')
adapt2=$(grep -v '#' $1 | grep 'adapter2' | cut -f 2 -d '=')

echo 'variables read, proceeding' 

echo 'run on ' date

datecode=$(date +"%Y%m%d-%H%M")
logfilename='adapter_removal_run_'${datecode}.log

echo 'adapter removal run' $datecode > $logfilename
date +"%Y %m %d-%H:%M" >> $logfilename
echo ' ' >> $logfilename

###########################################################################################################################################################################

# load modules
module load AdapterRemoval/2.3.1-GCCcore-7.4.0

##########################################################################################################################################################################

# check for input adapters, or go with defaults
if [ -z "$adapt1" ]; then
adapter1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
else
adapter1=$adapt1
fi

if [ -z "$adapt2" ]; then
adapter2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
else
adapter2=$adapt2
fi

##########################################################################################################################################################################

for samp in $samplist
do

echo 'Sample: ' $samp >> $logfilename
#mkdir for sample
mkdir ${samp}_trimmed
cd ${samp}_trimmed

###################
##adaptor removal##
###################

if [ $zipped_fastqs = "yes" ]; then

AdapterRemoval \
--file1 ${datadir}/${samp}${fq1} \
--file2 ${datadir}/${samp}${fq2}  \
--gzip \
--adapter1 ${adapter1} \
--adapter2 ${adapter2} \
--collapse \
--trimns  \
--trimqualities  \
--minlength 25  \
--mm 3  \
--minquality 20  \
--basename $samp

else

AdapterRemoval \
--file1 ${datadir}/${samp}${fq1} \
--file2 ${datadir}/${samp}${fq2}  \
--adapter1 ${adapter1} \
--adapter2 ${adapter2} \
--collapse \
--trimns  \
--trimqualities  \
--minlength 25  \
--mm 3  \
--minquality 20  \
--basename $samp


fi

grep -A 14 'Trimming statistics' ${samp}.settings >> ../${logfilename}

cd ../

done
