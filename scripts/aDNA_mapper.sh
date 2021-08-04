#!/bin/bash -e

#13-10-2016 #protocol by Catherine Collins #scripting by Sophia Cameron-Christie
#29-05-2020 # edited by Anna Gosling
#22-06-2020 # edited by Catherine to run on NeSI. All programs and parameters remain the same.
#31-07-2020 #edited by Catherine to map fastq reads from shotgun sequencing libraries from Wairau Bar lesion samples to human genome build hg37
# 31-05-2021 # edited by Hugh Cross to import variables from file, and allow for multiple adapters

#variables
#name of the reference file
reffile=$(grep -v '#' $1 | grep 'reffile' | cut -f 2 -d '=')
#directory with the reference file
refdir=$(grep -v '#' $1 | grep 'refdir' | cut -f 2 -d '=')
#list of samples (separated by spaces, with one pair of quotes around all)
samplist=$(grep -v '#' $1 | grep 'samplist' | cut -f 2 -d '=')
# additional variables for mapping
trim_dir=$(grep -v '#' $1 | grep 'trim_dir' | cut -f 2 -d '=')
# trim folder su
trim_folder_suffix=$(grep -v '#' $1 | grep 'trim_folder_suffix' | cut -f 2 -d '=')
# target trimmed data file
target_data=$(grep -v '#' $1 | grep 'target_data' | cut -f 2 -d '=')
# mapping output folder
map_dir=$(grep -v '#' $1 | grep 'map_dir' | cut -f 2 -d '=')
# seed length for bwa aln
seedlength=$(grep -v '#' $1 | grep 'seedlength' | cut -f 2 -d '=')
# error rate for bwa aln
error_rate=$(grep -v '#' $1 | grep 'error_rate' | cut -f 2 -d '=')
# gap open for bwa aln (-o)
gapopens=$(grep -v '#' $1 | grep 'gapopens' | cut -f 2 -d '=')
# variable for num threads
num_threads=$(grep -v '#' $1 | grep 'num_threads' | cut -f 2 -d '=')
# check if will run mapdamage
mapdamage=$(grep -v '#' $1 | grep 'mapdamage' | cut -f 2 -d '=')
# use DeDup to do deduplication
dedup=$(grep -v '#' $1 | grep 'dedup' | cut -f 2 -d '=')

# ref file 
ref=$refdir$reffile

# check for user inputs on bwa options
if [ -z "$seedlength" ]; then
seedlen=1024
else
seedlen=$seedlength
fi

if [ -z "$error_rate" ]; then
error=0.03
else
error=$error_rate
fi

if [ -z "$gapopens" ]; then
gapo=2
else
gapo=$gapopens
fi

if [ -z "$num_threads" ]; then
threads=8
else
threads=$num_threads
fi


echo 'variables read, proceeding' 

###########################################################################################################################################################################

#variables that shouldn't need changing - but may be worth checking if versions have updated

module load BWA/0.7.17-gimkl-2017a
module load SAMtools/1.12-GCC-9.2.0
module load Java/15.0.2
#module load GATK/4.1.4.1-gimkl-2018b
#module load picard/2.21.8-Java-11.0.4
module load mapDamage

##########################################################################################################################################################################

#index the reference fasta file
if [ ! -e $refdir$reffile.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi

##

cd ${map_dir}

###############################################
## create log file ############################

datecode=$(date +"%Y%m%d-%H%M")
logfilename='adna_mapping_run_'${datecode}.log

echo 'aDNA mapping run' $datecode > $logfilename
date +"%Y %m %d-%H:%M" >> $logfilename
echo ' ' >> $logfilename


#####################################################


for samp in $samplist
do

#${trim_dir}/${samp}${trim_folder_suffix}/${samp}.${target_data}

##################
## Alignment	##
##################

##Find the SA coordinates of the input reads, using .collapsed reads only 

bwa aln -t $threads -n $error -o $gapo -l $seedlen $ref ${trim_dir}/${samp}${trim_folder_suffix}/${samp}.${target_data} > ${samp}.sai

bwa samse $ref ${samp}.sai ${trim_dir}/${samp}${trim_folder_suffix}/${samp}.${target_data} -f ${samp}.sam


#get the total reads for calculating endogenous percent

collapsedtotal=$(samtools view -c ${samp}.sam)

echo "$samp $collapsedtotal" >> total_reads.txt


#output the alignments as bam, ignoring alignment with quality score lower than 20
samtools view -@ $threads -b -S -q 20 ${samp}.sam > ${samp}.bam

if [[ "$dedup" = 'yes' ]]; then
echo 'DeDup will be used for deduplication'

java -Xmx16G -jar /nesi/project/uoo02328/programs/dedup/DeDup-0.12.8.jar \
  --input ${samp}.bam \
  --merged \
  --output $PWD
# remove unmapped reads
samtools view -@ $threads -b -F 0x0004 ${samp}_rmdup.bam.bam -o ${samp}_maponly.bam

else
echo 'Samtools will be used for deduplication'
# add ms/MC tags and sorting as pipe

samtools sort -@ $threads -n ${samp}.bam | samtools fixmate -@ $threads -m - ${samp}_fixmate.bam

samtools sort -@ $threads ${samp}_fixmate.bam | samtools markdup -@ $threads -r - ${samp}_markdup.bam
# remove unmapped reads
samtools view -@ $threads -b -F 0x0004 ${samp}_markdup.bam -o ${samp}_maponly.bam

fi

samtools index ${samp}_maponly.bam

#get the number of unique, mapped reads for calculating endogenous percent

collapsedtotal=$(samtools view -c ${samp}_maponly.bam)

echo "$samp $collapsedtotal" >> uniq_mapped_reads.txt

###############
## mapDamage ## 
###############

#quantifies DNA damage patterns aDNA NGS sequencing reads
# leaving out for now, as takes a long time
if [[ "$mapdamage" = 'yes' ]]; then
echo 'map damage will be run on sample'
mapDamage -i ${samp}_maponly.bam -r $ref --rescale 
else
echo 'map damage will not be run on sample'
fi 


done

# get rid of sam files, and move intermediate files to intermediate folder
# todo: make option to create intermediate_files folder, in case it is there
if [ ! -d "${map_dir}/intermediate_files" ]; then
echo "creating intermediate_files directory"
mkdir intermediate_files
fi 

mv *.sai intermediate_files
mv *.bam intermediate_files
mv intermediate_files/*_maponly.bam ./

rm *.sam

echo 'total reads mapped:' >> $logfilename
cat total_reads.txt >> $logfilename

echo ' ' >> $logfilename
echo 'unique reads mapped:' >> $logfilename
cat uniq_mapped_reads.txt >> $logfilename