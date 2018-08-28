#!/bin/bash

# MinIon pipeline for structural variant calls v1.0
# Lots of credit to David Coffey

# Start by reading inputs from the user - [/path/to/data] [ID] [SAMPLE] [LAST vs. minimap2]

# /path/to/data
echo "Enter /path/to/data"
read FAST5

# ID
echo "Enter data ID"
read ID

# SAMPLE
echo "Enter Sample Name"
read SAMPLE

#/path/to/radich home
export RADICH_HOME="/fh/fast/radich_j"

#Working folder for the session is created with the date and timestamp. All processes will be run in this folder
DATE=$(date)
DATE_FORMATTED="${DATE// /_}".
export WORKING_FOLDER="$RADICH_HOME/ngs/Jobs/$DATE_FORMATTED"
mkdir $WORKING_FOLDER

export FLOWCELL="FLO-MIN106"
export KIT="SQK-LSK108"

#/path/to/Fastq/Original (Albacore Folder)
export FASTQ="$WORKING_FOLDER/albacore-fastq"

#/path/to/Fastq/Demultiplexed (Porechop Folder)
export DEMULTIPLEXED_DIRECTORY="$WORKING_FOLDER/porechop-demux"

#/path/to/NanoPlot
export NANOPLOT="$WORKING_FOLDER/nanoplot"

#/path/to/ReferenceGenomes
export REFERENCE="/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta"

#/path/to/LastHumandb
export LASTAL_HUMANDB="$RADICH_HOME/ngs/Last-Genome-Data" 

#/path/to/Aligned
export ALIGNMENT_DIRECTORY="$WORKING_FOLDER/last-aligned"

#/path/to/NanoSV
export NANOSV_DIRECTORY="$WORKING_FOLDER/nanosv-vcf"

# Basecall with albacore
mkdir $FASTQ
module load Python/3.6.5-foss-2016b-fh3
read_fast5_basecaller.py \
--flowcell $FLOWCELL \
--kit $KIT \
--worker_threads 12 \
--input $FAST5 \
--output_format fastq \
--save_path $FASTQ \
--recursive \
--resume

# Merge all fastq files
#cat $FASTQ/workspace/pass/barcode*/*.fastq $FASTQ/workspace/pass/unclassified/*.fastq > $FASTQ/workspace/pass/merged.fastq
cat $FASTQ/workspace/pass/*.fastq > $FASTQ/workspace/pass/merged.fastq

# Demultiplex
mkdir $DEMULTIPLEXED_DIRECTORY
porechop \
--input $FASTQ/workspace/pass/merged.fastq \
--barcode_dir $DEMULTIPLEXED_DIRECTORY \
--format fastq \
--threads 12 \
-v 2

# Create NanoPlot figures
# NanoStats report included
mkdir $NANOPLOT
NanoPlot \
--fastq $DEMULTIPLEXED_DIRECTORY/none.fastq.gz \
--outdir $NANOPLOT \
--plots hex dot

#Branch here MiniMap2 or LAST

###############################

# #MiniMap2
# # Align to human genome using MiniMap2 
# mkdir $ALIGNMENT_DIRECTORY
# module load minimap2/2.10-foss-2016b
# minimap2 -a $REFERENCE/genome.fa $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > $ALIGNMENT_DIRECTORY/minimap.alignment.sam

# # Sort SAM file, convert to BAM and create index
# module load samtools
# samtools view -bS $ALIGNMENT_DIRECTORY/minimap.alignment.sam | samtools sort > $ALIGNMENT_DIRECTORY/minimap.alignment.sorted
# samtools index $ALIGNMENT_DIRECTORY/minimap.alignment.sorted.bam

###############################

#LAST
# Add in directory with humandb stuff (radich_j/ngs/Last-...)
# Align to human genome using LAST

module load LAST/926-foss-2016b
module load Python/2.7.15-foss-2016b

# Only needed once - already stored now in $LASTAL_HUMANDB
# lastdb -v -P4 -uNEAR -R01 humandb $REFERENCE
# lastdb -v -P4 hg19.lastdb $REFERENCE

# Working way to maf file, but no last-train
# lastdb -v -P4 -uNEAR -R01 humandb $REFERENCE/genome.fa
# lastal -v -P4 -Q1 -D100 $RADICH_HOME/ngs/Last-Genome-data/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz | last-split > $ALIGNMENT_DIRECTORY/last.alignment.maf

# Last-train version
last-train -Q1 -P4 -v $LASTAL_HUMANDB/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > last-parameters
lastal -Q1 -P4 -vp last-parameters $LASTAL_HUMANDB/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz | last-split > last-alignment-maf

# Convert maf file into sam file for samtools
maf-convert -f $REFERENCE/hg19.dict sam -r 'ID:$ID PL:nanopore SM:$SAMPLE' $ALIGNMENT_DIRECTORY/last-alignment-maf > $ALIGNMENT_DIRECTORY/last-alignment-sam

# Sort SAM file, convert to BAM and create index
module load SAMtools/1.6-foss-2016b
samtools view -bS $ALIGNMENT_DIRECTORY/last-alignment-sam | samtools sort - $ALIGNMENT_DIRECTORY/last-alignment-sorted
samtools index $ALIGNMENT_DIRECTORY/last-alignment-sorted-bam

###############################

# Call structural variants
mkdir $NANOSV_DIRECTORY
module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $NANOSV_DIRECTORY/NanoSV.vcf $ALIGNMENT_DIRECTORY/last-alignment-sorted-bam