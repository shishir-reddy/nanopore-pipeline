#!/bin/bash

#/path/to/radich home
export RADICH_HOME="/fh/fast/radich_j/ngs"

#Working folder for the session is created with the date and timestamp. All processes will be run in this folder
DATE=$(date)
DATE_FORMATTED="${DATE// /_}".
export WORKING_FOLDER="$RADICH_HOME/Jobs/$DATE_FORMATTED"
mkdir $WORKING_FOLDER

# MinIon pipeline for structural variant calls v1.0
# Lots of credit to David Coffey

# Define variables later for input fields

#/path/to/fast5/pass (Original Data)
export FAST5="$RADICH_HOME/Data/fast5"

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

#/path/to/Aligned
export ALIGNMENT_DIRECTORY="$WORKING_FOLDER/last-aligned"

#/path/to/NanoSV
export NANOSV_DIRECTORY="$WORKING_FOLDER/nanosv-vcf"

export SAMPLE="BL"
export ID="033018"

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

----

#MiniMap2
# Align to human genome using MiniMap2 
mkdir $ALIGNMENT_DIRECTORY
module load minimap2/2.10-foss-2016b
minimap2 -a $REFERENCE/hg19.mmi $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > $ALIGNMENT_DIRECTORY/minimap.alignment.sam

# Sort SAM file, convert to BAM and create index
samtools view -bS $ALIGNMENT_DIRECTORY/minimap.alignment.sam | samtools sort - $ALIGNMENT_DIRECTORY/minimap.alignment.sorted
samtools index $ALIGNMENT_DIRECTORY/minimap.alignment.sorted.bam

----

#LAST
# Align to human genome using LAST
cd $ALIGNMENT_DIRECTORY

module load LAST/926-foss-2016b
module load Python/2.7.15-foss-2016b

lastdb -uNEAR -R01 humandb $REFERENCE -v
lastdb hg19.lastdb $REFERENCE -v

last-train -Q1 $REFERENCE/hg19.lastdb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > $ALIGNMENT_DIRECTORY/last.parameters
lastal -Q1 -p $ALIGNMENT_DIRECTORY/last.parameters $REFERENCE/hg19.lastdb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz | last-split > $ALIGNMENT_DIRECTORY/last.alignment.maf
maf-convert -f $REFERENCE/hg19.dict sam -r 'ID:$ID PL:nanopore SM:$SAMPLE' $ALIGNMENT_DIRECTORY/last.alignment.maf > $ALIGNMENT_DIRECTORY/last.alignment.sam

# Sort SAM file, convert to BAM and create index
samtools view -bS $ALIGNMENT_DIRECTORY/last.alignment.sam | samtools sort - $ALIGNMENT_DIRECTORY/last.alignment.sorted
samtools index $ALIGNMENT_DIRECTORY/last.alignment.sorted.bam

----

# Call structural variants
mkdir $NANOSV_DIRECTORY
module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
NanoSV $ALIGNMENT_DIRECTORY/last.alignment.sorted.bam -c $NANOSV_DIRECTORY/config.ini > $NANOSV_DIRECTORY/NanoSV.vcf
