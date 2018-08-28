#!/bin/bash

# MinIon pipeline for structural variant calls v7
# Lots of credit to David Coffey

###TBA###
# ADD IN INPUT FOR NAME + DATE OF TIME -v7
# Add both option for aligner -v7
# SPLIT EVERYTHING AFTER PORECHOP BY BARCODE
#   # Need to add names for separate barcodes AND vcf files should include sequence name?
# Add date for nanoplot folder when scp'ing over -v7
# Nanopolish
# Metapore
# Charlie Hill
#########

# Start by reading inputs from the user - [SEQUENCE_NAME] [/path/to/data] [LAST, minimap2, both]

# Sequence name
echo "Enter sequence name"
read SEQUENCE_NAME

# /path/to/data
echo "Enter /path/to/fast5"
read FAST5

# Aligner
echo "Enter Aligner to be used (LAST or minimap2 or both)"
read ALIGNER
# Convert input to lower-case
ALIGNER=${ALIGNER,,}

#/path/to/radich home
export RADICH_HOME="/fh/fast/radich_j"

#Working folder for the session is created with the date and timestamp. All processes will be run in this folder
DATE=$(date)
NAME_DATE_FORMATTED="${SEQUENCE_NAME}_${DATE// /_}"
export WORKING_FOLDER="$RADICH_HOME/nanopore/Jobs/$NAME_DATE_FORMATTED"
mkdir $WORKING_FOLDER

export FLOWCELL="FLO-MIN106"
export KIT="SQK-LSK108"

#/path/to/Fastq/Original (Albacore Folder)
export FASTQ="$WORKING_FOLDER/albacore-fastq"

#/path/to/Fastq/Demultiplexed (Porechop Folder)
export DEMULTIPLEXED_DIRECTORY="$WORKING_FOLDER/porechop-demux"

#/path/to/NanoPlot
export NANOPLOT="$WORKING_FOLDER/nanoplot_$NAME_DATE_FORMATTED"

#/path/to/ReferenceGenomes
export REFERENCE="/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta"

#/path/to/LastHumandb
export LASTAL_HUMANDB="$RADICH_HOME/nanopore/Last-Genome-Data" 

#/path/to/Aligned
export ALIGNMENT_DIRECTORY="$WORKING_FOLDER/aligned"

#/path/to/LastAligned
export LAST_ALIGNMENT_DIRECTORY="$ALIGNMENT_DIRECTORY/last-aligned"

#/path/to/Minimap2Aligned
export MINIMAP2_ALIGNMENT_DIRECTORY="$ALIGNMENT_DIRECTORY/minimap2-aligned"

#/path/to/NanoSV
export NANOSV_DIRECTORY="$WORKING_FOLDER/nanosv-vcf"

#/path/to/LastNanoSV
export LAST_NANOSV_DIRECTORY="$NANOSV_DIRECTORY/last-vcf"

#/path/to/Minimap2NanoSV
export MINIMAP2_NANOSV_DIRECTORY="$NANOSV_DIRECTORY/minimap2-vcf"

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
--verbose \
--threads 12 \
--plots hex dot

#Branch here MiniMap2 or LAST

###############################

mkdir $ALIGNMENT_DIRECTORY
mkdir $NANOSV_DIRECTORY

#LAST -Align to human genome using LAST

if [ "$ALIGNER" = "last" ]; then

    mkdir $LAST_ALIGNMENT_DIRECTORY
    module load LAST/926-foss-2016b
    module load Python/2.7.15-foss-2016b

    # Only needed once - already stored now in $LASTAL_HUMANDB
    # lastdb -v -P4 -uNEAR -R01 humandb $REFERENCE
    # lastdb -v -P4 hg19.lastdb $REFERENCE

    # Working way to maf file, but no last-train
    # lastdb -v -P4 -uNEAR -R01 humandb $REFERENCE/genome.fa
    # lastal -v -P4 -Q1 -D100 $RADICH_HOME/nanopore/Last-Genome-data/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz | last-split > $LAST_ALIGNMENT_DIRECTORY/last.alignment.maf

    # Last-train version
    last-train -Q1 -P4 -v $LASTAL_HUMANDB/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > $LAST_ALIGNMENT_DIRECTORY/last-parameters
    lastal -Q1 -P4 -vp $LAST_ALIGNMENT_DIRECTORY/last-parameters $LASTAL_HUMANDB/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz | last-split > $LAST_ALIGNMENT_DIRECTORY/last-alignment-maf

    # Convert maf file into sam file for samtools
    maf-convert -f $REFERENCE/genome.dict sam $LAST_ALIGNMENT_DIRECTORY/last-alignment-maf > $LAST_ALIGNMENT_DIRECTORY/last-alignment-sam

    # Sort SAM file, convert to BAM and create index
    module load samtools
    samtools view -bS $LAST_ALIGNMENT_DIRECTORY/last-alignment-sam | samtools sort - $LAST_ALIGNMENT_DIRECTORY/last-alignment-sorted
    samtools index $LAST_ALIGNMENT_DIRECTORY/last-alignment-sorted-bam

    # Call structural variants
    mkdir $LAST_NANOSV_DIRECTORY
    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $LAST_NANOSV_DIRECTORY/NanoSV.vcf $LAST_ALIGNMENT_DIRECTORY/last-alignment-sorted-bam

###############################

# MiniMap2 - Align to human genome using MiniMap2 

elif [ "$ALIGNER" = "minimap2" ]; then

    mkdir $MINIMAP2_ALIGNMENT_DIRECTORY
    module load minimap2/2.10-foss-2016b
    minimap2 -a $REFERENCE/genome.fa $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sam

    # Sort SAM file, convert to BAM and create index
    module load samtools
    samtools view -bS $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sam | samtools sort - $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sorted
    samtools index $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sorted-bam

    # Call structural variants
    mkdir $MINIMAP2_NANOSV_DIRECTORY
    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $MINIMAP2_NANOSV_DIRECTORY/NanoSV.vcf $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sorted-bam

# Both LAST and minimap2

else

    #LAST
    mkdir $LAST_ALIGNMENT_DIRECTORY
    module load LAST/926-foss-2016b
    module load Python/2.7.15-foss-2016b

    # Last-train version
    last-train -Q1 -P4 -v $LASTAL_HUMANDB/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > last-parameters
    lastal -Q1 -P4 -vp last-parameters $LASTAL_HUMANDB/humandb $DEMULTIPLEXED_DIRECTORY/none.fastq.gz | last-split > last-alignment-maf

    # Convert maf file into sam file for samtools
    maf-convert -f $REFERENCE/genome.dict sam $LAST_ALIGNMENT_DIRECTORY/last-alignment-maf > $LAST_ALIGNMENT_DIRECTORY/last-alignment-sam

    # Sort SAM file, convert to BAM and create index
    module load samtools
    samtools view -bS $LAST_ALIGNMENT_DIRECTORY/last-alignment-sam | samtools sort - $LAST_ALIGNMENT_DIRECTORY/last-alignment-sorted
    samtools index $LAST_ALIGNMENT_DIRECTORY/last-alignment-sorted-bam

    # Call structural variants
    mkdir $LAST_NANOSV_DIRECTORY
    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $LAST_NANOSV_DIRECTORY/NanoSV.vcf $LAST_ALIGNMENT_DIRECTORY/last-alignment-sorted-bam

    #######

    #minimap2
    mkdir $MINIMAP2_ALIGNMENT_DIRECTORY
    module load minimap2/2.10-foss-2016b
    minimap2 -a $REFERENCE/genome.fa $DEMULTIPLEXED_DIRECTORY/none.fastq.gz > $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sam

    # Sort SAM file, convert to BAM and create index
    module load samtools
    samtools view -bS $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sam | samtools sort - $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sorted
    samtools index $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sorted-bam

    # Call structural variants
    mkdir $MINIMAP2_NANOSV_DIRECTORY
    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $MINIMAP2_NANOSV_DIRECTORY/NanoSV.vcf $MINIMAP2_ALIGNMENT_DIRECTORY/minimap-alignment-sorted-bam

fi


#Next steps fusions 

#polyphen
#siftmutation
#mutation Taster
#gerp
#dbsnp
#1000 genomes
#archer software
#qiagen fusion software