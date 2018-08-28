#!/bin/bash

# MinIon pipeline for structural variant calls parallelized v2.0

###TBA###
### Nanoplot, LAST, and minimap2 all parallel
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

###############################

# Albacore

###############################

echo -e "\n-----------Beginning Albacore processing-----------\n"

# Basecall
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
#cat $FASTQ/workspace/pass/barcode*/*.fastq $FASTQ/workspace/pass/unclassified/*.fastq >> $FASTQ/workspace/pass/merged.fastq
cat $FASTQ/workspace/pass/*.fastq >> $FASTQ/workspace/pass/merged.fastq

echo -e "\n-----------Albacore processing completed-----------\n"

###############################

# Porechop

###############################

echo -e "\n-----------Beginning Porechop processing-----------\n"

# Demultiplex
mkdir $DEMULTIPLEXED_DIRECTORY
porechop \
--input $FASTQ/workspace/pass/merged.fastq \
--barcode_dir $DEMULTIPLEXED_DIRECTORY \
--format fastq \
--threads 12 \
-v 2

echo -e "\n-----------Porechop processing completed-----------\n"

###############################

# Nanoplot

###############################

echo -e "\n-----------Beginning Nanoplot processing-----------\n"

Create NanoPlot figures
NanoStats report included
mkdir $NANOPLOT

#Nanoplot function to run through barcodes
nanoplot_function () {
    local file=$1
    #Strip file name down to barcode name
    local file_name=$(basename "${file%.fastq.gz}")

    mkdir "$NANOPLOT/${file_name}"
    NanoPlot \
    --fastq $file \
    --outdir "$NANOPLOT/${file_name}" \
    --verbose \
    --threads 4 \
    --plots hex dot
}

#For
( for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do 
    nanoplot_function "$file" &
done 
wait
echo -e "\n-----------Nanoplot completed-----------\n" ) &



###############################

# Branch here LAST or minimap2

###############################

mkdir $ALIGNMENT_DIRECTORY
mkdir $NANOSV_DIRECTORY

#LAST -Align to human genome using LAST

( if [ "$ALIGNER" = "last" ] || [ "$ALIGNER" = "both" ]; then

    echo -e "\n-----------Beginning LAST/NanoSV processing-----------\n"

    mkdir $LAST_ALIGNMENT_DIRECTORY
    mkdir $LAST_NANOSV_DIRECTORY

    #LAST function to run in background within for loop
    last_function () {
        local file=$1
        #Strip file name down to barcode name
        local file_name=$(basename "${file%.fastq.gz}")

        # Last-train version
        last-train -Q1 -P4 -v $LASTAL_HUMANDB/humandb $file >> $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-parameters"
        lastal -Q1 -P4 -vp $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-parameters" $LASTAL_HUMANDB/humandb $file | last-split > $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment.maf"

        # Convert maf file into sam file for samtools
        maf-convert -f $REFERENCE/genome.dict sam $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment.maf" >> $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment.sam"
    }

    # Sort SAM file, convert to BAM and create index
    last_samtools_function () {
        local file=$1
        local file_name=$(basename "${file%.fastq.gz}")
        samtools view -bS $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment.sam" | samtools sort - $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment-sorted"
        samtools index $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment-sorted.bam"
    }

    # Call structural variants
    last_nanosv_function () {
        local file=$1
        local file_name=$(basename "${file%.fastq.gz}")
        NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $LAST_NANOSV_DIRECTORY/"${file_name}_last_NanoSV.vcf" $LAST_ALIGNMENT_DIRECTORY/"${file_name}_last-alignment-sorted.bam"
    }

    #Run LAST in it's own subshell to prevent module interference
    module load LAST/926-foss-2016b
    module load Python/2.7.15-foss-2016b
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do last_function "$file" & done

    module load samtools
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do last_samtools_function "$file" & done
    wait
    
    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do last_nanosv_function "$file" & done
    wait

    echo -e "\n-----------LAST/NanoSV processing completed-----------\n"

fi ) &

###############################

# MiniMap2 - Align to human genome using MiniMap2 


( if [ "$ALIGNER" = "minimap2" ] || [ "$ALIGNER" = "both" ]; then

    echo -e "\n-----------Beginning minimap2/NanoSV processing-----------\n"

    mkdir $MINIMAP2_ALIGNMENT_DIRECTORY
    mkdir $MINIMAP2_NANOSV_DIRECTORY

    minimap2_function () {
        local file=$1
        local file_name=$(basename "${file%.fastq.gz}")
        minimap2 -a $REFERENCE/genome.fa $file >> $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment.sam"
    }

    minimap2_samtools_function () {
        local file=$1
        local file_name=$(basename "${file%.fastq.gz}")
        samtools view -bS $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment.sam" | samtools sort - $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment-sorted"
        samtools index $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment-sorted.bam"
    }

    minimap2_nanosv_function () {
        local file=$1
        local file_name=$(basename "${file%.fastq.gz}")
        NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $MINIMAP2_NANOSV_DIRECTORY/"${file_name}_minimap_NanoSV.vcf" $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment-sorted.bam"
    }

    #Run minimap2, samtools, and NanoSV separated to prevent module interference. Each process is done in parallel accross all barcodes
    module load minimap2/2.10-foss-2016b
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do minimap2_function "$file" & done
    wait

    module load samtools
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do minimap2_samtools_function "$file" & done
    wait

    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do minimap2_nanosv_function "$file" & done
    wait

    echo -e "\n-----------minimap2/NanoSV processing completed-----------\n"

fi ) &
wait
echo "All processes finished"

#Next steps fusions 

#polyphen
#siftmutation
#mutation Taster
#gerp
#dbsnp
#1000 genomes
#archer software
#qiagen fusion software

# Initial framework built off David Coffey's pipeline