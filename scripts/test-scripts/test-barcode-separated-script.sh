#!/bin/bash

# SPLIT EVERYTHING AFTER PORECHOP BY BARCODE
#   # Need to add names for separate barcodes AND vcf files should include sequence name?

##### Time for Nanoplot/LAST/minimap2/NanoSV 3 Gbp read: 235m25.799s

# Aligner
echo "Enter Aligner to be used (LAST or minimap2 or both)"
read ALIGNER
# Convert input to lower-case
ALIGNER=${ALIGNER,,}

#/path/to/radich home
export RADICH_HOME="/fh/fast/radich_j"

export WORKING_FOLDER="$RADICH_HOME/nanopore/Jobs/barcode-separated-testing"

#/path/to/Fastq/Demultiplexed (Porechop Folder)
export DEMULTIPLEXED_DIRECTORY="$RADICH_HOME/nanopore/Jobs/Wed_Aug_22_13:53:12_PDT_2018./porechop-demux"

#/path/to/NanoPlot
export NANOPLOT="$WORKING_FOLDER/nanoplot_Wed_Aug_22_13:53:12_PDT_2018."

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
module load Python/3.6.5-foss-2016b-fh3
###############################

# Nanoplot

###############################

echo -e "\n-----------Beginning Nanoplot processing-----------\n"

# Create NanoPlot figures
# NanoStats report included
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
for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do 
    nanoplot_function "$file"
done

echo -e "\n-----------Nanoplot completed-----------\n"

###############################

# Branch here LAST or minimap2

###############################

mkdir $ALIGNMENT_DIRECTORY
mkdir $NANOSV_DIRECTORY

#LAST -Align to human genome using LAST

if [ "$ALIGNER" = "last" ] || [ "$ALIGNER" = "both" ]; then

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

    #Run LAST, samtools, and NanoSV separated to prevent module interference. Each process is done in parallel accross all barcodes
    module load LAST/926-foss-2016b
    module load Python/2.7.15-foss-2016b
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do last_function "$file"; done

    module load samtools
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do last_samtools_function "$file"; done
    
    module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do last_nanosv_function "$file"; done

    echo -e "\n-----------LAST/NanoSV processing completed-----------\n"

fi

###############################

# MiniMap2 - Align to human genome using MiniMap2 

if [ "$ALIGNER" = "minimap2" ] || [ "$ALIGNER" = "both" ]; then

    mkdir $MINIMAP2_ALIGNMENT_DIRECTORY
    module load minimap2/2.10-foss-2016b

    mkdir $MINIMAP2_NANOSV_DIRECTORY

    minimap2_function () {
        local file=$1
        local file_name=$(basename "${file%.fastq.gz}")

        minimap2 -a $REFERENCE/genome.fa $file >> $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment.sam"

        # Sort SAM file, convert to BAM and create index
        module load samtools
        samtools view -bS $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment.sam" | samtools sort - $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment-sorted"
        samtools index $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment-sorted.bam"

        # Call structural variants
        module load NanoSV/1.1.2-foss-2016b-Python-3.6.4
        NanoSV -s samtools -c config.ini -b $RADICH_HOME/osala/nanopore/hg19-refFlat.bed -o $MINIMAP2_NANOSV_DIRECTORY/"${file_name}_minimap_NanoSV.vcf" $MINIMAP2_ALIGNMENT_DIRECTORY/"${file_name}_minimap-alignment-sorted.bam"
    }

    #Run minimap2 on each barcoded entry in parallel
    for file in $DEMULTIPLEXED_DIRECTORY/*.fastq.gz; do 
        minimap2_function "$file"
    done

fi

echo -e "\n-----------All processes completed-----------\n"

#Next steps fusions 

#polyphen
#siftmutation
#mutation Taster
#gerp
#dbsnp
#1000 genomes
#archer software
#qiagen fusion software