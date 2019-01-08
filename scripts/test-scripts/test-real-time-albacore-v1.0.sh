#!/bin/bash

# Testing albacore runs on sequence uploaded in increments

# /path/to/data
export FAST5="/fh/fast/radich_j/nanopore/Data/test3/fast5"

#/path/to/radich home
export RADICH_HOME="/fh/fast/radich_j"

export WORKING_FOLDER="$RADICH_HOME/nanopore/Jobs/real-time-testing"

export FLOWCELL="FLO-MIN106"
export KIT="SQK-LSK108"

#/path/to/Fastq/Original (Albacore Folder)
export FASTQ="$WORKING_FOLDER/albacore-fastq"

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

echo -e "\n-----------Albacore processing completed-----------\n"