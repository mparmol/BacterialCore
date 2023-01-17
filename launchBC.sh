#!/bin/bash

## A wrapper for launching BacterialCore in batch on a cluster

## Directories
### SUBMITDIR: where all files and scripts are; will be copied temporarily to the workdir
### WORKDIR: temporary location of input files and scripts ($SUBMITDIR contents) during the job
export SUBMITDIR="/home/silviatm/micro/bc"
export WORKDIR="/temporal/silviatm/"

## Input data (relative to $SUBMITDIR)
export INPUTFA="tomate_subset.fa"

## Simulation parameters
export PERC="1"
export CUTOFFS="0.05 0.1 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5"
export TAXO_P="0.9"
export CORES="16"

for cutoff in $CUTOFFS # fíjate que 0.01 es: solo considerar nodos con abundancia del 1%
do
  echo $cutoff
  export outputfolder="tomate_BC_"$cutoff

  echo -e \#\!/bin/bash > temp$cutoff
  echo -e \#SBATCH \-o slurm.\%N.\%j_$cutoff.out \# STDOUT >> temp$cutoff
  echo -e \#SBATCH \-e slurm.\%N.\%j_$cutoff.err \# STDERR >> temp$cutoff

  echo '

  module load qiime/1.9.1
  module load R/3.5.0

  cp '$SUBMITDIR'/'$INPUTFA' '$WORKDIR'
  cp '$SUBMITDIR'/BacterialCore.py '$WORKDIR'
  
  # create output folder in working directory
  mkdir -p '$WORKDIR'/'$outputfolder'
  
  python BacterialCore.py -f '$WORKDIR'/'$INPUTFA' -o '$WORKDIR'/'$outputfolder'/ -p 3 -initial_level 0.99 -min_core_percentage '$PERC' -cutoff '$cutoff' -taxo_p '$TAXO_P' -tree_level 99 -tree_type_analysis 3 -threads '$CORES'
  
  # save and remove
  cp -prf '$WORKDIR'/'$outputfolder' '$SUBMITDIR'
  rm -rf '$WORKDIR'/'$outputfolder'
  rm '$WORKDIR'/'$INPUTFA'
  rm '$WORKDIR'/BacterialCore.py' >> temp$cutoff

  sbatch -A microbioma_serv -p biobis temp$cutoff
  
  # remove the temporal script
  rm temp$cutoff
done
