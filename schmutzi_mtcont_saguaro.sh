#!/bin/bash

#SBATCH -p cluster
#SBATCH -n 16
#SBATCH -t 0-72:00
##SBATCH -A mnievesc
#SBATCH -o /home/mnievesc/PRaDNA_MT_forschmutzi/schmutzi.%j.out
#SBATCH -e /home/mnievesc/PRaDNA_MT_forschmutzi/schmutzi.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mnievesc@asu.edu

## Script to run schmutzi contamination estimation program from mitochondrial reads.
## Written 10/16/17 by M. Nieves Colon for use on Saguaro cluster.
## Usage: schmutzi_mtcont_saguaro.sh

## Reference sequence and contaminant database
ref=/home/mnievesc/ref-seqs/rCRS.fasta
conts=/home/mnievesc/software/schmutzi/alleleFreqMT/197/freqs/

## Necessary modules for Saguaro
SAMTOOLS=/home/mnievesc/software/samtools-0.1.19/samtools
module load schmutzi/2017-10-13

## Set working directory
WD=/home/mnievesc/PRaDNA_MT_forschmutzi/

## Run script in directory containing all uniq.bamfiles. Must not run on rescaled files because
## Use for loop to run on all bamfiles at once.
cd $WD
for a in uniq.bam/*uniq.bam;
 do
  ID=$(basename $a)
  NAME=$(echo $ID |cut -d "." -f1)
  echo "############# STARTING ANALYSIS ################"
  echo "Running schmutzi pipeline on sample "${NAME}"..."

  ## Add MD tag (BAMfile must be sorted)
  $SAMTOOLS calmd -b uniq.bam/${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam $ref > ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam

  ## Index bamfiles
  $SAMTOOLS index ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam

  ## Run schmutzi deamination contamination estimate (used default length deam).
  ## contDeam.pl --lengthDeam [length] --library [library type] --out [output prefix] [mt reference] [input bam file]
  ## Output cont.est file is in format: estimate estimate_low estimate_high.
  contDeam.pl --library double --out ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam --ref $ref ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam

  ## Run schmutzi without predicting contaminant with reference dataset
  ## schmutzi.pl --notusepredC --uselength --ref [mt reference] --out [out prefix]_npred [output prefix] [path to schmutzi]/eurasian/freqs/ [input bam file]
  ## NOTE: HAS A BUG WHERE IT WANTS THE FULL BAM NAME AS OUTPUT PREFIX IN PREVIOUS STEP
  schmutzi.pl --notusepredC --uselength --ref $ref $conts --out ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam

  ## Depending on results you can now run schmutzi again with contaminant prediction for those samples with high contamination estimates.

  ## Organize output directory
  rm *bai  # delete index file
  mkdir ${NAME}.schmutzi
  mv ${NAME}* ${NAME}.schmutzi
  cd ${NAME}.schmutzi

  ## Batch rename schmutzi files and save in their own directory
  for i in ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam_*; do mv $i `echo $i | sed 's/trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam_//'`; done

  ## Batch rename contdeam files and save in their own directory
  for i in ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam.*; do mv $i `echo $i | sed 's/trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.md.bam.//'`; done

  ## Move and compress files, only leave the final contamination estimate in directory.
  mkdir schmutzi-iters.output
  mv ${NAME}.?_* schmutzi-iters.output
  tar -czvf schmutzi-iters.output.tar.gz schmutzi-iters.output
  rm -r schmutzi-iters.output

  mkdir contDeam.output
  mv ${NAME}.* contDeam.output
  mv contDeam.output/*final* . && mv contDeam.output/*bam .

  cd $WD

  echo " #### DONE WITH ESTIMATION ####"
done
