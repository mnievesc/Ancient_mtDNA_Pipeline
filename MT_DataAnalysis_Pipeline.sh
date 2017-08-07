#!/bin/bash

##################################################################################################################################################
# Analysis Pipeline for ancient mtDNA Capture Data
# Written by M. Nieves Colon - 9/16/2016 - Modified 8/7/2017 to fix bugs
# Usage: MT_DataAnalysis_Pipeline.sh
##################################################################################################################################################

# See README.md file for information on each step. Input R1 and R2 fastq.gz files must be in a directory named 1_RawFastq and named in the 
# following format: Sample1_R1_001.fastq.gz Sample1_R2_001.fastq.gz. Output file is one large text delimited file which can be opened and edited
# in excel. 

# Necessary dependencies
qmap=~/install/qualimap_v2.2/qualimap
ref=~/Maria/ref_seqs/rCRS.fasta 
contaMix=/home/user/install/contamMix/exec/estimate.R
mt311=~/Maria/ref_seqs/mt311.fa

echo " "
echo "****************************************************"
echo " ***   Welcome to the Ancient mtDNA Pipeline!     ***"
echo " ***                                              ***"
echo " ***   Script by M. Nieves-Colon - 9/16-20/16     ***"
echo " ***						***"	
echo " ***  Make sure fastq is in directory 1_RawFastq  ***"
echo "****************************************************"
echo " "
echo "Reference is: $ref"
echo " "

## 1. Set up sample information for output tab delimited file.
echo -e "Sample" > SampleList.txt   # echo -e interprets \t as tab

cd 1_Raw_Fastq
for a in *R1_001.fastq.gz;
 do 
  ID=$(basename $a)	
  NAME=$(echo $ID |cut -d "_" -f1)   
  echo "${NAME}" >> ../SampleList.txt		
done
cd ..
echo "***********************************"
echo "**** `grep -v "Sample" SampleList.txt | wc -l | cut -d " " -f1 `  Samples to be analyzed: ******"
cat SampleList.txt
echo "***********************************"
echo " "
echo " "



### 2A. Pre-Merging FastQC
echo "********************************************"
echo "**** Running fastQC on unmerged data ******"
mkdir -p 2_FastQC/{Pre-merging,Post-merging}   # creates multiple directories at once
fastqc -q 1_Raw_Fastq/*.fastq.gz		     # use quiet option to supress stdout messages	
mv 1_Raw_Fastq/*.html 2_FastQC/Pre-merging
mv 1_Raw_Fastq/*.zip 2_FastQC/Pre-merging
echo "********************************************"
echo " "
echo " "


### 3. SeqPrep for adapter removal and read merging. Save output to textfile as well as to stdout. Take values of interest and save as variables to
###    output to file.
echo "***********************************"
echo " ****** Running SeqPrep ******" 
mkdir 3_SeqPrep

cd 1_Raw_Fastq
echo -e "Sample \t "Premerge reads" \t "Postmerge reads" \t "% reads merged and kept" " > SeqPrepStats.txt

for a in *R1_001.fastq.gz;
 do 
  ID=$(basename $a)	
  NAME=$(echo $ID |cut -d "_" -f1) 
  echo "Running SeqPrep on sample "${NAME}"..."

  SeqPrep -f "${NAME}"_R1_001.fastq.gz -r "${NAME}"_R2_001.fastq.gz -o 11 -1 "${NAME}"_R1.trimmed.fastq.gz -2 "${NAME}"_R2.trimmed.fastq.gz -3 "${NAME}"_R1.discarded.fastq.gz -4 "${NAME}"_R2.discarded.fastq.gz -L 30 -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -s  "${NAME}".trimmed.merged.fastq.gz |& tee "${NAME}".out   # Need to use & sign to direct stderr (which has screen output) as well as stdout 

  allreads=`cat "${NAME}".out | awk "FNR==2" | awk '{print $3}'`    # Save number as variable
  mergedreads=`cat "${NAME}".out | awk "FNR==3" | awk '{print $3}'` 
  echo -e "$NAME \t $allreads \t $mergedreads \t `echo "scale=4;$mergedreads/$allreads" | bc` " >> SeqPrepStats.txt    # use scale to format output 
  rm  ${NAME}.out
done

mv *.discarded.fastq.gz ../3_SeqPrep
mv *.trimmed.* ../3_SeqPrep
mv SeqPrepStats.txt ../3_SeqPrep
cd ..

echo "***********************************"
echo " "



### 2B. Post-Merging FastQC
echo "**************************************************"
echo "**** Running fastQC on merged, trimmed data ******"
fastqc -q 3_SeqPrep/*.trimmed.merged.fastq.gz		     # use quiet option to supress stdout messages	
mv 3_SeqPrep/*.html 2_FastQC/Post-merging
mv 3_SeqPrep/*.zip 2_FastQC/Post-merging
echo "**************************************************"
echo " "
echo " "



### 4. Map reads and filter with samtools. See read me for details on what each command does
echo "****************************************"
echo "**** Mapping and filtering reads  ******"
mkdir 4_MappingFiltering   # change name to mapping and filtering
cd 4_MappingFiltering

echo -e "Sample \t Mapped Reads \t  % Mapped \t Q30 Mapped Reads \t Removed Duplicates \t Duplicate Rate (Duplicates per Q30 Mapped reads) \t \
Mapped Unique Reads \t % endogenous (Mapped Unique/Total Reads) \t Mapped Reads average length \t \
Cluster Factor (Total mapped reads / Unique reads)" > MappingFilteringStats.txt
cp ../3_SeqPrep/*.trimmed.merged.fastq.gz .

for a in *.fastq.gz;
 do 
  ID=$(basename $a)	
  NAME=$(echo $ID |cut -d "." -f1) 
  echo "Mapping data for ${NAME}"
  
  bwa aln -l 1000 -n 0.01 $ref ${NAME}.trimmed.merged.fastq.gz > ${NAME}.trimmed.merged.sai   # Map with seed disabled
  bwa samse $ref ${NAME}.trimmed.merged.sai ${NAME}.trimmed.merged.fastq.gz > ${NAME}.trimmed.merged.bwa.all.sam  # generate SAM alignment
  samtools view -bSh ${NAME}.trimmed.merged.bwa.all.sam > ${NAME}.trimmed.merged.bwa.all.bam  # Generate bam file
  samtools view -bh -F4 ${NAME}.trimmed.merged.bwa.all.bam > ${NAME}.trimmed.merged.mapped.bwa.bam  # Filter unmapped reads (F4 flag))
  samtools view -bh -q 30 ${NAME}.trimmed.merged.mapped.bwa.bam > ${NAME}.trimmed.merged.mapped.q30.bwa.bam  # Filter Q30 quality 
  samtools sort ${NAME}.trimmed.merged.mapped.q30.bwa.bam ${NAME}.trimmed.merged.mapped.q30.bwa.sort # Sort alignment by leftmost coordinate 
  samtools rmdup -s ${NAME}.trimmed.merged.mapped.q30.bwa.sort.bam ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.bam # Remove duplicates
  # Remove reads with multiple mappings	 
  samtools view -h ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.bam | grep -v 'XT:A:R'| grep -v 'XA:Z' |grep -v 'XT:A:M' | awk '{if($0~/X1:i:0/||$0~/^@/  )print $0}' | samtools view -bS - > ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam 

  # Print statistics per sample
  mapped=`samtools view -c ${NAME}.trimmed.merged.mapped.bwa.bam`
  q30=`samtools view -c ${NAME}.trimmed.merged.mapped.q30.bwa.sort.bam`
  rmdup=`samtools view -c ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.bam`
  uniq=`samtools view -c ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam`
  length=`samtools view ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam | awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'`  #fixed error in length estimation

  echo -e "$NAME \t $mapped \t `echo "calculate % mapped here"` \t $q30 \t $rmdup \t `echo "scale=4;($q30-$rmdup)/$rmdup" | bc` \t $uniq \t `echo "calculate endogenous here"` \t $length \t `echo "scale=4;$uniq/$q30" | bc` " >> MappingFilteringStats.txt
  echo " "
  echo " "
done

rm *.trimmed.merged.fastq.gz 
echo "***********************************"
echo " "
echo " "
cd ../



### 5. Generate deamination plots and rescale bam files with mapDamage (assumes reference is indexed)
echo "***********************************************".
echo "**** Estimate damage and rescale bamfile ******"
mkdir 5_mapDamage
cd 5_mapDamage
cp ../4_MappingFiltering/*uniq.bam .
echo -e "Sample \t Rescaled?" > mapDamage.RescaledList.txt

# Using parallel because it cuts processing time by at least half
find *.bam | parallel -j4 -k 'mapDamage -i {} -r ~/Maria/ref_seqs/rCRS.fasta -q --rescale'  # -quiet 

for a in *.uniq.bam;
 do 
  ID=$(basename $a)	
  NAME=$(echo $ID |cut -d "." -f1) 
  echo "Estimating Damage for ${NAME}"
  # Only erase uniq.bam file if the --rescaling exits with non zero exit status (i.e. rescaling works). Otherwise will get stats on uniq.bam file.
  [ -f results*/${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.rescaled.bam ] && echo -e "${NAME} \t Rescaled" >> mapDamage.RescaledList.txt || echo -e "${NAME} \t NotRescaled"  >> mapDamage.RescaledList.txt
  [ -f results*/${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.rescaled.bam ] && rm ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam || echo "No rescale performed"
done
mv */*rescaled.bam .
echo "***********************************************"
echo " "
echo " "
cd ..



## 6. Generate mapping stats report with Qualimap   
## 8-7-17: Fixed bug that would crash loop if a fastq file had no mapped reads (for example a library blank)
## Added check for zero exit status
echo "*****************************************"
echo "****  Generate mapping stats       ******"
mkdir 6_QualimapStats
cd 6_QualimapStats

echo -e "Sample \t Total mtDNA reads \t Mean Read Depth \t StdDev Read Depth \t % ref-seq covered >1x \t % ref-seq covered >2x " > QualimapStats.txt

cp ../5_mapDamage/*.bam .
ls *.bam   > bamlist
  
cat bamlist | while read line
 do
 sample=`echo $line | cut -d"." -f1`
 echo "Generating mapping stats for ${sample}"

$qmap bamqc -bam ${line} -outdir ${sample}.QualimapStats -outformat pdf  > /dev/null  # Do not print output to screen

 if [ $? -eq 0 ] 
  then
   echo " "
   cd ${sample}.QualimapStats
   nreads=`grep "number of reads" genome_results.txt |cut -d "=" -f 2`
   meanrd=`grep "mean coverageData" genome_results.txt| cut -d "=" -f2`
   sdrd=`grep "std coverageData" genome_results.txt| cut -d "=" -f2`
   cov1=`grep "reference with a coverageData >= 1X" genome_results.txt | awk '{print $4}'`
   cov2=`grep "reference with a coverageData >= 2X" genome_results.txt | awk '{print $4}'`
   echo -e "$sample \t $nreads \t $meanrd \t $sdrd \t $cov1 \t $cov2" >> ../QualimapStats.txt
   cd ../ 
   
  else
   echo " " && echo "*** Error: Sample has no mapped, unique reads **" && echo " "
   echo -e "${sample} \t 0 \t NA \t NA \t NA \t NA" >> QualimapStats.txt
  fi

done

rm bamlist *.bam
cd ..
echo "*****************************************"
echo " "



### 7. Generate pileup, call variants and produce vcf file   ## Must update and test in ubuntu computer
echo "***********************************"
echo "****  Generate pileup and vcf  ******"
mkdir 7_mpileup_VCF
cd 7_mpileup_VCF

echo -e "Variant sites with just 1x coverage \t Variant sites with >1x coverage" > temp.MpileupStats.txt 
echo -e "Sample" > samplelist
cp ../5_mapDamage/*.bam .
ls *.bam  > bamlist 

cat bamlist | while read line
 do
 NAME=`echo $line | cut -d"." -f1`

 echo "Generating pileup for ${NAME}"
 echo ${NAME} >> samplelist
 echo "$line	1" > ${NAME}.sample.txt
 
 # Generate vcf with just variant sites
 samtools mpileup -f $ref -C 50 -Bu ${line} | bcftools view -vcg -s ${NAME}.sample.txt - > ${NAME}.jv.vcf 
 
 # Generate vcf with variant sites >1x coverage
 samtools mpileup -f $ref -C 50 -Bu ${line} | bcftools view -vcg -s ${NAME}.sample.txt - | grep -v "DP=1;" > ${NAME}.jv_covfilt.vcf 
 
 # Print statistics per sample
 v1=`grep -v "^#" ${NAME}.jv.vcf | grep -c "DP=1;"`     # Variant sites covered by just 1 read
 v2=`grep -v "^#" ${NAME}.jv.vcf | grep -c -v "DP=1;"`  # Variant sites covered by just more than 1 read 
 echo -e "$v1 \t $v2" >> temp.MpileupStats.txt
done

paste samplelist temp.MpileupStats.txt > MpileupStats.txt
rm temp.MpileupStats.txt samplelist bamlist
rm *.bam
echo "***********************************"
echo " "
cd ..



###8. Generate consensus with mia and run contammix   
# ContamMix takes a really long time to run so we can comment this out if we want to run everything else separately.
echo "******************************************************************"
echo "****  Generate consensus and running contamination analyses ******"
mkdir 8_MIA_contamMix
cd 8_MIA_contamMix

echo -e "Reads Used \t MAP authentic \t 95% quantiles \t Pr reads matched other genome better than consensus (crude cont upper bound) \t error rate" > tempcmix
echo -e "Sample" > samplelist

cp ../5_mapDamage/*.bam .
ls *.bam  > bamlist 

cat bamlist | while read line
 do
 NAME=`echo $line | cut -d"." -f1`
 echo -e "$NAME" >> samplelist
 
 if [[ "$line" != *"rescaled"* ]]; then
    echo -e "No contamination assessment performed \t NA \t NA \t NA \t NA" >> tempcmix
    
 else
  echo "Creating MIA consensus for ${NAME}"
  bedtools bamtofastq -i ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.rescaled.bam -fq ${NAME}.q30.mapped.uniq.chrM.rescaled.fastq   # Trasform rescaled bam back to fastq using bedtools
  mia -r $ref -f ${NAME}.q30.mapped.uniq.chrM.rescaled.fastq -c -C -U -i -F -k 14 -m ${NAME}.q30.mapped.uniq.chrM.rescaled.fastq.maln  # create maln file
  ma -M ${NAME}.q30.mapped.uniq.chrM.rescaled.fastq.maln.? -f 5 -I ${NAME} > ${NAME}.q30.chrM.rescaled.mia.consensus.fasta   # create consensus

  echo "Aligning consensus with contaminating seqs for ${NAME}"
  cat ${NAME}.q30.chrM.rescaled.mia.consensus.fasta $mt311 > ${NAME}.mt311.fasta   # add sample to contaminant aln
  mafft --auto ${NAME}.mt311.fasta > ${NAME}.mt311.MAFFT.fasta  # Do not output stdout to screen  # realign
  bwa index ${NAME}.q30.chrM.rescaled.mia.consensus.fasta  # Do not output stdout to screen
  bwa aln -l 1000 -n 0.01 ${NAME}.q30.chrM.rescaled.mia.consensus.fasta ${NAME}.q30.mapped.uniq.chrM.rescaled.fastq > ${NAME}.remapped.q30.sai   # Do not output stdout to screen # realign to consensus
  bwa samse ${NAME}.q30.chrM.rescaled.mia.consensus.fasta ${NAME}.remapped.q30.sai ${NAME}.q30.mapped.uniq.chrM.rescaled.fastq  > ${NAME}.remapped.q30.sam  # Do not output stdout to screen
  samtools faidx ${NAME}.q30.chrM.rescaled.mia.consensus.fasta  # Do not output stdout to screen
  samtools view -bSh ${NAME}.remapped.q30.sam > ${NAME}.remapped.q30.bam  # new bam file of sample reads aligned to consensusD

  echo "Running contamMix for ${NAME}"
  $contaMix --samFn ${NAME}.remapped.q30.bam  --malnFn ${NAME}.mt311.MAFFT.fasta --consId ${NAME} --figure ${NAME}.contamMix_fig | tee ${NAME}.contamMixout.txt
 
  readsused=`grep "consist of" ${NAME}.contamMixout.txt | awk '{print $3}'`
  map=`grep "MAP authentic" ${NAME}.contamMixout.txt | cut -d":" -f2`
  quantiles=`awk 'NR==9 {print $0}' ${NAME}.contamMixout.txt`
  pr=`awk 'NR==4 {print $0}' ${NAME}.contamMixout.txt | cut -d"(" -f2 | cut -d")" -f1`	 
  err=`grep "error rate" ${NAME}.contamMixout.txt | awk '{print $9}' | sed 's/).//g'`
  echo -e "$readsused \t $map \t $quantiles \t $pr \t $err" >> tempcmix
 fi
done

paste samplelist tempcmix > ContamMixStats.txt

# Clean up a bit within folder
mkdir MIA_consensus
mv *.q30.chrM.rescaled.mia.consensus.fasta MIA_consensus
mkdir contamMix_output 
mv *.pdf contamMix_output 
mv *.contamMixout.txt  contamMix_output
rm samplelist tempcmix bamlist
mkdir intermediatefiles
mv * intermediatefiles
mv intermediatefiles/ContamMixStats.txt .
mv intermediatefiles/contamMix_output/ .
mv intermediatefiles/MIA_consensus/ .

echo "******************************************************************"
echo " "
cd ..



####9. Run PMDtools to look at different damage threshold
## 8-7-17: Fixed bug so that samples without any reads would not get analyzed.
echo "**************************"
echo "****  Run PMD tools ******"
mkdir 9_PMDtools
cd 9_PMDtools

cp ~/Maria/Analysis_scripts/pmdtools.py .
cp ~/Maria/Analysis_scripts/pmd_hist_v3.R .  # copy auxiliary scripts
cp ~/Maria/Analysis_scripts/plotPMD.R .
cp ../4_MappingFiltering/*.uniq.bam .  # Cannot used rescaled.bam for PMDtools

for a in *.uniq.bam;
do 
 ID=$(basename $a)	
 NAME=$(echo $ID |cut -d "." -f1) 

 #Check if sample has 0 mapped reads
 uniq=`samtools view -c ${NAME}.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam`
 
 if [ $uniq == 0 ]
 then 
  echo " "
  echo "No mapped reads, will not run PMDtools on sample "${NAME}"..."
  echo -e "${NAME} \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA \t NA" > ${NAME}.PMDscore.dist.txt 
  
 else
  echo " "
  echo "Running PMDtools on sample "${NAME}"..."
 
  samtools calmd -u "${NAME}".trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam $ref > "${NAME}".q30.md.bam   # Add MD tag
  samtools view "${NAME}".q30.md.bam | python pmdtools.py -p > "${NAME}".pmdscores.txt     #Produce distribution of pmd scores
  Rscript pmd_hist_v3.R "${NAME}"  # Call R script to produce histograms and temporary textfile

  #Run pmdtools at thresholds 0 and 3. Append information on reads excluded, passed and % data passed by threshold to the output from R. One tab delimited file.
  samtools view -h "${NAME}".q30.md.bam | python pmdtools.py --threshold 0 --header --stats 2> "${NAME}".pmd0.out | samtools view -Sb - > "${NAME}".q30.md.pmd0.bam 
  tot=`cat "${NAME}".pmd0.out | awk -v FS=":" '{print $2}' | awk "FNR==6"`
  ex0=`cat "${NAME}".pmd0.out | awk -v FS=":" '{print $2}' | awk "FNR==7"`
  pass0=`cat "${NAME}".pmd0.out | awk -v FS=":" '{print $2}' | awk "FNR==8"`

  samtools view -h "${NAME}".q30.md.bam | python pmdtools.py --threshold 3 --header --stats 2> "${NAME}".pmd3.out | samtools view -Sb - > "${NAME}".q30.md.pmd3.bam 
  ex3=`cat "${NAME}".pmd3.out | awk -v FS=":" '{print $2}' | awk "FNR==7"`
  pass3=`cat "${NAME}".pmd3.out | awk -v FS=":" '{print $2}' | awk "FNR==8"`

  echo -e "$tot \t 0 \t $ex0 \t $pass0 \t `echo "scale=2;$pass0/$tot" | bc` \t 3 \t $ex3 \t $pass3 \t `echo "scale=2;$pass3/$tot" | bc` " > tempnos
  paste "${NAME}".PMDscore.dist.temp.txt tempnos > "${NAME}".PMDscore.dist.txt   # paste R output and these results together in one file

  #Get deamination plots for data before subsetting, at thresholds 0 and threshold 3
  samtools view  "${NAME}".trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam | python pmdtools.py --deamination --range 30 > PMD_temp.txt
  Rscript plotPMD.R
  cp PMD_plot.pdf "${NAME}".q30.md.deamplot.pdf

  samtools view  "${NAME}".q30.md.pmd0.bam | python pmdtools.py --deamination --range 30 > PMD_temp.txt
  Rscript plotPMD.R
  cp PMD_plot.pdf "${NAME}".q30.md.pmd0_deamplot.pdf

  samtools view  "${NAME}".q30.md.pmd3.bam | python pmdtools.py --deamination --range 30 > PMD_temp.txt
  Rscript plotPMD.R
  mv PMD_plot.pdf "${NAME}".q30.md.pmd3_deamplot.pdf

  # Cleanup  and organize all files
  rm "${NAME}".PMDscore.dist.temp.txt tempnos
  mkdir pmdanalysis."${NAME}" 
  mv "${NAME}".* pmdanalysis."${NAME}"
  cp pmdanalysis."${NAME}"/"${NAME}".PMDscore.dist.txt  .
 
 fi
done

echo -e " Sample \t Min.PMD \t Max.PMD \t Mean.PMD \t Median.pmd \t Mode pmd \t Mapped.unique.reads \t Threshold \t Reads.excluded \t Reads.passed \t %data.kept \t Threshold \t Reads.excluded \t Reads.passed \t %data.kept" > header	
cat header *.PMDscore.dist.txt > PMDscore.distribution.Stats.txt
rm header
rm *.PMDscore.dist.txt
rm *uniq.bam
echo "**************************"
echo " "
cd ..
 


###10. Concatenate all textfiles into one large results file 
### 8-7-17: Fixed bug in paste command where tabbing was not correct.
cp */*Stats.txt .  # Copy all stats files to parent directory. Should be 7 in total.
cp 5_mapDamage/mapDamage.RescaledList.txt .

# Save columns we need
paste <(awk '{print $0}' SeqPrepStats.txt) <(cut -f 2-9 MappingFilteringStats.txt) <(awk '{print $2}' mapDamage.RescaledList.txt) <(cut -f 2-6 QualimapStats.txt ) <(cut -f 2-3  MpileupStats.txt )  <(cut -f 2-6 ContamMixStats.txt) <(cut -f 2-15 PMDscore.distribution.Stats.txt) > MT.DataAnalysis.Summary.txt
 
rm *Stats.txt
rm mapDamage.RescaledList.txt

echo " "
echo "****************************************************************************************"
echo " ***       You made it!    : )  Now go write a paper explaining what it all means!   ***"
echo "****************************************************************************************"




## References
# 1. http://askubuntu.com/questions/639877/tee-doesnt-get-whole-output-from-the-pipe
# 2. http://stackoverflow.com/questions/27388136/send-the-unix-output-to-a-csv-file
# 3.http://stackoverflow.com/questions/12745834/awk-extract-different-columns-from-many-different-files
# 4. http://stackoverflow.com/questions/1602035/how-to-print-third-column-to-last-column
# 5. http://www.cyberciti.biz/tips/find-out-if-file-exists-with-conditional-expressions.html
# 6. http://stackoverflow.com/questions/14174230/check-whether-files-in-a-file-list-exist-in-a-certain-directory
# 7. http://timmurphy.org/2013/05/13/string-contains-substring-in-bash/
# 8. http://stackoverflow.com/questions/617182/with-bash-scripting-how-can-i-suppress-all-output-from-a-command
# 9. http://www.thegeekstuff.com/2010/02/awk-conditional-statements
# 10. https://www.tutorialspoint.com/unix/if-else-statement.htm
# 11. http://bencane.com/2014/09/02/understanding-exit-codes-and-how-to-use-them-in-bash-scripts/
