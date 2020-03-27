**********************************************************************************************************************************
# Analysis Pipeline for ancient mtDNA Capture Data
### Written by M. Nieves Colon - 9/20/2016 - mnievesc@asu.edu / maria.nieves@cinvestav.mx
### Updated 9/30/2018
### Usage: ./MT_DataAnalysis_Pipeline.sh
**********************************************************************************************************************************

## READ ME

This pipeline processes Illumina reads from ancient mtDNA-captured samples in raw fastq format and performs QC, adapter trimming, genome mapping and filtering, 
calls variants, estimates contamination and damage. It is designed to work with multiple samples at once and output one final report in tab delimited txt format. 
Below an explanation of conditions used in each step and other important notes. Pipeline is designed to run within one parent directory already containing raw
fastq files (see section 1). Tab delimited stats files are created and stored in each sub-directory as well as in the parent directory.

**Update 9/30/2018:** Added script `schmutzi_mtcont_saguaro.sh` to run schmutzi contamination estimation for mitochondrial reads. This is not implemented in pipeline yet but can be run separately. Schmutzi is availabe on GitHub: https://github.com/grenaud/schmutzi.

**Update 3/27/2020:** Added notes for adapting pipeline to run on macOS. These notes were complied by Tre Blohm (tre.blohm@umconnect.umt.edu) PhD candidate at University of Montana between 2019 and 2020.

-----

###  Necessary dependencies and programs (versions noted)
1. Qualimap v2.2
2. SeqPrep 
3. FastQC v0.11.2
4. BWA v0.7.5
5. samtools v0.1.19
6. bcftools v0.1.19
7. mapDamage v2.0.2-12
8. pmdtools v0.50
9. mapping iterative assembler v1.0 (includes ma)
10. mafft v7.221
11. bedtools v2.17.0
12. gnu parallel v20130922
13. contamMix v1.0-10
14. R v3.3.1
15. Python 2.7.6

### Important: 
Reference needs to be specified at top of script, same for path to contamMix and Qualimap unless these are stored in PATH. 
In my case the reference is always rCRS.fasta *(gi|251831106|ref|NC_012920.1| Homo sapiens mitochondrion, complete genome)*. 
Pipeline assumes reference has been previously indexed using `samtools faidx` and `bwa index`.

In addition to the above listed programs there are also two helper scripts that need to be used with this pipeline: `pmd_hist_v3.R` (provided here), `plotPMD.R` (available at pmdtools website) and one reference file of contaminants which is necessary for contamMix: `mt311.fa`. This file can be obtained with the MIA program (see below).

----

### 1. Input files
Input files are raw, gzipped fastq files in format: `Sample1_R1_001.fastq.gz Sample1_R2_001.fastq.gz`. Rename accordingly. Files should be stored in a directory called
**1_Raw_Fastq**. The pipeline will make a list of samples to analyze by going into this directory and listing the basenames (separated by "_") for each sample. 
These should be the names you want the samples to have for the rest of the analysis. 

The sequencing center I use often names the files something like this: *Sample_S39_L001_R1_001.fastq.gz*
So before analysis I use the following commands to copy and rename the files. I always keep a copy of the raw data untouched and only perform analysis on the renamed files.
To get data in format for pipeline from raw basespace folders: 
```
cp */*/*.fastq.gz 1_Raw_Fastq/
rename 's/_S[0-9]+_L001_/_/' *.fastq.gz
```

----


### 2. FastQC 
FastQC is run using read 1 (R1) before merging and adapter trimming (Pre-merging) and after SeqPrep processing (Post-merging). See http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
for explanation of how to interpret FastQC reports. **Output from fastQC are .html and .zip files**.

----


### 3. SeqPrep Adapter trimming and merging
Step 3 takes fastq.gz as input for SeqPrep in order to trim off adapters and merge forward and reverse reads. Program also removes residual adapter sequences from 
NGS reads.  The program is run using default parameters. This can be changed depending on user requirements.  The -f and -r options specify the input R1 and R2 files
respectively. The -1 and -2 options specify the output R1 and R2 files to be written after trimming. The -s option specifies the final output file - the 
trimmed.merged.fastq.gz file. This file should be used for further analyses. The -L option ensures that all trimmed and merged sequences having a length of less
than 30 are discarded. The -A and -B options specify the forward and reverse adapters respectively.  The -o option specifies the amount of overlap between reads. 
For ancient DNA we modify this to 11 because our reads are short but most modern datasets use default value of 15. I use tee to save the output from screen to a textfile.
As mentioned above, if fastq files have extra numbers in filename the following sed one liner can serve to change basename and eliminate these (adjust number of
brackets according to how many digits need to be eliminated): `sed "s/_S[0-9][0-9]//" Sample1_S92_R1_001.fastq.gz > Sample1_R1_001.fastq.gz`
**Final output at this stage is a fastq file: `Sample.trimmed.merged.fastq.gz `**

----


### 4. Map reads to reference genome and filter for quality
Step 4 uses bwa to map reads to the reference and samtools to filter the reads. The following steps take place:
 1. bwa aln bwa aln  finds the SA coordinates of the input reads. Paramaters used include `-l 1000` which disables seed for usage with aDNA reads,
    `-n 0.01` which increases edit distance (makes alignment double as slow, remove if needed). 
 2. bwa samse generates alignments in the SAM format given single-end reads (ours are paired-ended, but already merged using SeqPrep)
 3. samtools view -bSh displays the previous output as a BAM file (b), input is SAM (S),and include the header (h). I also use samtools view -c
    to count alignments and save the total number as a variable for later printing to the stats file.
 4. samtols view -bh -F4 is for filtering out unmappped and low quality reads. -q displays the previous output as a BAM file (b), and including the header
    (h), but skipping alignments with MAPQ smaller than than 30 (-q 30), and alignments containing 4 flag ( 0x4 segment unmapped). 
 5. samtools sort sorts alignments by leftmost coordinates.
 6. samtools rmdup removes potential PCR duplicates: if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality.
    -s Removes duplicate for single-end reads. By default, the command works for paired-end reads only.
 7. Remove reads with multiple mappings using samtools and grep. In this step I keep single mapped reads from the rmdup file. grep -v means "invert the match", 
    it returns all non matching lines (all the reads that are not multiple mappings, so all the unique ones). The `XT:A:U` flag in the SAM file denotes unique 
    read, and `XT:A:R` and `XA:Z` denote multiple mappings for that read. When working with paired end reads, it may also be useful to consider filtering out 
    reads with the flag `XT:A:M` (one-mate recovered) which means that one of the pairs is uniquely mapped and the other isn't. Use awk to scan input file for
    lines that match the pattern: `{if($0~/X1:i:0/||$0~/^@/)print $0}` X1= Number of suboptimal hits found by BWA, if X1=0 then keep that read in file 
See the SAM format documentation for more information on each of these steps: https://samtools.github.io/hts-specs/SAMv1.pdf. After performing mapping and filtering
I save counts at different stages of filtering using variables (variable = `unix command that outputs desired value`) and I use echo -e to print these to the 
stats file for this step. Using `echo -e` we ensure that the shell understands we want elements in the file to be tab delimtied `\t`
**Final output at this stage is a BAM file: `Sample.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam`**
NOTE:Shoutout Drs. Rosa Fregel and Maria C. Avila Arcos who helped me in developing this part of the script. 

----


### 5. Estimate deamination patterns and rescale bam files
Step 5 uses mapDamage to estimate deamination patterns for ancient DNA samples. I also use the `--rescale` option to rescale mapping qualities while taking 
damage into account. This is important for later genotype calls. Since rescaling takes a considerable amount of time and computational resources, at this step 
I used gnu parallel to run 4 samples simultaneously. Further, since some samples will not have a lot of reads (ex. blanks or samples with low endogenous content) 
I created an exit status test. This ensures that samples with too few reads for rescaling can still be processed in the next step by copying the unrescaled .bam
file into the next process. For more information on mapDamage see: https://ginolhac.github.io/mapDamage/
**Final output at this stage is a  rescaled BAM file: `Sample.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.rescaled.bam`**
NOTE: Here there is an absolute path to the reference genome file, this must be changed if using another ubuntu machine.

----


### 6. Generate summary statistics per sample with Qualimap
Step 6 uses Qualimap to generate overall statistics per BAM file. These include overall read depth, sequence coverage, amount of reads, etc.. Only a subset of these
is reported in the final text file at the end of the pipeline. The rest of the information is stored in .pdf format in individual directories `Sample.QualimapStats`.
For more infromation on the summary data calculated by Qualimap see: http://qualimap.bioinfo.cipf.es/doc_html/index.html. Other programs can be used for generating
these summary statistics such as Geneious. In this step I used an if statement with another exit status test to  make sure that Qualimap is only run for samples that 
actually have mapped, unique reads. I added this step because if a sample has no reads (because its a library or extraction blanks or simply has no endogenous aDNA
mapping to the reference) Qualimap will error out and this will ruin the final output.

----


### 7. Call variants and produce VCF
Step 7 uses samtools and bcftools to generate genotype calls. By default this process is supposed to be for diploid organisms so to use it with mtDNA we need to 
specify that ploidy = 1 using a textfile containing the filename and ploidy:`echo "sample.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.rescaled.bam	1" > sample.txt`.
Use samtools mpileup to create vcf files from bam files and bcftools view for sub setting and filtering VCF/BCF files. 
We then run samtools with the following options used: -f input indexed reference, -u compute genotype likelihoods and output them in the binary call format (BCF)
uncompressed - C 50 Coefficient for downgrading mapping quality for reads containing excessive mismatches, 50 is recommended value for BWA alignments. 
Optional parameters: -d INT sets max number of read depth to call a SNP, you can use this if you have really high coverage (8000x+) limit memory usage 
(if you have aDNA you likely do not need this). It may also be necessary to disable BAQ computation (-B)  with aDNA data if some variants are skipped by samtools.
But careful, you still need to check bamfile by hand to confirm if some variant is missing. For more: http://samtools.sourceforge.net/samtools.shtml. 
BCFTOOLS options used: -c min/max count for non-reference, -g genotype, -s samples list, comma separated list of samples to include.
see above. Optional: -v output only variant sites. Also we could output as bcf vs vcf for compression (-b). VCF files produced can be used as input for
mtDNA haplotyping in Haplogrep: http://haplogrep.uibk.ac.at/. 
**Final output at this stage are two VCF file, one is filtered for >1x read depth and the other is not: `Sample.jv.vcf & Sample.jv_covfilt.vcf`**


----


### 8. Generate consensus with MIA and run contamMix
Step 8 uses MIA to generate a consensus file. To do this we use bedtools to transform the rescaled.bam file back to fastq. This ensures we only have the high quality,
uniq reads without duplicates for the consensus. The reads get reassembled to the reference. MIA requires a .maln file to be produced. I use the following
parameters: -r reference sequence, -f input fastq or fasta file, -c means reference/assembly is circular, -C collapse sequences with same start, end, strand info
into a single sequence, -U fasta database has repeat sequences, keep one based on sum of q-scores, -i iterate assembly until convergence (default), -F only output 
the FINAL assembly, each iteration -k use kmer filter with kmers of this length (parameter used by other researchers), -m file name for output file. Then I use ma 
to create the consensus Parameters: -M input file -f output format (see ma man page). '?' is wild card for any number maln file,  since we do not know what it will 
be. -I makes the consensus have same name as sample. For more information on MIA see: https://github.com/mpieva/mapping-iterative-assembler/tree/master/man
**Final output at this stage is a consensus fasta file: `Sample.q30.chrM.rescaled.mia.consensus.fasta`**

Next I use Philip Jonsson's program contamMix to estimate contamination along the mitochondrial genome compared to a dataset of 311 contaminant sequences from
around the world `mt311.fasta`. As input we need the consensus file, the fastq file of rescaled reads and the fasta file with contaminants. I concatenate
the sample consensus to the contaminant dataset and use mafft to perform a multiple alignment. I then remap the fastq reads generated from the rescaled bam
to the sample's own consensus using same parameters for original mapping with BWA and filtering with samtools. Lastly we run contamMix which takes as input
the remapped.bam, the mafft alignmnet between sample and contaminants, and the sample ID (--consId). 
**Final output at this stage is a pdf with contamMix analysis visuals as well as stats.txt file: `Sample.contamMix_fig`**

----


### 9. Generate postmortem damage profiles with PMDtools
Step 9 uses the python program pmdtools, developed by Pontus Skoglund, to evaluate postmortem damage (PMD) for all samples at two different thresholds (0 and 3). 
The more PMD a sample has the more ancient it likely is although this can vary depending on the amount of reads (i.e. the less reads the less useful the analysis is).
Here I use the *uniq.bam file as the *rescaled.bam file cannot be modified for the MD tag. So I use samtools to add an MD tag and then 
produce a distribution of PMD scores. I use Rscript `pmd_hist_v3.R` to produce a histogram based on these scores. I then run pmdtools at thresholds 0 and 3. 
This produces a new bam file with only those reads that have reads that pass the given thresholds. Append information on reads excluded, passed and % data passed by threshold to the output from R. 
Lastly I use the provided Rscript `plotPMD.R` to generate deamination plots at each threshold: no threshold, PMD = 0, PMD = 3. 
For more on pmdtools and to obtain the program go to: https://code.google.com/archive/p/pmdtools/ **Multiple outputs at this stage include histogram.pdf file, deamination plot pdfs, 
bamfiles at each threshold and stats.txt file. All store in: `pmdanalysis.Sample`**. In this step I also use a check to make sure to only run this analysis on
samples with >0 mapped, unique reads.


----

### 10. Generate summary tab delimited file
All the tab delimited files produced at each step are concatenated using paste. The final file called **MT.DataAnalysis.Summary.txt** provides summary of analyses
outputs and can be opened in Excel. Some statistics need to be calculated manually or in Excel such as % endogenous. These are noted in the file.

---- 

### This is what your parent directory structure should look like after pipeline has been run.
![Directory_Screenshot](https://github.com/mnievesc/Ancient_mtDNA_Pipeline/blob/master/Screenshot_ParentDirectory.png)

