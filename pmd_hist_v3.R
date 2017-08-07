#!/usr/bin/env Rscript

# R script for examining PMD scores in sequenced sample 
# Written by M. Nieves Colon - May 11, 2015
# Modified - Jul 19, 2016 to take commands from terminal 

# This script is to be used after running PMDtools v 0.5 as follows
# samtools view sample.md.bam | python pmdtools.py -p > pmdscore.txt

# Distribution of pmd scores has four tab delimited columns, column 4 is the pmdscore
# see Skoglund 2014 paper for explanation (p.2234), also see this line from pmdtools.py code
# if options.printDS:
#		print L_D,'\t',L_M,'\t',L_D/L_M,'\t',LR#,'\t',readlen,'\t',perc_identity,'\t',perc_identity*(math.log((L_D/L_M)))
# col 4 is the pmd ratio, columns 5 + are not present unless we choose % identity option 

# Import argument from command line: sample name
args=commandArgs(trailingOnly=TRUE)


# Mode function  
Mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
 }


# 1. Read in the pmd scores as a table
pmd=read.table(paste(args[1],".pmdscores.txt", sep=""), header=F)
head(pmd)



# 2. To understand how the values are distributed we can obtain the min value, the max value, 
#    the mean and the mode. Append these to a dataframe and then output as a tab delimited text file.

valuevec=c(args[1],min(pmd$V4),max(pmd$V4),mean(pmd$V4),median(pmd$V4),Mode(pmd$V4))
df=as.data.frame(rbind(valuevec))
#colnames(df)=c("Min pmd", "Max pmd", "Mean pmd", "Median pmd", "Mode pmd")
write.table(df, file=paste(args[1],".PMDscore.dist.temp.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



# 3. Make a boxplot 

pdf(paste(args[1],".PMDscore_dist.pdf", sep=""))
boxplot(pmd$V4, main=paste("PMD score distribution for sample ",args[1]," (Q30 reads)", sep=""), ylab="PMD scores", xlab=args[1])

# 4. Make a histogram 
hist(pmd$V4, main=paste("PMD score distribution for sample ",args[1]," (Q30 reads)", sep=""), xlab="PMD scores", ylab="Frequency", col="purple")
dev.off()

# If we want to make a better boxplot we can use the information from above for setting axes to make another plot and append to the pdf file.
#Set values for axes
#minpmd=ceiling(min(pmd$V4)-1)
#maxpmd=ceiling(max(pmd$V4)+1)

#pdf("PMDscore_dist_PI437D_hist2.pdf")
#hist(pmd$V4, main="PMD score distribution for sample PI437D (Q25 reads)", axes=F, xlab="PMD scores", ylab="Frequency", col="red")
#axis(1,at=seq(minpmd,maxpmd,5)) # x axis
#axis(2, at=seq(0,1000000,200000)) # y axis
#dev.off()


## References
# 1. http://www.statmethods.net/advgraphs/axes.html
# 2. http://www.inside-r.org/r-doc/graphics/hist
# 3. http://www.dummies.com/how-to/content/how-to-round-off-numbers-in-r.html
# 4. http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode
# 5. http://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
# 6. http://www.r-bloggers.com/paste-paste0-and-sprintf/
