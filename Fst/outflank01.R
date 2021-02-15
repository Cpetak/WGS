#for i in $(ls); do sed '1d' $i > ${i}_fixed; done
#cat *csv_fixed > combined.csv
#example of how to use outflank. copied from https://rpubs.com/lotterhos/outflank
library(OutFLANK)
library(vcfR)

#inputs I need to provide. example files were downloaded from github.
vcf <- read.vcfR("sim1a.vcf", verbose=FALSE)
ind <- read.table("Pop.txt", header=TRUE) #only ind$pop is used

#Convert VCF format to SNP data format
convertVCFtoCount3 <- function(string){
  # This function assumes 0 for reference
  # and 1 for alternate allele
  a <- as.numeric(unlist(strsplit(string, split = c("[|///]"))))
  odd = seq(1, length(a), by=2)
  a[odd] + a[odd+1]
}

#setting up
all.vcf.gen <- vcf@gt[,-1]
system.time(gen_table <- matrix(convertVCFtoCount3(all.vcf.gen), ncol=ncol(all.vcf.gen)))

locinames <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_")

SNPdata <- t(gen_table)
temp <- SNPdata[rowSums(is.na(SNPdata))>0,]
length(temp)


k <- max(ind$pop)

SNPdata[is.na(SNPdata)] <- 9

FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind$pop)

#SAVE TO FILE BEFORE PROCEEDING
write.csv(FstDataFrame, file = "my_example_data.csv")
#uncomment following if needed
FstDataFrame <- read.csv(file = 'my_data.csv', header=TRUE,row.names=1)

#subset data for plotting -otherwise too many points crush my computer
reduced_df<-FstDataFrame[seq(1,nrow(FstDataFrame),100000),]

#should be on a line
pdf("line.pdf")
plot(reduced_df$FST, reduced_df$FSTNoCorr, 
     xlim=c(-0.01,0.3), ylim=c(-0.01,0.3),
     pch=20)
abline(0,1)
dev.off()

#example of how a SNP with missing data would look like that should be removed
#shown with star on plot
#SNPdata_missing <- SNPdata
#missing <- sample(1:nrow(SNPdata_missing), 500, replace=FALSE)
#SNPdata_missing[missing,1] <- 9
#FstDataFrame_missing <- MakeDiploidFSTMat(SNPdata_missing,locinames,ind$pop)
#plot(FstDataFrame_missing$FST, FstDataFrame_missing$FSTNoCorr, xlim=c(-0.01,0.3), ylim=c(-0.01,0.3), pch=20)
#points(FstDataFrame_missing$FST[1], FstDataFrame_missing$FSTNoCorr[1], col="blue", pch=8, cex=1.3)
#abline(0,1)

pdf("dots.pdf")
plot(reduced_df$He, reduced_df$FSTNoCorr, pch=20, col="grey") #outliers
# Note the large FST values for loci with low heterozygosity (He < 0.1)
dev.off()

pdf("hist.pdf")
#demonstration of how we need to remove low heterozygocity loci 
#hist(reduced_df$FSTNoCorr, breaks=50)
hist(reduced_df$FSTNoCorr[reduced_df$He>0.05], breaks=50)
#this is how the distribution should look like instead:
#hist(FstDataFrame$FSTNoCorr[FstDataFrame$He>0.1], breaks=seq(0,0.3, by=0.001))
dev.off()

#FstDataFrame_corr <- FstDataFrame[FstDataFrame$He>0.1,] don't need this as it's part of next step
#By default, OutFLANK removes from consideration all loci with expected heterozygosity less than 10%.

#check if there is any NA. They come to be when everyone in the vcf if heterozygous to a location.
temp <- FstDataFrame[rowSums(is.na(FstDataFrame))>0,] #967
#len of temp is = number of rows where everyone is heterozygous to the locus
#FstDataFrame2 <- na.omit(FstDataFrame)
#FstDataFrame3 <- FstDataFrame2[!(FstDataFrame2$FST==0),] #this and the next line were written to deal with a bug but in the end adjusting the righttrimfraction was the solution so these might not be useful
#row.names(FstDataFrame3) <- (1:nrow(FstDataFrame3))

k=7

#running OutFlank
out1 <- OutFLANK(FstDataFrame, NumberOfSamples=k)#RightTrimFraction = 0.01) #see Rstudio "help" for options to this function
OutFLANKResultsPlotter(out1, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL) #see Rstudio "help" for options to this function
#this is a good fit. but if not>
#tweak qthreshold, RightTrimFraction

#looking at outliers
outlier = OutFLANK(FstDataFrame,NumberOfSamples = k, 
                   RightTrimFraction = 0.06, LeftTrimFraction = 0.35,
                   qthreshold = 0.05, Hmin = 0.1)
sum(outlier$results$qvalues<0.01, na.rm=TRUE)
plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
points(outlier$results$He[outlier$results$qvalues<0.01], outlier$results$FST[outlier$results$qvalues<0.01], pch=21, col="blue")

top_candidates <- outlier$results$qvalues<0.01 & outlier$results$He>0.1
topcan <- outlier$results[top_candidates,]
topcan[order(topcan$LocusName),]

write.csv(topcan, file = "top_fst.csv")

plot(vcf@fix[,"POS"], outlier$results$FST, 
     col=grey((as.numeric(vcf@fix[,"CHROM"])%%2+1)/3),
     pch=20
) 
points(vcf@fix[top_candidates,"POS"], outlier$results$FST[top_candidates], 
       pch=21, cex=2, col=2)

# clean up by removing low freq variants
keep <- outlier$results$He>0.1 & !is.na(outlier$results$He)
plot(vcf@fix[keep,"POS"], outlier$results$FST[keep], 
     col=grey((as.numeric(vcf@fix[keep,"CHROM"])%%2+1)/3),
     pch=20
) 
points(vcf@fix[top_candidates,"POS"], outlier$results$FST[top_candidates], 
       pch=21, cex=2, col=2)
