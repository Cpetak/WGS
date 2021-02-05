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

k <- max(ind$pop)

FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind$pop)

#SAVE TO FILE BEFORE PROCEEDING
write.csv(FstDataFrame, file = "my_example_data.csv")
