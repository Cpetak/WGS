# Running Fst

## Make vcf file
### Make bcf using ANGSD
```

ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"

./angsd -b /users/c/p/cpetak/WGS/make_vcf/list_for_angsd.txt \
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/make_vcf/all_pop_angsd \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 119 \
-setMinDepthInd 3 \
-skipTriallelic 1 \
-dobcf 1 \
-GL 1 \
-doPost 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-doHWE 1 \
-SNP_pval 1e-6


```
### bcf to vcf
```
bcftools view all_pop_angsd.bcf > all_pop_angsd.vcf
```

### adding GT and preparing for OutFlank
Adding GT:
```
vcfglxgt all_pop_angsd.vcf > fixed_all_pop_angsd.vcf
```
If all values in vcf (GT:GP:PL:GL:DP) are 0 -> due to missing data, replace with NA
```
sed -i 's/0\/0\:0\,0\,0\:0\,0\,0\:0\,0\,0\:0/NA\:0\,0\,0\:0\,0\,0\:0\,0\,0\:0/g' fixed_all_pop_angsd_copy.vcf
```
Keep only GT information
```
awk -v FS="\t" -v OFS="\t" '{for(i=9;i<=NF;i++) {split($i, gt, ":"); $i=gt[1]} print}' fixed_all_pop_angsd_copy.vcf > fixed_all_pop_angsd_copy_onlyGT.vcf
```
Separate header a split rest for parallelisation \\
these are the files I copied that I am starting over with. \\ 
ATTENTION fixed_all_pop_angsd.vcf and fixed_all_pop_angsd_onlyGT.vcf in the make_vcf folder are the files prior the the above 2 steps!!!
```
head -916 fixed_all_pop_angsd_copy_onlyGT.vcf > vcf_head.vcf
cat fixed_all_pop_angsd_copy_onlyGT.vcf | grep -v "#" > vcf_tail.vcf
split -l 100000 vcf_tail.vcf

```
### using OutFlank to get per-site Fst

```

library(OutFLANK)
library(vcfR)

#inputs I need to provide. example files were downloaded from github.
vcf <- read.vcfR("fixed_all_pop_angsd_onlyGT.vcf", verbose=FALSE)
ind <- read.table("my_pops.txt", header=TRUE) #only ind$pop is used

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
write.csv(FstDataFrame, file = "all_fst_data.csv")
```
