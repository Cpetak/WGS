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
these are the files I copied that I am starting over with in the new_results folder. \\ 
ATTENTION fixed_all_pop_angsd.vcf and fixed_all_pop_angsd_onlyGT.vcf in the make_vcf folder are the files prior the the above 2 steps!!!
```
head -916 fixed_all_pop_angsd_copy_onlyGT.vcf > vcf_head.vcf
cat fixed_all_pop_angsd_copy_onlyGT.vcf | grep -v "#" > vcf_tail.vcf
split -l 100000 vcf_tail.vcf
mkdir partial_vcfs
mv x* partial_vcfs/
cp vcf_head.vcf new_results/partial_vcfs
cd new_results/partial_vcfs
#!/bin/bash
pfiles=$(ls | grep "x")
for i in $pfiles
do
        cat vcf_head.vcf $i > topped_${i}.vcf
done

```
## using OutFlank to get per-site Fst
Make a separate R script for each partial vcf to analyse
```
#!/bin/sh
pfiles=$(ls | grep "topped")

for i in $pfiles
do
	echo $i
	script_name=${i}script.R
	echo -e "library(OutFLANK)" >> $script_name
	echo -e "library(vcfR)" >> $script_name
	echo -e "vcf <- read.vcfR(\"$i\", verbose=FALSE)" >> $script_name
	echo -e "ind <- read.table(\"Pop.txt\", header=TRUE)" >> $script_name
	echo -e "convertVCFtoCount3 <- function(string){" >> $script_name
	echo -e "\ta <- as.numeric(unlist(strsplit(string, split = c(\"[|///]\"))))" >> $script_name
	echo -e "\todd = seq(1, length(a), by=2)" >> $script_name
	echo -e "\ta[odd] + a[odd+1]" >> $script_name
	echo -e "}" >> $script_name
	echo -e "all.vcf.gen <- vcf@gt[,-1]" >> $script_name
	echo -e "system.time(gen_table <- matrix(convertVCFtoCount3(all.vcf.gen), ncol=ncol(all.vcf.gen)))" >> $script_name
	echo -e "locinames <- paste(vcf@fix[,\"CHROM\"], vcf@fix[,\"POS\"], sep=\"_\")" >> $script_name
	echo -e "SNPdata <- t(gen_table)" >> $script_name
  echo -e "SNPdata[is.na(SNPdata)] <- 9" >> $script_name
	echo -e "k <- max(ind\$pop)" >> $script_name
	echo -e "FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind\$pop)" >> $script_name
	echo -e "write.csv(FstDataFrame, file = \"${i}data.csv\")" >> $script_name
done
```
