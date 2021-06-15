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

Fist, fix all topped partial files such as
```
#!/bin/bash
pfiles=$(ls | grep "topped")

 for i in $pfiles
 do
    echo $i
    sed -i 's/NA/NA\/NA/g' $i
 done
```
because outflank expects NA/NA instead of just NA

Then, make a separate R script for each partial vcf to analyse, in partial_vcfs, which uses outflank to calculate Fst for each site
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
-> example output: topped_xaa.vcfscript.R
Pop.txt is a text file where the first line is the word "pop", then each line contains the population number, here it is just 20 1s, 20 2s, etc.
```
ls | grep "vcfscript" > scripts.txt

```

Next, launch each R script as a separate job, the input to the following code is scripts.txt

```
#!/bin/sh
while read line ; do
        FILE=$(mktemp)
        cat header.txt >> $FILE
        echo "Rscript $line" >> $FILE
        sbatch $FILE
        sleep 0.5
        rm $FILE
done < $1

```
-> example output: topped_xaa.vcfdata.csv
Combining small csvs together after taking out header:
```
mkdir csvs
mv *data.csv csvs
cd csvs
for i in $(ls); do sed '1d' $i > ${i}_fixed; done
cat *csv_fixed > combined.csv
cut -f 2-10 -d , combined.csv | nl -w 1 -p -s , > fixed_combined.csv #reindexing csv
vim fixed_combined.csv -> insert as first line: "","LocusName","He","FST","T1","T2","FSTNoCorr","T1NoCorr","T2NoCorr","meanAlleleFreq" 
awk -F, '$3 > 0.1' fixed_combined.csv > fixed_combined_goodhe.csv #getting only sites with He > 0.1
cat fixed_combined_goodhe.csv | cut -d ',' -f2,7 > twocol.csv #keeping only position and FSTNoCorr columns
sort -k 2 -t , -n -r twocol.csv > sorted_twocol.csv #sort by Fst
cat sorted_twocol.csv | grep -v "e" > test.csv remove first few lines that were incorrectly sorted due to eg e-10
head -6823 test.csv > reduced_sorted_twocol.csv #keep only high Fst pos, the threshold will depend on the distribution
```
(note to self: number of outliers slightly different from before, whereas old fixed_combined and new are the same -> must have been an error during processing, specifically due to rounding, now it is correct)
reduced_sorted_twocol.csv is the input to this python notebook -> https://colab.research.google.com/drive/1KBItMrL52mDsmvZ_2efnrm9F-XaXS459?usp=sharing
output: list of genes or genes in which high Fst SNPs fall

## using OutFlank to plot Fst distribution
Inside new_results/partial_vcfs/csvs/results

```
library(OutFLANK)
library(vcfR)

FstDataFrame <- read.csv(file = 'fixed_combined.csv', header=TRUE,row.names=1)

print("read file")

reduced_df<-FstDataFrame[seq(1,nrow(FstDataFrame),1000),] #print every 1000th datapoint, otherwise the output is huge

print("reduced file")

pdf("line.pdf")
plot(reduced_df$FST, reduced_df$FSTNoCorr, xlim=c(-0.01,0.3), ylim=c(-0.01,0.3), pch=20)
abline(0,1)
dev.off()

print("first plot finished")

pdf("dots.pdf")
plot(reduced_df$He, reduced_df$FSTNoCorr, pch=20, col="grey")
dev.off()

print("second plot finished")

pdf("hist.pdf")
hist(reduced_df$FSTNoCorr[reduced_df$He>0.1],xlim=c(0,0.3), breaks=50)
dev.off()

print("third plot finished")
```
