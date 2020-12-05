# This is a markdown to report all steps taken to get from the raw sequence data to Bayenv factors

### Checking quality of sequencing data
```
pip install multiqc
spack load fastqc@0.11.7

-----------
#!/bin/bash

yourfilenames=`ls /users/c/p/cpetak/WGS/all_fastqs/18170X*.fastq`

for file in $yourfilenames

do
	fastqc $file -o /users/c/p/cpetak/WGS/fastqc_output/
done
-----------

cd /users/c/p/cpetak/WGS/fastqc_output
multiqc .

```

### Mapping to the reference genome

```
spack load bwa@0.7.17
spack load samtools@1.10
bwa index GCF_000002235.5_Spur_5.0_genomic.fna

-----------
while read line ; do
        F1=$(cut -d ' ' -f1 <<< $line)
        F2=$(cut -d ' ' -f2 <<< $line)
        echo "$F1 -- $F2"
        FILE=$(mktemp)
        cat header.txt >> $FILE
        echo "spack load samtools@1.10" >> $FILE
        echo "spack load bwa@0.7.17" >> $FILE
        ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"
        out_name=$(cut -d '.' -f1 <<< $F1)
        echo "bwa mem -t 1 -M $ref /users/c/p/cpetak/WGS/all_fastqs/$F1 /users/c/p/cpetak/WGS/all_fastqs/$F2 | samtools view -S -b > /users/c/p/cpetak/WGS/BWA_out/$out_name.bam" >> $FILE
          sbatch $FILE
          sleep 0.5
          rm $FILE
done < $1
-----------
```
#### Checking mapping statistics
```
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools sort /users/c/p/cpetak/WGS/BWA_out/$line -o /users/c/p/cpetak/WGS/BWA_out/$out_name.sorted.bam" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools rmdup /users/c/p/cpetak/WGS/BWA_out/$line /users/c/p/cpetak/WGS/BWA_out/$out_name.rmdup.bam" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools flagstat /users/c/p/cpetak/WGS/BWA_out/$line | awk 'NR>=6&&NR<=13 {print \$1}' | column -x >> /users/c/p/cpetak/WGS/$out_name.flagstats.txt" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools depth /users/c/p/cpetak/WGS/BWA_out/$line | awk '{sum+=\$3} END {print sum/NR}' >> /users/c/p/cpetak/WGS/$out_name.coverage.txt" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
```

### Getting genotype likelihoods and SNPs using ANGSD

```
#RUN THIS CODE FOR EACH POP
#INPUT POP_NAME.rmdup.bam

cd /users/c/p/cpetak/WGS/angsd

ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"

./angsd -b /users/c/p/cpetak/WGS/FOG_rmdups_jo.txt \
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/angsd_new/FOG_angsd_allsites \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 17 \ #85%
-setMinDepthInd 3 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-doHWE 1 \
# -SNP_pval 1e-6
```

```
#RUN THIS CODE FOR ALL POPS AT ONCE FOR PCA
#FILTERING IS STRICTER HERE TO REDUCE DATA FOR PCA

cd /users/c/p/cpetak/WGS/angsd

ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"

./angsd -b /users/c/p/cpetak/WGS/all_rmdups_jo.txt \
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 119 \ #85%
-setMinDepthInd 4 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doGlf 2 \
-SNP_pval 1e-6

-----------

python /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/pcangsd_covmatrix -threads 16

-----------
#R
C <- as.matrix(read.table("pcangsd_covmatrix.cov"))
ids <- read.table("~/Downloads/pca_pops.txt")
e <- eigen(C)
# base R
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2", bg=ids$V1, pch=21)
#ggplot
library(ggplot2)
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])

df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()
```
