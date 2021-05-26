# Running Fst

### Make vcf using ANGSD
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
### Creating an individual file for each SNP
```
mkdir temp
cd temp

split -a 10 -l 2 $SNPFILE snp_batch
```
