# WGS

## Step 1
From raw sequence data to SNPs -> [Sequence processing](https://github.com/Cpetak/WGS/blob/main/sequence_processing.md)

## Step 2
From temporal pH and temperature data to scaled mean and variability -> [Environmental data processing](https://github.com/Cpetak/WGS/blob/main/env_data_processing.R)

## Step 3
Converting the output of ANGSD with the SNPs into the format Bayenv expects -> [Converting for Bayenv](https://github.com/Cpetak/WGS/blob/main/bayenv_conversion.py)

## Step 4
Estimating covarience matrix and Bayenv factors -> [Running Bayenv](https://github.com/Cpetak/WGS/blob/main/Running_bayenv.md)

## Step 5
Generating dataframe with Bayenv factor, number of neighbors, expression level and SNP number -> [Running Bayenv](https://github.com/Cpetak/WGS/blob/main/get_gene_from_snp.ipynb)
