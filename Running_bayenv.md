# Running Bayenv

### Estimating covariance matrix
```
cd /users/c/p/cpetak/bayenv2_public
./bayenv2 -i /users/c/p/cpetak/WGS/angsd_new_noout/snpsfile_100 -p 6 -k 100000 -r 63479 > /users/c/p/cpetak/WGS/angsd_new_noout/matrix_100.out
```
### Creating an individual file for each SNP
```
mkdir temp
cd temp

split -a 10 -l 2 $SNPFILE snp_batch
```

### Running Bayenv for each SNP and appending output to single file
```
#!/bin/bash

ENVFILE=$1
MATFILE=$2
POPNUM=$3
ITNUM=$4
ENVNUM=$5

for f in $(find . -type f -name "use2_snp*")
do
./bayenv2 -i $f -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -o filter_100
done
```
ITNUM was set to 10.000.
