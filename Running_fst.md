# Running Fst

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
