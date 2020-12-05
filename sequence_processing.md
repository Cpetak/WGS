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

### Processing environmental data

```
#R
library(ggplot2)
library(varhandle)
library(arules)
library(dplyr)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

setwd("/Users/csengepetak/Dropbox/Research/ph_data/ph_data_Mooring")
phdata <- read.csv(file="Fogarty_Creek.csv", header=FALSE)
phdata <- phdata[, c("V1","V2","V4","V5")]
colnames(phdata) <- c("Location","Date","Temp","pH")
phdata[,1] <- "A/Fogarty_Creek"

csvs <- c("Cape_Blanco.csv", "Kibesillah_Hill.csv", "Bodega.csv", "Terrace_Point.csv", "Lompoc_Landing.csv")
locs <- c("D/Cape_Blanco", "F/Kibesillah_Hill", "H/Bodega", "J/Terrace_Point", "M/Lompoc_Landing")

i <- 1

setwd("/Users/csengepetak/Dropbox/Research/ph_data/ph_data_Mooring/csvs")

for (f in csvs){
  temp <- read.csv(file=f, header=FALSE)
  temp <- temp[, c("V1","V2","V4","V5")]
  temp[,1] <- locs[i]
  colnames(temp) <- c("Location","Date","Temp","pH")
  phdata <- rbind(temp, phdata)
  i <- i+1
}

#phdata$pH <- unfactor(phdata$pH)
phdata$pH <- as.numeric(phdata$pH)
phdata <- na.omit(phdata)

phdata_2011 <- subset(phdata, Date < 20120000)
phdata_2011[phdata_2011$Date < 20110500, "Date"] <- "4"
phdata_2011[phdata_2011$Date < 20110600 & phdata_2011$Date >= 20110500, "Date"] <- "5"
phdata_2011[phdata_2011$Date < 20110700 & phdata_2011$Date >= 20110600, "Date"] <- "6"
phdata_2011[phdata_2011$Date < 20110800 & phdata_2011$Date >= 20110700, "Date"] <- "7"
phdata_2011[phdata_2011$Date < 20110900 & phdata_2011$Date >= 20110800, "Date"] <- "8"
phdata_2011[phdata_2011$Date < 20111000 & phdata_2011$Date >= 20110900, "Date"] <- "9"
#phdata_2011$Date<- as.Date(as.character(phdata_2011$Date) , format = "%Y%m%d",origin = "19601001")


ggplot(phdata_2011, aes(x=Date, y=pH, fill=Location)) + 
  geom_boxplot()
ggplot(phdata_2011, aes(x=Location, y=pH)) + 
  geom_boxplot()

tgc <- summarySE(phdata_2011, measurevar="pH", groupvars=c("Location"))
ggplot(tgc, aes(x=Location, y=pH)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=pH-ci, ymax=pH+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(7.8,8.2))

res.aov <- aov(pH ~ Location, data = phdata_2011)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

phdata_2012 <- subset(phdata, Date < 20130000 & Date > 20120000)
phdata_2012[phdata_2012$Date < 20120500, "Date"] <- "4"
phdata_2012[phdata_2012$Date < 20120600 & phdata_2012$Date >= 20120500, "Date"] <- "5"
phdata_2012[phdata_2012$Date < 20120700 & phdata_2012$Date >= 20120600, "Date"] <- "6"
phdata_2012[phdata_2012$Date < 20120800 & phdata_2012$Date >= 20120700, "Date"] <- "7"
phdata_2012[phdata_2012$Date < 20120900 & phdata_2012$Date >= 20120800, "Date"] <- "8"
phdata_2012[phdata_2012$Date < 20121000 & phdata_2012$Date >= 20120900, "Date"] <- "9"
phdata_2012[phdata_2012$Date < 20121100 & phdata_2012$Date >= 20121000, "Date"] <- NA
phdata_2012 <- na.omit(phdata_2012)

ggplot(phdata_2012, aes(x=Date, y=pH, fill=Location)) + 
  geom_boxplot()
ggplot(phdata_2012, aes(x=Location, y=pH)) + 
  geom_boxplot()

tgc <- summarySE(phdata_2012, measurevar="pH", groupvars=c("Location"))
ggplot(tgc, aes(x=Location, y=pH)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=pH-ci, ymax=pH+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(7.8,8.2))

res.aov <- aov(pH ~ Location, data = phdata_2012)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

phdata_2013 <- subset(phdata, Date < 20140000 & Date > 20130000)
phdata_2013[phdata_2013$Date < 20130500, "Date"] <- "4"
phdata_2013[phdata_2013$Date < 20130600 & phdata_2013$Date >= 20130500, "Date"] <- "5"
phdata_2013[phdata_2013$Date < 20130700 & phdata_2013$Date >= 20130600, "Date"] <- "6"
phdata_2013[phdata_2013$Date < 20130800 & phdata_2013$Date >= 20130700, "Date"] <- "7"
phdata_2013[phdata_2013$Date < 20130900 & phdata_2013$Date >= 20130800, "Date"] <- "8"
phdata_2013[phdata_2013$Date < 20131000 & phdata_2013$Date >= 20130900, "Date"] <- "9"
phdata_2013[phdata_2013$Date < 20131100 & phdata_2013$Date >= 20131000, "Date"] <- NA
phdata_2013 <- na.omit(phdata_2013)

ggplot(phdata_2013, aes(x=Date, y=pH, fill=Location)) + 
  geom_boxplot()
ggplot(phdata_2013, aes(x=Location, y=pH)) + 
  geom_boxplot()

tgc <- summarySE(phdata_2013, measurevar="pH", groupvars=c("Location"))
ggplot(tgc, aes(x=Location, y=pH)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=pH-ci, ymax=pH+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(7.8,8.2))

res.aov <- aov(pH ~ Location, data = phdata_2013)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

phdata_2011$Year <- "2011"
phdata_2012$Year <- "2012"
phdata_2013$Year <- "2013"

all_ph <- rbind(phdata_2011, phdata_2012, phdata_2013)

ggplot(all_ph, aes(x=Location, y=pH)) + 
  geom_boxplot() 
ggplot(all_ph, aes(x=Date, y=pH, fill=Location)) + 
  geom_boxplot()

tgc <- summarySE(all_ph, measurevar="Temp", groupvars=c("Location"))
temp_means <- tgc["Temp"]
tgc <- summarySE(all_ph, measurevar="pH", groupvars=c("Location"))
ph_means <- tgc["pH"]

ggplot(tgc, aes(x=Location, y=Temp)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Temp-ci, ymax=Temp+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(9,13)) +
  ylab("Mean temperature")

res.aov <- aov(Temp ~ Location, data = all_ph)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

#geom_ribbon(aes(ymax = pH + sd, ymin = pH - sd, fill= Location))
test <- all_ph %>% group_by(Location) %>% group_map(~ sum(.x$pH < 7.8)/length(.x$pH))
plot(seq(length(test)),test)

newdf <- data.frame("temp" = 1:6)

newdf$xval <- seq(length(test))
newdf$yval <- c(0.1780976, 0.07188418, 0.003375068, 0.05026657, 0,0.00729927 )

ggplot(newdf, aes(x=temp, y=yval)) +
  geom_point(size=3,color="blue") +
  ylab("Frequency of pH under 7.8") +
  xlab("Location") +
  scale_x_continuous(breaks=seq(1, 6, 1), labels = c("Fogarty_Creek", "Cape_Blanco", "Kibesillah_Hill", "Bodega", "Terrace_Point", "Lompoc"))

test2 <- all_ph %>% group_by(Location) %>% group_map(~ sd(.x$Temp))
plot(seq(length(test2)),test2)

newdf <- data.frame("temp" = 1:6)

newdf$xval <- seq(length(test))
newdf$yval <- c(1.713702,1.344381, 1.205099, 1.549824,1.210292 , 1.098031)

ggplot(newdf, aes(x=temp, y=yval)) +
  geom_point(size=3,color="blue") +
  ylab("Standard deviation of temperature") +
  xlab("Location") +
  scale_x_continuous(breaks=seq(1, 6, 1), labels = c("Fogarty_Creek", "Cape_Blanco", "Kibesillah_Hill", "Bodega", "Terrace_Point", "Lompoc"))

Fy <- subset(phdata, Location == "M/Lompoc_Landing")
sd(Fy$pH)
min(Fy$pH)
mean(Fy$pH)
Sh <- subset(phdata, Location == "D/Cape_Blanco")
sd(Sh$pH)

ph_freq <- c(0.1780976, 0.07188418, 0.003375068, 0.05026657, 0,0.00729927 )
temp_var <- c(1.713702,1.344381, 1.205099, 1.549824,1.210292 , 1.098031)
ph_mean <- c(7.983623,8.009176,8.022380,7.996050,8.177942,8.076118)
temp_mean <- c(10.35057,10.17522,10.15551,11.41514,12.67584,12.47098)

final <- as.data.frame(ph_freq)
final["temp_var"] <- temp_var
final["ph_mean"] <- ph_mean
final["temp_mean"] <- temp_mean
df.scaled <- as.data.frame(scale(final))
locs <- c("Fogarty Creek","Cape Blanco", "Kibesilah Hill", "Bodega Bay", "Terrace Point", "Lompoc Landing")
df.scaled["locs"] <- locs

level_order <- c("Fogarty Creek","Cape Blanco", "Kibesilah Hill", "Bodega Bay", "Terrace Point", "Lompoc Landing")

ggplot(df.scaled, aes(x= factor(locs, level = level_order),group = 1)) + 
  geom_line(aes(y = temp_var,color = "darkred")) +
  geom_point(aes(y = temp_var),colour = 'darkgrey', size = 2) +
  geom_line(aes(y = temp_mean,color="steelblue")) +
  geom_point(aes(y = temp_mean),colour = 'darkgrey', size = 2) +
  geom_line(aes(y = ph_mean,color="orange")) +
  geom_point(aes(y = ph_mean),colour = 'darkgrey', size = 2) +
  geom_line(aes(y = ph_freq,color="darkgreen")) +
  geom_point(aes(y = ph_freq),colour = 'darkgrey', size = 2)+
  scale_color_identity(name = "Environmental variable",
                       breaks = c("darkgreen", "darkred", "orange","steelblue"),
                       labels = c("Freq. pH < 7.8", "Temp var", "pH mean", "Temp mean"), guide = "legend")+
  xlab("Location")+
  ylab("Scaled value")+
  theme_bw()

```
