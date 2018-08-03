source("http://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
library(SNPRelate)
library(tidyverse)

getwd()
setwd("/Users/kevinxchan/Documents/UBC/YEAR3/biol525D/Topic_8-9/")

#Set up file names
vcf_filename<- c("biol525d.snps.vcf.gz")
gds_filename<- c("biol525d.snps.gds")
sampleinfo_filename <- c("biol525d.popinfo.txt")
#Convert your vcf to gds for use with snprelate
snpgdsVCF2GDS(vcf_filename, gds_filename,  method="biallelic.only",ignore.chr.prefix="Chr")

#Load the gds file
genofile <- snpgdsOpen(gds_filename)
#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

snpset.id <- unlist(snpset_pruned)
#Run the PCA
pca <- snpgdsPCA(genofile, num.thread = 2, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.10, maf = 0.05,autosome.only = F)

#Lets take a look at the percent variance explained
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#Load your sample information for plotting purposes.
sampleinfo <- read_tsv(sampleinfo_filename)

#Make a dataframe of your PCA results
tab <- data.frame(name = pca$sample.id,
                  EV1 = pca$eigenvect[,3],    # the first eigenvector
                  EV2 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)

#Merge the sampleinfo into that
tab <- merge(tab, sampleinfo)

#Plot a PCA image
ggplot(data=tab,aes(EV1,EV2)) + geom_point()
#There are only a few points because so many of the samples are identical once we've pruned it down to the 12 SNPs. To check to make sure that the points are actually there, we can add a slight jitter to the positions.
ggplot(data=tab,aes(EV1,EV2)) + geom_jitter(width=0.01,height=0.01)
#Next lets color code by population and add axis labels
ggplot(data=tab,aes(EV1,EV2)) + geom_jitter(aes(color=as.factor(population)),width=0.01,height=0.01) + ylab("Principal component 2") + xlab("Principal component 1")
#We can make that look nicer
plot <- ggplot(data=tab,aes(EV1,EV2)) + geom_jitter(aes(color=as.factor(lat)),
                                            width=0.01,height=0.01) + ylab("Principal component 2") + xlab("Principal component 1") +
  theme_classic() + scale_color_discrete(name="Population") +
  theme(panel.border = element_rect(fill = NA, colour = "grey50")) +
  ggtitle("PCA plot") + theme(plot.title = element_text(hjust = 0.5))
pdf("pca.pdf")
plot(plot)
dev.off()







