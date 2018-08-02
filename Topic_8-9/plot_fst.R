install.packages("qqman", repos = "http://cran.us.r-project.org")
library(qqman)
library(tidyverse)

getwd()
setwd("/Users/kevinxchan/Documents/UBC/YEAR3/biol525D/Topic_8-9/")

fst.filename <- "biol525D.weir.fst"
data <- read.delim(fst.filename, header=T)
#Now one problem with plotting this is that the chromosomes are not intergers
summary(data$CHROM)

#This strips the "group" from the chromosome name
data$CHROM <- gsub("Chr", "", data$CHROM)

#This converts it to numeric
data$CHROM <- as.numeric(data$CHROM)

#Lets also remove values that are NA
data %>% filter(WEIR_AND_COCKERHAM_FST != "NaN") -> data

#It's important to make sure that the chromosomes are numeric instead of character
data$CHROM <- as.numeric(data$CHROM)

#Lets do a basic plot using the manhattan tool in qqman. This is generally designed for plotting pvalues from GWAS, but it works here.
manhattan(data, chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST", snp="POS",
          logp=FALSE,ylab="Fst")
