library(reshape)
library(tidyverse)

getwd()
setwd("/Users/kevinxchan/Documents/UBC/YEAR3/biol525D/Topic_8-9/")

sampleinfo <- read_tsv("biol525d.popinfo.txt")
name <- "biol525d"
#We're going to loop through each k value, so we need a dataframe to save those values
data.full <- data.frame(name=character(),
                        population=character(),
                        lat=numeric(),
                        long=numeric(),
                        variable=character(),
                        value=numeric(),
                        k=numeric())
#Now we loop through each K output
for (k in 1:4){
  #Load the data file
  data <- read.table(paste(name, k, "meanQ",sep="."),colClasses="numeric")
  #Label the columns (one for each group)
  Qnames <-paste("Q",seq(1:k),sep = "")
  colnames(data) <- Qnames
  #Add sample info to Q values
  data.info <- cbind(data, sampleinfo)
  #Melt the data into long format
  data.melt <- melt(data.info, id.var=c("name","population","lat","long"))
  #We have to make sure to include the K value for each file
  data.melt$k <- k
  #Now rbind it to the data frame outside the loop
  data.full <- rbind(data.full, data.melt)
}

#Lets try plotting for k=2
data.full %>% filter(k == 2) %>% #This selects only the k=2 lines out of the full set
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack")

#facet wrap
data.full %>%
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack") +
  facet_wrap(~k, nrow=2)

#How about if we want to order samples by the reverse of latitude.
data.full %>%
  mutate(name = factor(name, levels = name[order(-lat)])) %>%
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack") +
  facet_wrap(~k, nrow=2)

#The order works, but lets try to make it look nicer

data.full %>%
  mutate(name = factor(name, levels = name[order(-lat)])) %>%
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack") +
  facet_wrap(~k, nrow=5) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        strip.background = element_blank()) +
  ylab("q-value")+
  guides(fill = guide_legend(title = "Group", title.position = "left"))
