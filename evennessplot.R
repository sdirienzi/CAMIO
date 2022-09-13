#to plot Renyi number over the course of the experiment to look for changes in evenness.

library("dplyr")
library("ggplot2")
library("reshape2")
library("data.table")
library("gridExtra")
library("vegan")
library("tidyr")
library("BiodiversityR")
library("stringr")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
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

#Figure 3E
data<-read.delim(file="SequencingData.txt",header=TRUE, sep="\t")
dataincombined<-subset(data,data$time!="inoculate" & data$fraction!="NA")

wide<-reshape2::dcast(dataincombined,experiment+sample+Sugar+time+generation+samplename~strain, value.var="correctedfraction2")

Renyi.2 <- renyiresult(wide[,7:12]*100,method="s")  #3:8 for MiSeq
Renyimerge<-data.table(wide, Renyi.2)
names(Renyimerge)<-c(names(wide)[1:12],"R0","R0.25","R0.5","R1","R2","R4","R8","infinite") #MiSeq

shapevalues<-c(seq(0,12,1))
names(shapevalues)<-levels(as.factor(Renyimerge$sample))

#get stats
mydata<-subset(Renyimerge, Renyimerge$time=="6" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 0.73077, df = 2, p-value = 0.6939
mydata<-subset(Renyimerge, Renyimerge$time=="7" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 0.80769, df = 2, p-value = 0.6677
mydata<-subset(Renyimerge, Renyimerge$time=="9" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 4.1923, df = 2, p-value = 0.1229

mydata<-subset(Renyimerge, Renyimerge$time=="11" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 0.96154, df = 2, p-value = 0.6183

mydata<-subset(Renyimerge, Renyimerge$time=="13" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 3.2308, df = 2, p-value = 0.1988

mydata<-subset(Renyimerge, Renyimerge$time=="15" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 2.3462, df = 2, p-value = 0.3094

mydata<-subset(Renyimerge, Renyimerge$time=="18" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 3.0385, df = 2, p-value = 0.2189


mydata<-subset(Renyimerge, Renyimerge$time=="takedown" )
kruskal.test(as.numeric(infinite) ~ Sugar, data=mydata)
# Kruskal-Wallis chi-squared = 6.9615, df = 2, p-value = 0.03078
library("DescTools")
DunnettTest(as.numeric(infinite) ~ as.factor(Sugar), data=mydata, control="Sucrose")
# diff     lwr.ci      upr.ci   pval    
# HFCS-Sucrose      -0.5531136 -0.8671140 -0.23911326 0.0024 ** 
#   Trehalose-Sucrose -0.3930804 -0.7070808 -0.07908002 0.0175 *  


#boxplot of final evenness

Renyimerge$Sugar<-factor(Renyimerge$Sugar,c("Sucrose","HFCS","Trehalose"))
#make box plot
pdf(file="Finalevennessboxplotcolored.pdf",width=4,height=4)
ggplot(subset(Renyimerge, Renyimerge$time=="takedown" ), aes(x=Sugar,y=infinite,fill=Sugar,colour=Sugar)) +
  geom_boxplot(size=0.8,colour="black",alpha=0.4) +  #, lty=FAMI
  geom_point(position=position_jitterdodge(jitter.width =.7),colour="black",pch=21,size=4)+
  scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  scale_fill_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  ylab("Evenness")+
  scale_y_continuous(expand=expansion(mult= c(0.07,.3)))+
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", size=1.1),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=18,colour="black"),
    axis.title.y= element_text(size=18,colour="black"),
    axis.text.x= element_text(size=14,colour="black"),
    axis.text.y= element_text(size=14,colour="black"),
    axis.ticks = element_line(colour="black"),  
    legend.position="none"
  ) 
dev.off()


