require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)
library(ggsci)
library(DescTools)

#Supplemental Figure 4
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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

hormonedata<-read.delim(file="selectedanalyteshormonedata.txt",header=TRUE,sep="\t",check.names=FALSE)
limits<-read.delim(file="analytelimites.txt")
limitssel<-limits[,c(3,4,5,6,11,12,13,14,15)]

sugarcolours<-c("#6A5ACD","#FFA500","violetred1")
names(sugarcolours)<-c("SDRBRMSucrose","SDRBRMHFCS","SDRBRMTrehalose")
media.labs <- c( "HFCS","Sucrose","Trehalose")
names(media.labs) <- c("SDRBRMHFCS","SDRBRMSucrose","SDRBRMTrehalose" )
analytes<-levels(as.factor(hormonedata$Analyte))
hormonedata$media<-factor(hormonedata$media,levels=c("SDRBRMSucrose","SDRBRMHFCS","SDRBRMTrehalose"))
hormonedata<-as.data.table(hormonedata)
hormonedataReactorReplicate<-factor(hormonedata$ReactorReplicate,c("M1","R1", "R2", "R3", "M2", "R4"))

analyteplots <- list()
#excluded serotonin here, different axes needed
for (i in 1:(length(analytes)) ) {
  analyte<-analytes[i]
  if (analyte!="Serotonin") {
    analyteplots[[i]]<-ggplot(hormonedata[Analyte==analyte], aes(x=ReactorReplicate,y=Concentration,fill=media)) +    #shape=factor(media)
      geom_hline(yintercept = as.numeric(limitssel[1,i]),colour="gray")+
      geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+
      #geom_dotplot(binaxis = "y",stackdir = "center",colour="black",dotsize=0.5)+
      geom_point(position=position_jitterdodge(jitter.width =.7),size=3,pch=21)+
       facet_wrap(~media,strip.position="bottom",scales = "free_x",labeller = labeller(media = media.labs)) +
      ylab("pg/ml") +
      labs(title = analyte)+
      scale_fill_manual(name = "Media",values=sugarcolours) +
      theme_classic() +
      theme(
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.title.x= element_blank(), #element_text(size=1, colour="white"),
        axis.title.y= element_text(size=14),
        axis.text.x= element_text(size=12, colour = "black"),
        axis.text.y= element_text(size=12, colour = "black"), 
        #axis.text.x=element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        plot.title = element_text(size=14,face="bold")
      )  
  }
  else if (analyte=="Serotonin") {
    analyteplots[[i]]<-ggplot(hormonedata[Analyte==analyte], aes(x=ReactorReplicate,y=Concentration,fill=media)) +    #shape=factor(media)
      geom_hline(yintercept = as.numeric(limitssel[1,i]),colour="gray")+
      geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+
      #geom_dotplot(binaxis = "y",stackdir = "center",colour="black",dotsize=0.5)+
      geom_point(position=position_jitterdodge(jitter.width =.7),size=3,pch=21)+
      facet_wrap(~media,strip.position="bottom",scales = "free_x",labeller = labeller(media = media.labs)) +
      ylab("ng/ml") +
      labs(title = analyte)+
      scale_fill_manual(name = "Media",values=sugarcolours) +
      theme_classic() +
      theme(
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.title.x= element_blank(), #element_text(size=1, colour="white"),
        axis.title.y= element_text(size=14),
        axis.text.x= element_text(size=12, colour = "black"),
        axis.text.y= element_text(size=12, colour = "black"), 
        #axis.text.x=element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        plot.title = element_text(size=14,face="bold")
      ) 
  }
  }

pdf(file="rawhormoneplots1.pdf",width=10,height=8.5)
grid.arrange(analyteplots[[1]],analyteplots[[2]],analyteplots[[3]],ncol=1)
dev.off()

pdf(file="rawhormoneplots2.pdf",width=10,height=8.5)
grid.arrange(analyteplots[[4]],analyteplots[[5]],analyteplots[[6]],ncol=1)
dev.off()

pdf(file="rawhormoneplots3.pdf",width=10,height=8.5)
grid.arrange(analyteplots[[7]],analyteplots[[8]],analyteplots[[9]],ncol=1)
dev.off()

#Figure 5A-C
#media alone analysis
#to plot just the media results:
analyteplotsmedia <- list()
for (i in 1:(length(analytes)) ) {
  analyte<-analytes[i]
  if (analyte!="Serotonin") {
    analyteplotsmedia[[i]]<-ggplot(hormonedata[Blank=="TRUE"&Analyte==analyte], aes(x=media,y=Concentration,fill=media)) +    #shape=factor(media)  #&rawplotkeep=="yes" <-use this to exclude points below lod
      geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+
      geom_point(position=position_jitterdodge(jitter.width =.7),size=3,pch=21)+
      scale_y_continuous(expand=c(0.2,3))+
      ylab("pg/ml") +
      labs(title = analyte)+
      scale_x_discrete(labels=c("Sucrose", "HFCS","Trehalose"))+
      scale_fill_manual(name = "Media",values=sugarcolours) +
      theme_classic() +
      theme(
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.title.x= element_blank(), #element_text(size=1, colour="white"),
        axis.title.y= element_text(size=14),
        axis.text.x= element_text(size=12, colour = "black"),
        axis.text.y= element_text(size=12, colour = "black"), 
        legend.position="none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        plot.title = element_text(size=14,face="bold")
      )  
  }
}

pdf(file="mediaalonerawvalues.pdf",width=9,heigh=3)
grid.arrange(analyteplotsmedia[[2]],analyteplotsmedia[[6]],analyteplotsmedia[[7]],ncol=3)
dev.off()

kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="C.Peptide"])
#Kruskal-Wallis chi-squared = 1.0686, df = 2, p-value = 0.5861
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="GIP.total"])
#Kruskal-Wallis chi-squared = 4.4461, df = 2, p-value = 0.1083
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="GLP.1.total"])
#Kruskal-Wallis chi-squared = 0.67283, df = 2, p-value = 0.7143
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="MCP.1"])
#Kruskal-Wallis chi-squared = 1.731, df = 2, p-value = 0.4208
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="Serotonin"])
#Kruskal-Wallis chi-squared = 3.8246, df = 2, p-value = 0.1477
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="TNF.alpha"])
#Kruskal-Wallis chi-squared = 3.5387, df = 2, p-value = 0.1704
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="Ghrelin"])
#Kruskal-Wallis chi-squared = 2.4569, df = 2, p-value = 0.2927
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="PP"])
#Kruskal-Wallis chi-squared = 4.6696, df = 2, p-value = 0.09683
kruskal.test(Concentration~as.factor(media), data=hormonedata[Blank=="TRUE"&Analyte=="PYY"])
#Kruskal-Wallis chi-squared = 11.732, df = 2, p-value = 0.002834
DunnTest(Concentration~as.factor(media), data=selectedanalytes[Blank=="TRUE"&Analyte=="PYY"],method="BH")
# mean.rank.diff   pval    
# SDRBRMHFCS-SDRBRMSucrose           8.5833333 0.0072 ** 
#   SDRBRMTrehalose-SDRBRMSucrose     -0.8333333 0.7843    
# SDRBRMTrehalose-SDRBRMHFCS        -9.4166667 0.0059 ** 

analytesold<-levels(as.factor(selectedanalytes$Analyte))

#to plot just the media results:
analyteplotsoldmedia <- list()
#excluded serotonin here, different axes needed
for (i in 1:(length(analytesold)) ) {
  analyte<-analytesold[i]
  if (analyte!="Serotonin") {
    # mydatatablep=olddatatable[Analyte==analyte]
    analyteplotsoldmedia[[i]]<-ggplot(selectedanalytes[Blank=="TRUE"&Analyte==analyte], aes(x=media,y=lowerlimitadjusted,fill=media)) +    #shape=factor(media)  #&rawplotkeep=="yes" <-use this to exclude points below lod
      #   geom_hline(yintercept = as.numeric(oldlimitssel[1,i]),colour="gray")+
      geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+
      geom_point(position=position_jitterdodge(jitter.width =.7),size=3,pch=21)+
      #   facet_wrap(~media,strip.position="bottom",scales = "free_x",labeller = labeller(media = media.labs)) +
      scale_y_continuous(expand=c(0.2,3))+
      ylab("pg/ml") +
      labs(title = analyte)+
      scale_x_discrete(labels=c("Sucrose", "HFCS","Trehalose"))+
      scale_fill_manual(name = "Media",values=sugarcolours) +
      theme_classic() +
      theme(
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.title.x= element_blank(), #element_text(size=1, colour="white"),
        axis.title.y= element_text(size=14),
        axis.text.x= element_text(size=12, colour = "black"),
        axis.text.y= element_text(size=12, colour = "black"), 
        #axis.text.x=element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        plot.title = element_text(size=14,face="bold")
      )  
  }
}

pdf(file="mediaalonerawvaluesv2.pdf",width=9,heigh=3)
grid.arrange(analyteplotsoldmedia[[2]],analyteplotsoldmedia[[6]],analyteplotsoldmedia[[7]],ncol=3)
dev.off()


#Figure 5D-F
analytemediaplots<- list()
for (i in 1:(length(analytes)) ) {
  analyte<-analytes[i]
  analytemediaplots[[i]]<-ggplot(hormonedata[Analyte==analyte&treatment=="nothing"], aes(x=media,y=Concentration,fill=media)) + #Concentration_fc_removal
  geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+
  geom_point(position=position_jitterdodge(jitter.width =.7),size=3,pch=21)+
  ylab("pg/ml") +  #change to ng/ml for serotonin
  labs(title = analyte)+
  scale_y_continuous(expand=c(0.2,3))+
  scale_fill_manual(name = "Media",values=sugarcolours) +
  scale_x_discrete(labels=c("Sucrose", "HFCS","Trehalose"))+
  theme_classic() +
  theme(
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=14),
    axis.text.x= element_text(size=12, colour = "black"),
    axis.text.y= element_text(size=12, colour = "black"), 
    #axis.text.x=element_blank(),
    legend.position="none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size=14,face="bold")
  )  
}


analyte="Ghrelin"
ghrelin<-ggplot(hormonedata[Analyte==analyte&treatment!="nothing"], aes(x=ReactorReplicate,y=Concentration_fc,fill=media)) + #Concentration_fc_removal
  geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+
  #geom_dotplot(binaxis = "y",stackdir = "center",colour="black",dotsize=0.5)+
  geom_point(position=position_jitterdodge(jitter.width =.7),size=3,pch=21)+
  #geom_hline(data=oldlimits, yintercept = oldlimits[1,i+1])+
  facet_wrap(~media,strip.position="bottom",scales = "free_x",labeller = labeller(media = media.labs)) +
  ylab("Fold change") +
  labs(title = analyte)+
  ylim(c(0.7,1.3))+  #use c(0.3,1.3) for PYY; #ghrelin = 0.7,1.3; PP = c(0.1,2.3)
  scale_fill_manual(name = "Media",values=sugarcolours) +
  theme_classic() +
  theme(
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=14),
    axis.text.x= element_text(size=12, colour = "black"),
    axis.text.y= element_text(size=12, colour = "black"), 
    #axis.text.x=element_blank(),
    legend.position="none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size=14,face="bold")
  )  

pdf(file="HormoneFCplots.pdf",width=12,height=3.5)
grid.arrange(ghrelin,PP,PYY,ncol=3)
dev.off()

kruskal.test(Concentration_fc_removal~as.factor(media), data=selectedanalytes[Analyte=="Ghrelin"&treatment!="nothing"])
DunnTest(Concentration_fc_removal~as.factor(media), data=selectedanalytes[Analyte=="Ghrelin"&treatment!="nothing"], method="BH")
# mean.rank.diff    pval    
# SDRBRMHFCS-SDRBRMSucrose           13.744444 0.00081 ***
#   SDRBRMTrehalose-SDRBRMSucrose       1.188889 0.75277    
# SDRBRMTrehalose-SDRBRMHFCS        -12.555556 0.00178 **

kruskal.test(Concentration_fc_removal~as.factor(media), data=selectedanalytes[Analyte=="PP"&treatment!="nothing"])
DunnTest(Concentration_fc_removal~as.factor(media), data=selectedanalytes[Analyte=="PP"&treatment!="nothing"], method="BH")
# mean.rank.diff    pval    
# SDRBRMHFCS-SDRBRMSucrose          -13.833333  0.0019 ** 
#   SDRBRMTrehalose-SDRBRMSucrose       4.583333  0.2858    
# SDRBRMTrehalose-SDRBRMHFCS         18.416667 5.4e-05 ***

kruskal.test(Concentration_fc_removal~as.factor(media), data=selectedanalytes[Analyte=="PYY"&treatment!="nothing"])
DunnTest(Concentration_fc_removal~as.factor(media), data=selectedanalytes[Analyte=="PYY"&treatment!="nothing"], method="BH")
# mean.rank.diff    pval    
# SDRBRMHFCS-SDRBRMSucrose         -16.6666667 0.00032 ***
#   SDRBRMTrehalose-SDRBRMSucrose     -0.3333333 0.93818    
# SDRBRMTrehalose-SDRBRMHFCS        16.3333333 0.00032 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1




