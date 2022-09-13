library("ggplot2")
library("RColorBrewer")
library("dplyr") #for filtering dataframes
library("Polychrome")
library("gridExtra")
library("data.table")
library("ggsci")
library("reshape2")
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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

reactorinfo<-read.delim(file="reactorinfo.txt",sep='\t',header=TRUE)
seqdata<-read.delim(file="Rawsequencingdataforprocessing.txt",sep="\t",header=TRUE)
seqdatamelt<-reshape2::melt(seqdata,id.vars=names(seqdata)[1],variable.name="strain",value.name="fraction")
mergeddata<-merge(reactorinfo,seqdatamelt,by=c("samplename","strain"),all.x=TRUE)
mergeddata$correctedfraction<-mergeddata$fraction*4.7/mergeddata$CN #where 4.7 is the median/mean of the 16S CN for all 6 strains
mergeddata$sum<-ave(mergeddata$fraction,mergeddata$experiment,mergeddata$time, mergeddata$generation,mergeddata$sample, FUN=sum)
mergeddata$correctedsum<-ave(mergeddata$correctedfraction,mergeddata$experiment,mergeddata$time, mergeddata$generation,mergeddata$sample, FUN=sum)
mergeddata$correctedfraction2<-mergeddata$correctedfraction*1/mergeddata$correctedsum
mergeddata$correctedsum2<-ave(mergeddata$correctedfraction2,mergeddata$experiment,mergeddata$time, mergeddata$generation,mergeddata$sample, FUN=sum)

mergedtab<-as.data.table(mergeddata)
results<-mergedtab[time!="inoculate" & samplename!="S3R4_0"]  #remove inoculates and S3R4_0 which didn't pass rarefaction
results<-droplevels(results)
inoculates<-mergedtab[time=="inoculate"]
results$sample<-factor(results$sample,c("S2R1", "S2R2","S2R3","S2R4","S3R1","S3R2","S3R3","S3R4","S1R1","S1R2","S1R3","S1R4"))


#SuppFig2 of all reactors; if want to see the takedown sample, remove the !=takedown
#plots of all reactors individually, excluding the takedown; if want to include takedown, change the conditional
shapevalues<-rep(c(1,0,5,2),3)
names(shapevalues)<-levels(factor(results$sample))

reactors<-levels(as.factor(results$sample))
plots <- list()
sugar<-c("Sucrose","Sucrose","Sucrose","Sucrose","HFCS","HFCS","HFCS", "HFCS","Trehalose","Trehalose","Trehalose","Trehalose")
for (i in 1:length(reactors) ){
  datasub<-results[sample==reactors[i]]
  reactorname<-str_trunc(reactors[i], 2,"left",ellipsis ="")
  if(datasub$Sugar== "HFCS"| datasub$Sugar=="Trehalose") {
  plots[[i]]<-ggplot(datasub[time!="takedown"], aes(x=as.numeric(generation), y=correctedfraction2, color=strain, shape=sample)) +  
    geom_rect(aes(xmin=0, xmax=57, ymin=0, ymax=0.94), fill="light gray",color=NA, alpha=0.5)+
    geom_line(size=1) +
    geom_point(size=3.0) +
    xlab("Generations") +
    ylab("Fraction") +
    ylim(0,1)+
     scale_color_d3()+
    scale_shape_manual(name = "reactor",values=shapevalues)+
    labs(title = paste(sugar[i],":",reactorname))+
    theme_classic() +
    theme(
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.title.x= element_text(size=15, colour = "black"),
      axis.title.y= element_text(size=15, colour="black"),
      axis.text.x= element_text(size=12, colour="black"),
      axis.text.y= element_text(size=12,colour="black"), 
      axis.ticks.y = element_line(colour = "black"),
      axis.ticks.x = element_line(colour = "black"),
      legend.position="none"
    )
  }
  else {
    plots[[i]]<-ggplot(datasub[time!="takedown"], aes(x=as.numeric(generation), y=correctedfraction2, color=strain, shape=sample)) +  
    #  geom_rect(aes(xmin=0, xmax=54, ymin=0, ymax=0.9), fill="gray",color=NA, alpha=0.5)+
      geom_line(size=1) +
      geom_point(size=3.0) +
      xlab("Generations") +
      ylab("Fraction") +
      ylim(0,1)+
      # geom_vline(xintercept=57,colour="black",lwd=0.8)+
      # geom_vline(xintercept=17,colour="black",lwd=1,linetype="dotted")+
      #scale_colour_manual(name = "strain",values=colours)+
      scale_color_d3()+
      scale_shape_manual(name = "reactor",values=shapevalues)+
      #scale_linetype_manual(name="flow",values=flowvalues)+
      labs(title = paste(sugar[i],":",reactorname))+
      theme_classic() +
      theme(
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.title.x= element_text(size=15, colour = "black"),
        axis.title.y= element_text(size=15, colour="black"),
        axis.text.x= element_text(size=12, colour="black"),
        axis.text.y= element_text(size=12,colour="black"), 
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        legend.position="none"
      )
  }
}

pdf(file=paste("AllreactorsCNcorrectedcomposition.pdf"),width=12,height=8)
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
             plots[[5]],plots[[6]],plots[[7]],plots[[8]],
             plots[[9]],plots[[10]],plots[[11]],plots[[12]],
             ncol=4)
dev.off()

#Fig 2A
averageresults<-summarySE(results,measurevar = "correctedfraction2",groupvars = c("strain","time"))
avefirst50<-ggplot(subset(averageresults,averageresults$time!="takedown" & as.numeric(averageresults$time)<7),aes(x=as.numeric(time),y=correctedfraction2,color=strain))+
  geom_line(size=1) +
  geom_errorbar(aes(ymin=correctedfraction2-sd,ymax=correctedfraction2+sd),width=.5,alpha=0.5,size=0.8)+ 
  geom_point(size=3.0) +
  xlab("Day") +
  ylab("Fraction") +
   scale_color_d3()+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=15, colour = "black"),
    axis.title.y= element_text(size=15, colour="black"),
    axis.text.x= element_text(size=12, colour="black"),
    axis.text.y= element_text(size=12,colour="black"), 
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    legend.position="none"
  )

pdf(file="first50genaverage.pdf",width=3, height=2.5)
avefirst50
dev.off()

averageresultssugar<-summarySE(results,measurevar = "correctedfraction2",groupvars = c("strain", "Sugar","time"))
averagegens<-summarySE(results,measurevar = "generation",groupvars = c("strain", "Sugar","time"))
mergedave<-merge(averageresultssugar,averagegens,by=c("strain", "Sugar","time"))

mergedave$Sugar<-factor(mergedave$Sugar, c("Sucrose","HFCS","Trehalose"))


#Fig2B
pdf(file="Generation50barplot.pdf",width=2,height=2)
ggplot(results[time==6], aes(x=sample, y=correctedfraction2, fill=strain)) +  
  geom_col()+
  geom_vline(xintercept=4.5,  colour = "black")+ #linetype
  geom_vline(xintercept=8.5, colour = "black")+
  # xlab("Days post inoculation") +
  ylab("Fraction") + xlab("Sugar")+
  ggtitle(paste("Generation ","50",sep=""))+
  scale_color_d3(palette = "category20")+
  scale_fill_d3(palette = "category20")+
  geom_text(aes(x=2.5, y=-.05, label="Sucrose",size=50))+
  geom_text(aes(x=6.5, y=-.05, label="HFCS",size=50))+
  geom_text(aes(x=10.5, y=-.05, label="Trehalose",size=50))+
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=18,colour="black"),
    axis.title.y= element_text(size=18,colour="black"),
    # axis.text.x= element_text(size=12,colour="black"),
    axis.text.x= element_blank(),
    axis.text.y= element_text(size=12,colour="black"),
    axis.ticks.y = element_line(colour="black"),  
    axis.ticks.x = element_blank(), 
    
    legend.position="none"
  )
dev.off()


#Fig 3A
aveallbysugar<-ggplot(subset(mergedave,mergedave$time!="takedown"),aes(x=as.numeric(generation),y=correctedfraction2,color=strain))+
  geom_line(size=1) +
  geom_errorbar(aes(ymin=correctedfraction2-sd.x,ymax=correctedfraction2+sd.x),width=.5,alpha=0.5,size=0.8)+ 
  geom_point(size=3.0) +
  facet_wrap(~Sugar)+
  xlab("Generation") +
  ylab("Fraction") +
  #  ylim(0,1)+
  # geom_vline(xintercept=57,colour="black",lwd=0.8)+
  # geom_vline(xintercept=17,colour="black",lwd=1,linetype="dotted")+
  scale_color_d3()+
  #labs(title = paste(sugar[i],":",reactorname))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=15, colour = "black"),
    axis.title.y= element_text(size=15, colour="black"),
    axis.text.x= element_text(size=12, colour="black"),
    axis.text.y= element_text(size=12,colour="black"), 
    axis.ticks.y = element_line(colour = "black", ),
    axis.ticks.x = element_line(colour = "black"),
    legend.position="none"
  )


pdf(file="avebysugar.pdf",width=7, height=3)
aveallbysugar
dev.off()

#what strains increased/decreased
results$sugarchangestatus
for (i in 1:length(results$generation) ){
  gen<-results$generation[i]
 if (gen<=57) {
   results$sugarchangestatus[i] ="before"
 } 
  else if (gen>57) {
    results$sugarchangestatus[i] ="after"
  } 
}

library("DescTools")

strainlist<-levels(as.factor(results$strain))
pvaluesafter<-c()
for(i in 1:length(strainlist)){
  mystrain<-strainlist[i]
modeal<-kruskal.test(correctedfraction2~as.factor(Sugar), data=results[strain==mystrain & generation>57])
pval<-modeal$p.value
pvaluesafter<-c(pvaluesafter,pval)
}
names(pvaluesafter)<-strainlist
pafteradjust<-p.adjust(pvaluesafter,method="BH")
# Lcas         Lreu         Phis         Smit         Ssal         Vaty 
# 0.0766489955 0.0001395372 0.3590513330 0.0002840456 0.0046605700 0.0047421836

pvaluesbefore<-c()
for(i in 1:length(strainlist)){
  mystrain<-strainlist[i]
  modeal<-kruskal.test(correctedfraction2~as.factor(Sugar), data=results[strain==mystrain & generation<57])
  pval<-modeal$p.value
  pvaluesbefore<-c(pvaluesbefore,pval)
}
names(pvaluesbefore)<-strainlist
pbeforeadjust<-p.adjust(pvaluesbefore,method="BH")
# Lcas        Lreu        Phis        Smit        Ssal        Vaty 
# 0.886835766 0.006388234 0.897794863 0.006388234 0.366053084 0.884576072 

DunnettTest(correctedfraction2~as.factor(Sugar), data=results[strain=="Vaty" & generation>57], control="Sucrose")

# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# 
# $Sucrose
#                         diff      lwr.ci    upr.ci    pval    
# HFCS-Sucrose      0.06640258 -0.07197703 0.2047822 0.45404    
# Trehalose-Sucrose 0.23387949  0.09549988 0.3722591 0.00054 ***
# 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

DunnettTest(fraction~as.factor(Sugar), data=results[strain=="Vaty" & generation<=57], control="Sucrose")

# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# 
# $Sucrose
#                           diff      lwr.ci    upr.ci   pval    
# HFCS-Sucrose      0.0004018636 -0.13046527 0.1312690 1.0000    
# Trehalose-Sucrose 0.0536330985 -0.07435748 0.1816237 0.5301    
# 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#just for Vatyp before and after the sugar change

DunnettTest(correctedfraction2~as.factor(Sugar), data=results[strain=="Ssal" & generation>57], control="Sucrose")
# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# 
# $Sucrose
#                          diff     lwr.ci      upr.ci   pval    
# HFCS-Sucrose      -0.03982661 -0.1027412  0.02308802 0.2681    
# Trehalose-Sucrose -0.07660580 -0.1395204 -0.01369118 0.0143 *  
# 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

DunnettTest(correctedfraction2~as.factor(Sugar), data=results[strain=="Ssal" & generation<=57], control="Sucrose")
# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# 
# $Sucrose
#                          diff     lwr.ci    upr.ci   pval    
# HFCS-Sucrose      -0.04305592 -0.1880500 0.1019382 0.7186    
# Trehalose-Sucrose  0.01386113 -0.1279459 0.1556682 0.9635    
# 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

DunnettTest(correctedfraction2~as.factor(Sugar), data=results[strain=="Smit" & generation>57], control="Sucrose")
# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# 
# $Sucrose
#                           diff        lwr.ci      upr.ci   pval    
# HFCS-Sucrose      0.0004265656 -0.0011254644 0.001978596 0.7613    
# Trehalose-Sucrose 0.0021228797  0.0005708497 0.003674910 0.0055 ** 
# 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

DunnettTest(correctedfraction2~as.factor(Sugar), data=results[strain=="Smit" & generation<=57], control="Sucrose")
# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# 
# $Sucrose
#                           diff       lwr.ci     upr.ci   pval    
# HFCS-Sucrose      0.0008761258 -0.030725558 0.03247781 0.9970    
# Trehalose-Sucrose 0.0291271579 -0.001779898 0.06003421 0.0667 .  
# 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

results$sugarchangestatus<-factor(results$sugarchangestatus,c("before", "after"))
results$Sugar<-factor(results$S,c("Sucrose", "HFCS", "Trehalose"))
shapvaluesbeforeafter<-c(21,22)
names(shapvaluesbeforeafter)<-c("before","after")

#Fig 3B and C
sugarchangestrainsstatus<-ggplot(results[strain=="Vaty" | strain == "Ssal" ], aes(x=Sugar,y=fraction,fill=Sugar,shape=sugarchangestatus)) +
  geom_boxplot(size=.8,outlier.shape = NA,alpha = 0.4,colour="black")+
  geom_point(position=position_jitterdodge(jitter.width =.8),size=4)+
  # scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  scale_fill_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  facet_wrap(~strain + sugarchangestatus,ncol=2,scales="free")+
  ylab("Fraction")+
  # labs(x="Sugar",y=expression(paste("Fraction of ",italic("V. atypica"))))+
  scale_shape_manual(name = "sugarchangestatus",values=shapvaluesbeforeafter)+
 # ylim(c(0,1))+
  scale_y_continuous(expand=expansion(mult= c(0.07,.35)))+
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black", size=1),
    # axis.line = element_line(colour = "black". size=1.2),
    axis.title.x= element_text(size=16,colour="black"),
    axis.title.y= element_text(size=16,colour="black"),
    axis.text.y= element_text(size=12,colour="black"),
    axis.ticks = element_line(colour="black"),  
    panel.border = element_blank(),
    axis.text.x= element_text(size=12, colour = "black"), #angle=45,vjust=0.6
    strip.text = element_text(size=16, margin = margin()),
    strip.background = element_blank(),
    plot.title=element_blank(),
    legend.position="none")

pdf(file="BeforeandAfterSugarChange.pdf",width=4.5,height=5)
sugarchangestrainsstatus
dev.off()


#Fig 3D
pdf(file="GenerationTakedownbarplot.pdf",width=2,height=2)
ggplot(results[time=="takedown"], aes(x=sample, y=correctedfraction2, fill=strain)) +  
  geom_col()+
  geom_vline(xintercept=4.5,  colour = "black")+ #linetype
  geom_vline(xintercept=8.5, colour = "black")+
  # xlab("Days post inoculation") +
  ylab("Fraction") + xlab("Sugar")+
  ggtitle(paste("Generation ","50",sep=""))+
  scale_color_d3(palette = "category20")+
  scale_fill_d3(palette = "category20")+
  geom_text(aes(x=2.5, y=-.05, label="Sucrose",size=50))+
  geom_text(aes(x=6.5, y=-.05, label="HFCS",size=50))+
  geom_text(aes(x=10.5, y=-.05, label="Trehalose",size=50))+
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=18,colour="black"),
    axis.title.y= element_text(size=18,colour="black"),
    # axis.text.x= element_text(size=12,colour="black"),
    axis.text.x= element_blank(),
    axis.text.y= element_text(size=12,colour="black"),
    axis.ticks.y = element_line(colour="black"),  
    axis.ticks.x = element_blank(), 
    
    legend.position="none"
  )
dev.off()


#Supplemental Figure 6B
#run the inoculates
#reorder the plotting of the inoculates
inoculationorder<-data.frame("samplename"=c("Phis","LactoStreps","Atyp"), "order"=c(1,2,3))
inoculatesorder<-merge(x=inoculates,y=inoculationorder,by="samplename",all.x=TRUE)

pdf(file="inoculates.pdf",width=3,height=4)
ggplot(inoculatesorder, aes(x=as.character(order), y=correctedfraction2, fill=strain)) +  
  geom_col()+
  xlab("Inoculate") +
  scale_x_discrete(labels=c("1"="Day -2", "2"="Day -1", "3" = "Day 2"))+
  ylab("Fraction") +
  scale_fill_d3()+
  theme_classic()+
  theme(
    axis.line = element_line(colour = "black", size=1),
    # axis.line = element_line(colour = "black". size=1.2),
    axis.title.x= element_text(size=16,colour="black"),
    axis.title.y= element_text(size=16,colour="black"),
    axis.text.y= element_text(size=12,colour="black"),
    axis.ticks = element_line(colour="black"),  
    panel.border = element_blank(),
    axis.text.x= element_text(size=12, colour = "black",angle=45,vjust=0.6))
  
dev.off()











