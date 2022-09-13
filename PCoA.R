
library("ggplot2")
library("Biostrings")
library("vegan")
library("phyloseq")
library("gridExtra")
library("RColorBrewer")
library("reshape2")
library("onewaytests")
library("data.table")

#Supplemental Fig 1
mergeddata<-read.delim(file="SequencingData.txt",header=TRUE, sep="\t")
wide<-dcast(mergeddata,experiment+sample+time+generation+Sugar~strain, value.var="correctedfraction2")
wide["sampletime"]<-paste(wide$sample,wide$time,sep="_")
widesamples<-wide

otumatrix<-data.matrix(cbind(widesamples[,6:11]*100))
rownames(otumatrix)<-widesamples$sampletime
otumatrix<-otumatrix[complete.cases(otumatrix), ]
sampletable<-otu_table(otumatrix,taxa_are_rows=FALSE)

#make metadata:
#sampletime, experiment, sample, time, generation, sugar (add later)
widemetadata<-cbind.data.frame(widesamples$sampletime,widesamples$experiment,widesamples$sample,widesamples$time, widesamples$generation,widesamples$Sugar)
names(widemetadata)<-c("sampletime","experiment","sample","time","generation","Sugar")
widemetadata$sugartime<-paste(widemetadata$Sugar,widemetadata$time,sep=".")
metadata<-sample_data(widemetadata)
sample_names(metadata)<-widemetadata$sampletime
#metadata$time<-as.numeric(as.character(gsub("takedown", "21", metadata$time)))  #to make it easier to work with the takedown for plotting
phyloseq<-merge_phyloseq(sampletable,metadata)
bray_pcoa = ordinate(phyloseq, "PCoA", "bray")
ordplot<-plot_ordination(phyloseq ,bray_pcoa, justDF=TRUE )
ordtable<-as.data.table(ordplot)

#make plot of each timepoint separately over sugars
sugarcolours<-c("#6A5ACD","#FFA500","violetred1")
names(sugarcolours)<-c("Sucrose","HFCS","Trehalose")
ordtable$time<-factor(ordtable$time,c("0","3","6","7","9","11","13","15","18","takedown"))

days<-levels(as.factor(ordtable$time))
i<-6
ordplotsbytime <- list()
for (i in 1:(length(days))) {
  day<-days[i]
  if (day =="0" | day =="9" | day=="15"){
  ordplotsbytime[[i]]<-ggplot(ordtable[time==day], aes(x=Axis.1,y=Axis.2,colour=Sugar))+
    geom_point(data=ordtable[time!=day],aes(x=Axis.1,y=Axis.2),colour="gray",alpha=0.2)+
    geom_point(size= 1) +
    stat_ellipse()+
  scale_colour_manual(name = "SugarDay",values=sugarcolours) +
 #   scale_fill_manual(name = "SugarDay",values=sugarcolours) +
   #  ylim(c(-0.3,0.7))+
    #xlim(c(-0.6,0.35))+
    labs(title = paste("Day ",day) )+
    theme_classic() + theme(
      axis.line.x = element_line(colour = "black", size = 1), 
      axis.line.y = element_line(colour = "black", size = 1), legend.position="none"
    )
    #
  }
  else {
    ordplotsbytime[[i]]<-ggplot(ordtable[time==day], aes(x=Axis.1,y=Axis.2,colour=Sugar))+
      geom_point(data=ordtable[time!=day],aes(x=Axis.1,y=Axis.2),colour="gray",alpha=0.2)+
      geom_point(size= 1) +
      scale_colour_manual(name = "SugarDay",values=sugarcolours) +
      #   scale_fill_manual(name = "SugarDay",values=sugarcolours) +
      #  ylim(c(-0.3,0.7))+
      #xlim(c(-0.6,0.35))+
      labs(title = paste("Day ",day) )+
      theme_classic() + theme(
        axis.line.x = element_line(colour = "black", size = 1), 
        axis.line.y = element_line(colour = "black", size = 1), legend.position="none"
      )
  }
}

pdf(file="BrayCurtisbyTime.pdf",width=12,height=8)
grid.arrange(ordplotsbytime[[1]],ordplotsbytime[[2]],ordplotsbytime[[3]],ordplotsbytime[[4]],
             ordplotsbytime[[5]],ordplotsbytime[[6]],ordplotsbytime[[7]],ordplotsbytime[[8]],
             ordplotsbytime[[9]],ordplotsbytime[[10]],ncol=5)
dev.off()

#try adonis on each time
phyloseqnoi<-subset_samples(phyloseq, Sugar!="inoculate")
dm<-phyloseq::distance(phyloseqnoi,method = "bray", type="samples")
df<-data.frame(sample_data(phyloseqnoi))
dft<-as.data.table(df)
timept<-adonis(dm ~ Sugar + time, data = dft,permutations = 9999) #everything but the inoculates

# adonis(formula = dm ~ Sugar + time, data = dft, permutations = 9999) 
# 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Sugar       2    0.5941 0.29706  4.1309 0.03931 0.0015 ** 
#   time        9    6.8265 0.75851 10.5479 0.45164 0.0001 ***
#   Residuals 107    7.6945 0.07191         0.50906           
# Total     118   15.1151                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# adonis(formula = dm ~ Sugar, data = dft, permutations = 9999) 
# 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Sugar       2    0.5941 0.29706   2.373 0.03931 0.0363 *
#   Residuals 116   14.5210 0.12518         0.96069         
# Total     118   15.1151                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

day<-15
pvalues<-c()
for (i in 1:(length(days))) {
  day<-days[i]
dm<-phyloseq::distance(subset_samples(phyloseq,time==day),method = "bray", type="samples")
df<-data.frame(sample_data(subset_samples(phyloseq,time==day)))
dft<-as.data.table(df)
timept<-adonis(dm ~ Sugar, data = dft,permutations = 9999)
pval<-timept$aov.tab$`Pr(>F)`[1]
pvalues<-c(pvalues,pval)
}
names(pvalues)<-days

# 0        3        6        7        9       11       13       15       18 takedown 
# 0.0107   0.3233   0.6343   0.1932   0.0217   0.1004   0.1642   0.0064   0.4763   0.3014 

p.adjust(pvalues,method="BH")
# 0          3          6          7          9         11         13         15         18   takedown 
# 0.05350000 0.40412500 0.63430000 0.32200000 0.07233333 0.25100000 0.32200000 0.05350000 0.52922222 0.40412500 

#timept 0
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Sugar      2   0.37098 0.185489  13.914 0.77671 0.0107 *
#   Residuals  8   0.10665 0.013331         0.22329         
# Total     10   0.47763                  1.00000         

#timept 9
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Sugar      2   0.24398 0.121989  6.6489 0.59637 0.0228 *
#   Residuals  9   0.16513 0.018347         0.40363         
# Total     11   0.40911                  1.00000 

#timept 15
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
# Sugar      2   0.65027 0.32514  5.1898 0.5356 0.0055 **
#   Residuals  9   0.56384 0.06265         0.4644          
# Total     11   1.21411                 1.0000  


