#to linear model the analytes and the strain composition
#try to do a model like analyte~(microbes)*sugar

library("ggplot2")
library("dplyr") #for filtering dataframes
library("Polychrome")
library("gridExtra")
library("data.table")
library("ggdendro")
library("scales")
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

#load in microbe data
mergeddata<-read.delim(file="SequencingData.txt",header=TRUE, sep="\t")
microbedata<-subset(mergeddata,mergeddata$time=="takedown")

#load in analyte data
analytedata<-read.delim(file="selectedanalyteshormonedata.txt",header=TRUE,sep="\t")
limits<-read.delim(file="analytelimites.txt",header=TRUE,sep="\t")

#check the analytefile for problems
#does blankgroup have 2 if exp 3?
analytedatatable<-as.data.table(analytedata)
which(analytedatatable$Exp==3)==which(str_sub(analytedatatable$Blankgroup,-1,-1)==2)
#HFCS
which(analytedatatable$Blankgroup =="HFCSmedia" | analytedatatable$Blankgroup =="HFCSmedia2") ==which(analytedatatable$media=="SDRBRMHFCS")
#Suc
which(analytedatatable$Blankgroup =="sucmedia" | analytedatatable$Blankgroup =="sucmedia2") ==which(analytedatatable$media=="SDRBRMSucrose")
#Tre
which(analytedatatable$Blankgroup =="tremedia" | analytedatatable$Blankgroup =="tremedia2") ==which(analytedatatable$media=="SDRBRMTrehalose")
#media == "SDRBRMHFCS"
unique(analytedatatable[media == "SDRBRMHFCS"]$SampleName)
unique(analytedatatable[media == "SDRBRMSucrose"]$SampleName)
unique(analytedatatable[media == "SDRBRMTrehalose"]$SampleName)
#treatment
unique(analytedatatable[treatment == "nothing"]$SampleName)
#SampleGroup
problems<-which(str_sub(analytedatatable$SampleName,1,-2)!=analytedatatable$SampleGroup)
str_sub(analytedatatable$SampleName[problems],1,-3) ==analytedatatable$SampleGroup[problems]
analytedatatable[Blank == TRUE] == analytedatatable[treatment == "nothing"]

names(microbedata)
levels(as.factor(microbedata$reactorname))
levels(as.factor(analytedata$treatment))
analytedatatable$sample<-substring(analytedatatable$treatment,1,4)

#add blank spots for microbedata
blankrows<-c(0,NA,"Lcas", "noth", NA,NA,NA,NA,NA,"blank",NA,0,0,0,0,0,0,0)
blankrows2<-c(0,NA,"Lreu", "noth", NA,NA,NA,NA,"blank",NA,0,0,0,0,0,0,0)
blankrows3<-c(0,NA,"Ssal", "noth", NA,NA,NA,NA,"blank",NA,0,0,0,0,0,0,0)
blankrows4<-c(0,NA,"Smit", "noth", NA,NA,NA,NA,"blank",NA,0,0,0,0,0,0,0)
blankrows5<-c(0,NA,"Phis", "noth", NA,NA,NA,NA,"blank",NA,0,0,0,0,0,0,0)
blankrows6<-c(0,NA,"Vaty", "noth", NA,NA,NA,NA,"blank",NA,0,0,0,0,0,0,0)

microdatablank<-rbind(microbedata,blankrows,blankrows2,blankrows3,blankrows4,blankrows5,blankrows6)

#merge
alldatamerge<-merge(analytedatatable,microdatablank,by="sample",allow.cartesian = TRUE)
alldatatable<-as.data.table(alldatamerge)
#need to add no strains for the media alone controls
alldatatable[sample=="noth"]$time<-'takedown'
alldatatable[treatment=="nothing"]$sample<-alldatatable[treatment=="nothing"]$Blankgroup
mylinear<-list()

analytes<-names(limits)[2:15]
analyteset<-c(3,11,12,13)
i<-12
myanalyte<-analytes[i]
#put loop here

mysub<-alldatatable[Analyte==myanalyte  & SampleName!="S2R1D18A"]  
#add zeroes to blanks
mysub$correctedfraction2<-as.numeric(mysub$correctedfraction2)

dcasted<-dcast(mysub, sample + Analyte + Concentration + Concentration_fc + media+ time+ replica + Plate~strain, value.var = "correctedfraction2")
dcasted$media<-as.factor(dcasted$media)
dcasted$media <- relevel(dcasted$media, ref = 2)  #set ref level as sucrose

#change and explore the modeling
full<-lm(Concentration~(Lcas+Vaty
                        )*media,data=dcasted )

reduced<-lm(Concentration~media+Lcas+Vaty,data=dcasted )
print(myanalyte)
print(summary(full))
print(summary(reduced))
anova(full,reduced,test="LRT")
#see separate text file for how the models care about the different 
#microbial abundances
emmeans(full,pairwise~Ssal*media,adjust="BH")



sugarcolours<-c("#6A5ACD","#FFA500","violetred1")
names(sugarcolours)<-c("SDRBRMSucrose","SDRBRMHFCS","SDRBRMTrehalose")

dcastmelt<-melt.data.table(dcasted,id.vars=names(dcasted)[c(1:8)], variable.name = "Strain",value.name="RAbundance")
dcastmeltsummary<-summarySE(dcastmelt,measurevar = c("Concentration"),groupvars = c("Analyte","Strain","sample"))
dcastmeltsummary$N

strainnames<-c('Lcas' = "L. casei", 'Ssal' = "S. salivarius", 'Vaty' = "V. atypica", 'Lreu' = "L. reuteri")

models<-read.delim(file="models.txt",header=TRUE,sep="\t")

ghrelinplot<-ggplot (dcastmelt[ Strain=="Lcas" |Strain=="Vaty" |Strain=="Ssal"], aes(x=RAbundance, y =Concentration, colour=media)) +
  facet_wrap(.~Strain,ncol=4,scales="free", labeller = as_labeller(strainnames))+
  geom_point()+
#  geom_smooth(method="lm")+
  geom_abline(data = subset(models,models$Hormone=="Ghrelin" & models$Sig=="y"),aes(intercept=Intercept, slope=Coeff,colour=media),lwd=1)+
  geom_abline(data = subset(models,models$Hormone=="Ghrelin"& models$Sig=="y"),aes(intercept=UpperI, slope=UpperC,colour=media),lwd=1)+
  geom_abline(data = subset(models,models$Hormone=="Ghrelin"& models$Sig=="y"),aes(intercept=LowerI, slope=LowerC,colour=media),lwd=1)+
  scale_colour_manual(name = "media",values=sugarcolours)+
  scale_fill_manual(name = "media",values=sugarcolours) +
  ylab("Concentration") +
  xlab("Microbe abundance") +
  labs(title = myanalyte)+
  theme_classic() +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=14),
    axis.title.y= element_text(size=14),
    axis.text.x= element_text(size=12),
    axis.text.y= element_text(size=12), 
    plot.title = element_text(size=14,face="bold"),
    strip.background = element_blank(),
    strip.text =element_text(size=16,face="italic"),
    legend.position="none"
  )

pdf(file="ghrelinmodel.pdf",width=12,height=4)
ghrelinplot
dev.off()



PPplot<-ggplot (dcastmelt[ Strain=="Ssal" ], aes(x=RAbundance, y =Concentration, colour=media)) +
  facet_wrap(.~Strain,ncol=4,scales="free", labeller = as_labeller(strainnames))+
  geom_point()+
#  geom_smooth(method="lm")+
  geom_abline(data = subset(models,models$Hormone=="PP" & models$Sig=="y"),aes(intercept=Intercept, slope=Coeff,colour=media),lwd=1)+
  geom_abline(data = subset(models,models$Hormone=="PP"& models$Sig=="y"),aes(intercept=UpperI, slope=UpperC,colour=media),lwd=1)+
  geom_abline(data = subset(models,models$Hormone=="PP"& models$Sig=="y"),aes(intercept=LowerI, slope=LowerC,colour=media),lwd=1)+
  scale_colour_manual(name = "media",values=sugarcolours)+
  scale_fill_manual(name = "media",values=sugarcolours) +
  ylab("Concentration") +
  xlab("Microbe abundance") +
  labs(title = myanalyte)+
  theme_classic() +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=14),
    axis.title.y= element_text(size=14),
    axis.text.x= element_text(size=12,colour="black"),
    axis.text.y= element_text(size=12,colour="black"), 
    plot.title = element_text(size=14,face="bold"),
    strip.background = element_blank(),
    strip.text =element_text(size=16,face="italic"),
    legend.position="none"
  )

pdf(file="PPmodel.pdf",width=4,height=4)
PPplot
dev.off()

PYYplot<-ggplot (dcastmelt[ Strain=="Lcas" |Strain=="Vaty" ], aes(x=RAbundance, y =Concentration, colour=media)) +
  facet_wrap(.~Strain,ncol=4,scales="free", labeller = as_labeller(strainnames))+
  geom_point()+
  geom_abline(data = subset(models,models$Hormone=="PYY" & models$Sig=="y"),aes(intercept=Intercept, slope=Coeff,colour=media),lwd=1)+
  geom_abline(data = subset(models,models$Hormone=="PYY"& models$Sig=="y"),aes(intercept=UpperI, slope=UpperC,colour=media),lwd=1)+
  geom_abline(data = subset(models,models$Hormone=="PYY"& models$Sig=="y"),aes(intercept=LowerI, slope=LowerC,colour=media),lwd=1)+
  #geom_smooth(method="lm")+
  scale_colour_manual(name = "media",values=sugarcolours)+
  scale_fill_manual(name = "media",values=sugarcolours) +
  ylab("Concentration") +
  xlab("Microbe abundance") +
  labs(title = myanalyte)+
  theme_classic() +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=14),
    axis.title.y= element_text(size=14),
    axis.text.x= element_text(size=12,colour="black"),
    axis.text.y= element_text(size=12,colour="black"), 
    plot.title = element_text(size=14,face="bold"),
    strip.background = element_blank(),
    strip.text =element_text(size=16,face="italic"),
    legend.position="none"
  )

pdf(file="PYYmodel.pdf",width=8,height=4)
PYYplot
dev.off()



