#to plot the sugar specific growth curves

library("dplyr")
library("ggplot2")
library("reshape2")
library("data.table")
library("gridExtra")
library("stringr")

#Figure 4C and Supplemental Figure 3
allplates1<-read.delim(file="regrowplate1.txt",header=TRUE,sep="\t")
allplates2<-read.delim(file="regrowplate2.txt",header=TRUE,sep="\t")
allplates1$run<-1
allplates2$run<-2
allplates<-rbind(allplates1,allplates2)

platesmelt<-melt(allplates,id.var=c("plate","well","strip", "reactor","day","sugar","nativesugar","run"),variable.name = "time",value.name = "OD")
platesmelt$time<-as.numeric(gsub('X','',platesmelt$time))
#plot blanks
blanks<-subset(platesmelt,platesmelt$strip=="blank")
#get averages
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

averageblanks<-summarySE(blanks,measurevar = "OD",groupvars = c("plate","strip","reactor","day","sugar","nativesugar","time","run"))
averageblanks$N
min(averageblanks$sd)
max(averageblanks$sd)
#substract blanks
#merge blank with number

blankinfo<-data.frame(averageblanks$run,averageblanks$plate,averageblanks$sugar,averageblanks$time,averageblanks$OD)
names(blankinfo)<-c("run","plate","sugar","time","OD")
#wells with data
datawells<-subset(platesmelt,platesmelt$strip !="")
blanked<-merge(x=datawells, y=blankinfo,by=c("sugar","time","plate"),all.x=TRUE)
blanked$ODblanked<-blanked$OD.x-blanked$OD.y

#get averages and sd
averageblanked<-summarySE(blanked,measurevar = "ODblanked",groupvars = c("plate","strip","reactor","day","sugar","nativesugar","time"))
min(averageblanked$N)  

averageblanked$name1=paste(averageblanked$strip,averageblanked$reactor,averageblanked$day,sep="_")


averageblankedtable<-as.data.table(averageblanked)

reactorsdays<-levels(as.factor(averageblanked$name1))
averageblankedtable$daynum=as.numeric(gsub('D','',averageblankedtable$day))
averageblankedtable$name2=paste(averageblanked$strip,averageblanked$reactor,sep="")
#linetypes<-c("dotted","dashed","solid")

reactors<-levels(factor(averageblankedtable$name2))
reactorsplots <- list()
nativesugarreactors<-c("blank",rep("Trehalose",4),rep("Sucrose",4),rep("HFCS",4))
for(j in 2:length(reactors) ) {
  reactorsplots[[j]]<-ggplot(averageblankedtable[name2==reactors[j]], aes(x=time/60,y=ODblanked,shape=as.factor(daynum))) +  
  #geom_errorbar(aes(ymin=ODblanked-sd, ymax=ODblanked+sd), width=.1,alpha=0.5) +
  geom_ribbon(aes(ymin=ODblanked-sd, ymax=ODblanked+sd,fill=sugar),alpha=0.1)+
  geom_line(aes(colour=sugar,linetype=as.factor(daynum)),size=1,) +
  #geom_point(aes(colour=sugar),size=1) +
 #geom_text(data=averageblankedtable[name2==reactors[j]& time==max(averageblankedtable$time)], aes(x= time/60+.3,label=daynum),size=6)+  #toggle to get annotations
  #geom_text(averageblankedtable[name2==reactors[j]&&time==max(averageblankedtable$time)],aes(x=time/60,y=ODblanked,label=daynum) )+
  xlab("Hours") +
  ylab(bquote(paste('OD'['600']))) +
  xlim(0,19)+
  labs(title = paste(nativesugarreactors[j],substring(reactors[j],3,4),sep=" : "))+
  scale_colour_manual(values = c("trehalose" = "violetred1", "sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  scale_fill_manual(values = c("trehalose" = "violetred1", "sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  scale_linetype_manual(values=c("6" = "dotted","12" = "dashed","18" = "solid"))+
  theme_classic() +
    theme(
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.title.x= element_text(size=16,colour="black"),
      axis.title.y= element_text(size=16,colour="black"),
      axis.text.x= element_text(size=12,colour="black"),
      axis.text.y= element_text(size=12,colour="black"),
      legend.position="none"
  )
}
pdf(file="growthcurvesovergenerations.pdf",width=12,height=10)
grid.arrange(reactorsplots[[6]],reactorsplots[[7]],reactorsplots[[8]],reactorsplots[[9]],
reactorsplots[[10]],reactorsplots[[11]],reactorsplots[[12]],reactorsplots[[13]],
             reactorsplots[[2]],reactorsplots[[3]],reactorsplots[[4]],reactorsplots[[5]],ncol=4)

dev.off()

pdf(file="S2R2S1R3a.pdf",width=12,height=10)
grid.arrange(reactorsplots[[7]],reactorsplots[[4]],ncol=2)

dev.off()
