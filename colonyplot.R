require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)



data<-read.delim(file="Colonycounts.txt",sep="\t",header = TRUE, check.names=FALSE)

data$Sugar<-factor(data$Sugar,c("Sucrose","HFCS","Trehalose"))

plot<-ggplot(subset(data), aes(x=Generation, y=log10(ColonyCount),fill=as.factor(Sugar),shape=as.factor(Reactor))) +  #colour=Diet,fill=Experiment, shape=Rep
  geom_line()+  #
  geom_point(size=2,show.legend = FALSE)+
  facet_wrap(~Sugar,ncol=3)+
  scale_fill_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  scale_shape_manual(values = c("R1" = 21, "R2" = 22, "R3"=23, "R4" = 24))+
  xlab("Day") +
  ylab(bquote('Colony Counts (log'['10']*')')) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=9),
    axis.title.x= element_text(size=9),
    axis.text.x= element_text(size=8, colour = "black"),
    axis.text.y= element_text(size=8, colour = "black"),
    panel.border = element_blank(),
    strip.text = element_text(size=8, margin = margin()),
    strip.background = element_blank(),
    plot.title=element_blank())

pdf(file="ColonyPlot.pdf",width=4.5,height=2)
plot
dev.off()
