library("ggplot2")
library("dplyr") #for filtering dataframes
library("Polychrome")
library("gridExtra")
library("data.table")
library("ggdendro")
library("scales")
library("stringr")
library("reshape")

#Supplemental Figure 5
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


Strains <- read.delim("Straingrowth.txt", header = TRUE)

sample_summary <- summarySE(Strains, measurevar = "value", groupvars = c("timemins", "Sugar", "strain"))

Lcas<-subset(sample_summary, sample_summary$strain == "Lcas")
Vatp<-subset(sample_summary, sample_summary$strain == "Vatp")
Lreu<-subset(sample_summary, sample_summary$strain == "Lreu")
Smit<-subset(sample_summary, sample_summary$strain == "Smit")
Ssal<-subset(sample_summary, sample_summary$strain == "Ssal")
Pjej<-subset(sample_summary, sample_summary$strain == "Pjej")


#pdf("Lcas.pdf", width =7 , height =4)
u <-ggplot(Lcas, aes(x=timemins, y=value, colour=Sugar)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=4.5) +
  xlab("Minutes") +
  ylab("Growth Value") +
  scale_colour_manual(values = c("Trehalose" = "#00FFFF", "Sucrose" = "#6A5ACD", "NoSugar"= "black", "HFCS"="#191970"))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


v <- ggplot(Lreu,aes(x=timemins, y=value, colour=Sugar)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=4.5) +
  xlab("Minutes") +
  ylab("Growth Value") +
  scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "NoSugar"= "black", "HFCS"="#FFA500"))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


#pdf("Smit.pdf", width =7 , height =4)
w <- ggplot(Smit, aes(x=timemins, y=value, colour=Sugar)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=4.5) +
  xlab("Minutes") +
  ylab("Growth Value") +
  scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "NoSugar"= "black", "HFCS"="#FFA500"))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )

#pdf("Ssal.pdf", width =7 , height =4)
x <-ggplot(Ssal, aes(x=timemins, y=value, colour=Sugar)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=4.5) +
  xlab("Minutes") +
  ylab("Growth Value") +
  scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "NoSugar"= "black", "HFCS"="#FFA500"))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )
#pdf("Pjej.pdf", width =7 , height =4)
y <-  ggplot(Pjej, aes(x=timemins, y=value, colour=Sugar)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=4.5) +
  xlab("Minutes") +
  ylab("Growth Value") +
  scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "NoSugar"= "black", "HFCS"="#FFA500"))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )
#pdf("Vatp.pdf", width =7 , height =4)
z <-  ggplot(Vatp, aes(x=timemins, y=value, colour=Sugar)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=4.5) +
  xlab("Minutes") +
  ylab("Growth Value") +
  scale_colour_manual(values = c("Trehalose" = "violetred1", "Sucrose" = "#6A5ACD", "NoSugar"= "black", "HFCS"="#FFA500"))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


#plots:
u <- u + ylim(0,0.8)
v <- v + ylim(0,0.8)
w <- w + ylim(0,0.8)
x <- x + ylim(0,0.8)
y <- y + ylim(0,0.8)
z <- z + ylim(0,0.8)
