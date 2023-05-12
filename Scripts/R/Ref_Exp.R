#Sets working directory if using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
setwd("..")

#Sets working directory for other R IDE's (don't run in RStudio)
setwd(getSrcDirectory(function(){})[1])
setwd("..")
setwd("..")

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)
library(ggplot2)
library(rstatix)
library(ggpubr)

# Functions
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

ref_data <- read.csv("Data\\Sco-GEM reference model RiPP results.csv")
ref_data$Core.Class[ref_data$Core.Class == "lassopeptide"] <- "Lasso Peptide" 
ref_data$Core.Class[ref_data$Core.Class == "thiopeptide"] <- "Thiopeptide" 
ref_data$Core.Class[ref_data$Core.Class == "lanthipeptide"] <- "Lanthipeptide" 

perf_data <- read.csv("RiPP Pathway Output\\RiPP reconstruction summary.csv")
#perf_data_mature <- perf_data[grepl("cyclic_peptide_c",perf_data$Metabolite.list,fixed=TRUE)|
                                #grepl("macrocyclic_thiopeptide_c",perf_data$Metabolite.list,fixed=TRUE)|
                                #grepl("macrolactam_peptide_c",perf_data$Metabolite.list,fixed=TRUE),]
ref_data$Core.Number <- perf_data$Core.Number

ref_data$Metabolite.list <- perf_data$Metabolite.list[match(ref_data$BGC.ID, perf_data$Name)&
                                                        match(ref_data$Core.Number,perf_data$Core.Number)]
ref_data_mature <- ref_data[grepl("cyclic_peptide_c",ref_data$Metabolite.list,fixed=TRUE)|
                            grepl("macrocyclic_thiopeptide_c",ref_data$Metabolite.list,fixed=TRUE)|
                            grepl("macrolactam_peptide_c",ref_data$Metabolite.list,fixed=TRUE),]


ref_data$Conc.SubClass <- paste(ref_data$Core.Class, ref_data$Core.Subclass)
ref_data <- ref_data[!ref_data$Core.Class=="sactipeptide",]
ref_data_mature$Conc.SubClass <- paste(ref_data_mature$Core.Class, ref_data_mature$Core.Subclass)
ref_data_mature <- ref_data_mature[!ref_data_mature$Core.Class=="sactipeptide",]

ref_data_high_g <- ref_data_mature[ref_data_mature$Gradient>1.0,]

nrow(ref_data_mature[ref_data_mature$Conc.SubClass=="Thiopeptide Type III",])

#Boxplot chart for gradient values distribution by subclass
p <- ref_data_mature %>% ggplot(aes(Core.Class, Gradient)) +
  geom_boxplot(varwidth = T, fill="plum") +
  labs(title="", 
       subtitle="",
       caption="",
       x="RiPP class",
       y="Growth-Production Gradient") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=30, hjust = 0.85, vjust = 0.85))
show(p)

ref_data_filt <- ref_data[!(ref_data$Conc.SubClass=="Thiopeptide Type I" | ref_data$Conc.SubClass=="Thiopeptide Type II"|ref_data$Core.Class=="sactipeptide"),]

ref_data_mature <- na.omit(ref_data_mature)
#Violin plot of gradient values distribution by subclass
p <- ref_data_mature[ref_data_mature$Core.Class=="Lasso Peptide",] %>% ggplot(aes(Conc.SubClass, Gradient)) +
  geom_violin() +
  labs(title="Lanthipeptide gradients", 
       subtitle="Violin plot of gradient values",
       caption="",
       x="RiPP Subclass",
       y="Gradient")
show(p)

#Pie chart of subclasses
p <- ref_data_mature %>% ggplot(aes(x="", fill=factor(Conc.SubClass))) +
  geom_bar(width=1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(title="RiPP subclasses", 
       subtitle="Pie chart for distribution of RiPP subclass for cores in dataset",
       caption="",
       x="",
       y="") +
  scale_fill_discrete(name="RiPP Subclass") +
  coord_polar(theta = "y", start=0)
show(p)

#Pie chart of classes
p <- ref_data_mature %>% ggplot(aes(x="", fill=factor(Core.Class))) +
  geom_bar(width=1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(title="RiPP classes", 
       subtitle="Pie chart for distribution of RiPP class for cores in dataset",
       caption="",
       x="",
       y="") +
  scale_fill_discrete(name="RiPP Class") +
  coord_polar(theta = "y", start=0)
show(p)

library(ggplot2)
library(webr)
library(dplyr)

p <- PieDonut(ref_data_mature, aes(Core.Class, Core.Subclass), showPieName = FALSE,showRatioThreshold = F, explode = 3, donutLabelSize = 4,
              pieLabelSize = 5)
show(p)

# Subclass data summary

ref_data_mature.summary <- data_summary(ref_data_mature2, varname="Gradient", groupnames=c("Core.Class","Core.Subclass"))
ref_data_mature.summary$Core.Subclass=as.factor(ref_data_mature.summary$Core.Subclass)

ref_data_mature2 <- ref_data_mature[ref_data_mature$Gradient<1.0,]

#Bar plot of gradients
p <- ref_data_mature.summary %>% ggplot(aes(x=Core.Class, y=Gradient, fill=Core.Subclass)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Gradient-sd, ymax=Gradient+sd), width=.2,
                position=position_dodge(.9)) +
  xlab("Core RiPP Class") +
  ylab("Growth-Production Gradient") +
  scale_fill_discrete(name="Core RiPP Subclass")
show(p)

# Filtered subclass data summary

ref_data_filt.summary <- data_summary(ref_data_filt, varname="Gradient", groupnames=c("Core.Class","Core.Subclass"))
ref_data_filt.summary$Core.Subclass=as.factor(ref_data_filt.summary$Core.Subclass)

#Bar plot of gradients (subclasses)
p <- ref_data_filt.summary %>% ggplot(aes(x=Core.Class, y=Gradient, fill=Core.Subclass)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Gradient-sd, ymax=Gradient+sd), width=.2,
                position=position_dodge(.9)) 
show(p)

# Class data summary
ref_data_mature.summary <- data_summary(ref_data_mature, varname="Gradient", groupnames=c("Core.Class"))

#Bar plot of gradients (classes)
p <- ref_data_mature.summary %>% ggplot(aes(x=Core.Class, y=Gradient, fill=Core.Class)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Gradient-sd, ymax=Gradient+sd), width=.2,
                position=position_dodge(.9)) +
  xlab("") +
  ylab("Growth-Production Gradient") +
  scale_fill_discrete(name="Core RiPP Class")
show(p)

#ANOVA and T.test for RiPP Class

ref_data_t.test <- ref_data %>% group_by(Core.Class) %>% get_summary_stats(Gradient, type = "mean_sd")
#ANOVA
res.aov <- ref_data %>% anova_test(Gradient ~ Core.Class)
res.aov
#Pairwise T.test
pwc <- pairwise.t.test(ref_data_mature2$Gradient, ref_data_mature2$Conc.SubClass, p.adjust.method="bonferroni")
pwc


table(ref_data_mature$Conc.SubClass)

#Peptide length metrics
mean(ref_data$Peplength)
sd(ref_data$Peplength)
median(ref_data$Peplength)

#Sets up methionine and glycine gradient data
metgly_data <- read.csv("Data\\Methionine-Glycine Data.csv")
metgly_long <- metgly_data %>% 
  pivot_longer(
    cols = "Methionine":"Glycine",
    names_to = "Type",
    values_to = "Gradient"
  )

metgly_long$X <- NULL
colnames(metgly_long)[which(names(metgly_long) == "Peptide.Length")] <- "Peptide Length"
metgly_long$Core.Class <- "None"
metgly_long$Core.Subclass <- "None"

ref_data$Peplength <- perf_data$Precursor.Length[match(ref_data$BGC.ID, perf_data$Name)&match(ref_data$Core.Number, perf_data$Core.Number)]
ref_data$Metabolite.list <- perf_data$Metabolite.list[match(interaction(ref_data$BGC.ID, perf_data$Name),interaction(ref_data$Core.Number,perf_data$Core.Number))]
ref_data <- ref_data[grepl("cyclic_peptide_c",ref_data$Metabolite.list,fixed=TRUE)|
                       grepl("macrocyclic_thiopeptide_c",ref_data$Metabolite.list,fixed=TRUE)|
                       grepl("macrolactam_peptide_c",ref_data$Metabolite.list,fixed=TRUE),]
ref_data <- ref_data[ref_data$Gradient < 1.0,]

ref_subset <- subset(ref_data, select=c(Gradient, Core.Class, Core.Subclass, Peplength))
colnames(ref_subset)[which(names(ref_subset) == "Peplength")] <- "Peptide Length"
ref_subset$Type <- "None"

bind_data <- rbind(metgly_long, ref_subset)
bind_data <- bind_data[bind_data$`Peptide Length`>15,]


ref_data$Core.Class[ref_data$Core.Class == "lassopeptide"] <- "Lasso Peptide" 
ref_data$Core.Class[ref_data$Core.Class == "thiopeptide"] <- "Thiopeptide" 
ref_data$Core.Class[ref_data$Core.Class == "lanthipeptide"] <- "Lanthipeptide" 
p <- ref_data %>%
  ggplot(aes(x=Peplength, y=Gradient, color=Core.Class, shape=Core.Subclass, group=interaction(Core.Class, Core.Subclass))) +
  geom_point(size=2) +
  theme(text = element_text(size=15)) +
  xlab("Precursor Peptide Length") +
  ylab("S. coelicolor Growth-Production Gradient")
show(p)  

bind_data$Type[bind_data$Type == "None"] <- "RiPP"
p <- bind_data %>%
  ggplot(aes(x=`Peptide Length`, y=Gradient, color=Type)) +
  geom_point(size=2) +
  theme(text = element_text(size=15)) +
  scale_color_discrete(name="Peptide") +
  xlab("Peptide Length") +
  ylab("S. coelicolor Growth-Production Gradient")
show(p)
