#Sets working directory if using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
setwd("..")

#Sets working directory for other R IDE's (don't run in RStudio)
setwd(getSrcDirectory(function(){})[1])
setwd("..")
setwd("..")

#Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(plyr)
library(tidyverse)
library(viridis)
library(forcats)

#Loads reference and carveme GEM results
ref_data <- read.csv("Data\\Sco-GEM reference model RiPP results.csv")
src_data <- read.csv("Data\\CarveMe extended model RiPP results.csv")

src_data <- na.omit(src_data)

src_data <- src_data[src_data$Max.NP!=0,]

src_data$ref.gradient <- ref_data$Gradient[which(ref_data$BGC.ID %in% src_data$BGC.ID & ref_data$Core.Number %in% src_data$Core.Number)]
src_data$peplength <- ref_data$Peplength[which(ref_data$BGC.ID %in% src_data$BGC.ID & ref_data$Core.Number %in% src_data$Core.Number)]
src_data$Core.Class[src_data$Core.Class == "lassopeptide"] <- "Lasso Peptide" 
src_data$Core.Class[src_data$Core.Class == "thiopeptide"] <- "Thiopeptide" 
src_data$Core.Class[src_data$Core.Class == "lanthipeptide"] <- "Lanthipeptide" 
src_data$Conc.Subclass <- paste(src_data$Core.Class, src_data$Core.Subclass)

#Determine clades
perf_data <- read.csv("RiPP Pathway Output\\RiPP reconstruction summary.csv")
clades <- read.csv("Data\\clade_dataframe.csv")

src_data$Organism <- perf_data$Organism[match(src_data$BGC.ID, perf_data$Name)]
src_data$Clade1 <- clades$Clade.Level.1[match(src_data$Organism, clades$Organism.Name)]
src_data$Clade2 <- clades$Clade.Level.2[match(src_data$Organism, clades$Organism.Name)]
src_data$Clade3 <- clades$Clade.Level.3[match(src_data$Organism, clades$Organism.Name)]

#Unique organisms
length(unique(src_data$BGC.Source.Organism))

#Names organisms in uncategorized clades
src_data[is.na(src_data)] = "Streptomyces cattleya NRRL 8057 = DSM 46488"
src_data$Clade3[src_data$Clade3 == ""] <- "Catenulispora acidiphila DSM 44928"

#Metrics of gradient values in native host (Gradient) and reference (ref.gradient)
mean(src_data$Gradient)
sd(src_data$Gradient)
mean(src_data$ref.gradient)
sd(src_data$ref.gradient)

#T-test for determining of there is a significant difference
var.test(src_data$Gradient, src_data$ref.gradient)
t.test(src_data$Gradient, src_data$ref.gradient, var.equal=TRUE)

#Scatter plot by clades
p <- src_data %>%
  ggplot(aes(x=Gradient, y=ref.gradient, color=Clade3, shape=Clade2, group=interaction(Clade2, Clade3))) +
  geom_point(size=2.5) +
  geom_abline(intercept = 0, slope = 1, size = 1.25) +
  xlab("Native Host Growth-Production Gradient") +
  ylab("S. coelicolor Growth-Production Gradient") + 
  scale_color_discrete(name="Clade Level 3") +
  scale_shape_discrete(name="Clade Level 2") +
  theme(text = element_text(size=15),legend.position= c(0.8,0.4)) +
  expand_limits(x = 0, y = 0)
show(p)

#Scatter plot by RiPP class
p <- src_data %>%
  ggplot(aes(x=Gradient, y=ref.gradient, color=Core.Subclass, shape=Core.Class, group=interaction(Core.Class, Core.Subclass))) +
  geom_point(size=2.5) +
  geom_abline(intercept = 0, slope = 1, size = 1.25) +
  xlab("Native Host Growth-Production Gradient") +
  ylab("S. coelicolor Growth-Production Gradient") + 
  scale_color_discrete(name="Core Subclass") +
  scale_shape_discrete(name="Core Class") +
  theme(text = element_text(size=15),legend.position= c(0.8,0.4)) +
  expand_limits(x = 0, y = 0)
show(p) 

#Peptide length metrics
mean(src_data$peplength)
sd(src_data$peplength)
median(src_data$peplength)
