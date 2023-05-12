#Sets working directory if using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
setwd("..")

#Sets working directory for other R IDE's
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

#Reads heterologous experiment results
het_data <- read.csv("Data\\CarveMe extended heterologous model RiPP results.csv")

#Gives number of unique BGCs
length(unique(het_data$Host.Model.BGC))

#Gives data of native host results
src_data <- read.csv("Data\\CarveMe extended model RiPP results.csv")
src_data <- na.omit(src_data)

#Merges so that we get the gradient of the native and heterologous host
merge <- left_join(het_data,src_data, by=c("BGC.ID"="BGC.ID", "Core.Number"="Core.Number"))
het_data$Source.Gradient <- src_data$Gradient[match(interaction(het_data$BGC.ID, het_data$Core.Number), interaction(src_data$BGC.ID, src_data$Core.Number))]
het_data <- na.omit(het_data)
het_data <- het_data[!het_data$Core.Class=="thiopeptide",]

#Reads in summary file to get organism names
perf_data <- read.csv("RiPP Pathway Output\\RiPP reconstruction summary.csv")
clades <- read.csv("Data\\clade_dataframe.csv")

#Defines source organism by BGC id, and 3 first clades of common tree
het_data$Source.Organism <- perf_data$Organism[match(het_data$BGC.ID, perf_data$Name)]
het_data$Source.Clade1 <- clades$Clade.Level.1[match(het_data$Source.Organism, clades$Organism.Name)]
het_data$Source.Clade2 <- clades$Clade.Level.2[match(het_data$Source.Organism, clades$Organism.Name)]
het_data$Source.Clade3 <- clades$Clade.Level.3[match(het_data$Source.Organism, clades$Organism.Name)]

#Defines host organism by BGC id, and 3 first clades of common tree
het_data$Host.Organism <- perf_data$Organism[match(het_data$Host.Model.BGC, perf_data$Name)]
het_data$Host.Clade1 <- clades$Clade.Level.1[match(het_data$Host.Organism, clades$Organism.Name)]
het_data$Host.Clade2 <- clades$Clade.Level.2[match(het_data$Host.Organism, clades$Organism.Name)]
het_data$Host.Clade3 <- clades$Clade.Level.3[match(het_data$Host.Organism, clades$Organism.Name)]
het_data <- na.omit(het_data)


#Determines groups of phylogenetic similarity between native and heterologous host
het_data <- het_data %>%
  mutate(Clade2sim = if_else(
    Host.Clade2 == Source.Clade2, "Same", "Different"
  ))
het_data <- het_data %>%
  mutate(Clade3sim = if_else(
    Host.Clade3 == Source.Clade3, "Same", "Different"
  ))

#Scatter plot
p <- het_data %>%
  ggplot(aes(x=Source.Gradient, y=Gradient, color=Clade2sim, shape=Clade3sim, group=interaction(Clade2sim, Clade3sim))) +
  geom_point(size=2.5) +
  geom_abline(intercept = 0, slope = 1, size=1.25) +
  scale_color_discrete(name="Clade 2") +
  scale_shape_discrete(name="Clade 3") +
  theme(text = element_text(size=15),legend.position="top", legend.box = "vertical") +
  xlab("Native Host Growth-Production Gradient") +
  ylab("Heterologous Host Growth-Production Gradient") +
  expand_limits(x = 0, y = 0)
show(p)

#Remove S. humidus outlier
het_data <- het_data[het_data$Source.Organism != "Streptomyces humidus",]

#For each Welch Two sample T-test, a Fischer's homoscedasticity test is performed, if p-value < 0.05 the equal variance parameter for the T.test is set to FALSE, and TRUE otherwise.
#T. value tests for gradients
var.test(het_data$Source.Gradient, het_data$Gradient)
t.test(het_data$Source.Gradient, het_data$Gradient, var.equal=TRUE)

#T. value tests for gradients of same second clade
var.test(het_data$Source.Gradient[het_data$Clade2sim == "Same"], het_data$Gradient[het_data$Clade2sim == "Same"])
t.test(het_data$Source.Gradient[het_data$Clade2sim == "Same"], het_data$Gradient[het_data$Clade2sim == "Same"], var.equal=TRUE)

#T. value tests for gradients of same second and third clade
var.test(het_data$Source.Gradient[het_data$Clade3sim == "Same"], het_data$Gradient[het_data$Clade3sim == "Same"])
t.test(het_data$Source.Gradient[het_data$Clade3sim == "Same"], het_data$Gradient[het_data$Clade3sim == "Same"], var.equal=TRUE)

#T. value tests for gradients of different second and third clade
var.test(het_data$Source.Gradient[het_data$Clade2sim == "Different"], het_data$Gradient[het_data$Clade2sim == "Different"])
t.test(het_data$Source.Gradient[het_data$Clade2sim == "Different"], het_data$Gradient[het_data$Clade2sim == "Different"], var.equal=TRUE)

#Metrics of source and host gradient values for group of different clades (different clade 2 => different clade 3)
mean(het_data$Source.Gradient[het_data$Clade2sim == "Different"])
sd(het_data$Source.Gradient[het_data$Clade2sim == "Different"])
mean(het_data$Gradient[het_data$Clade2sim == "Different"])
sd(het_data$Gradient[het_data$Clade2sim == "Different"])

#T. value tests for heterologous gradient of different second and third clade and heterologous gradient of same second and third clade
var.test(het_data$Gradient[het_data$Clade3sim == "Same"], het_data$Gradient[het_data$Clade2sim == "Different"])
t.test(het_data$Gradient[het_data$Clade3sim == "Same"], het_data$Gradient[het_data$Clade2sim == "Different"], var.equal=FALSE)

#T. value tests for native gradient of different second and third clade and native gradient of same second and third clade
var.test(het_data$Source.Gradient[het_data$Clade3sim == "Same"], het_data$Source.Gradient[het_data$Clade2sim == "Different"])
t.test(het_data$Source.Gradient[het_data$Clade3sim == "Same"], het_data$Source.Gradient[het_data$Clade2sim == "Different"], var.equal=FALSE)

