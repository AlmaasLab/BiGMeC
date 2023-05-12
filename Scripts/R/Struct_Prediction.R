library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(plyr)
library(tidyverse)
library(viridis)
library(forcats)
library(data.table)

PRISM_T <- read.csv("Data\\PRISM SMILES Tanimoto Scores.csv")
Thesis_T <- read.csv("Data\\RiPP performance summary.csv")
Thesis_T <- Thesis_T %>% group_by(Name) %>% slice(which.max(Prediction.Score))

Thesis_T$PRISM.Median.T <- PRISM_T$Median.J[match(Thesis_T$Name, PRISM_T$BGC.ID)]
Thesis_T$PRISM.Max.T <- PRISM_T$X90th.percentile.J[match(Thesis_T$Name, PRISM_T$BGC.ID)]
Thesis_T <- na.omit(Thesis_T)

length(Thesis_T$Prediction.Score[Thesis_T$Prediction.Score > 0.85])

length(Thesis_T$PRISM.Median.T[Thesis_T$PRISM.Median.T > 0.85])

length(Thesis_T$PRISM.Max.T[Thesis_T$PRISM.Max.T > 0.85])
length(Thesis_T$PRISM.Max.T[Thesis_T$Core.Class == "Lasso Peptide"&Thesis_T$PRISM.Max.T > 0.85])


length(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "lassopeptide"&Thesis_T$Prediction.Score > 0.85])
length(Thesis_T$Core.Class[Thesis_T$Core.Class == "lanthipeptide"])

length(Thesis_T$PRISM.Median.T[Thesis_T$Core.Class == "lassopeptide"&Thesis_T$PRISM.Median.T > 0.85])

mean(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "lassopeptide"])
mean(Thesis_T$Base.Prediction[Thesis_T$Core.Class == "lassopeptide"])
(mean(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "lassopeptide"])-mean(Thesis_T$Base.Prediction[Thesis_T$Core.Class == "lassopeptide"]))/mean(Thesis_T$Base.Prediction[Thesis_T$Core.Class == "lassopeptide"])

mean(Thesis_T$Prediction.Score)
mean(Thesis_T$PRISM.Median.T)
(mean(Thesis_T$Prediction.Score)-mean(Thesis_T$PRISM.Median.T))/mean(Thesis_T$PRISM.Median.T)

mean(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "lassopeptide"])
mean(Thesis_T$PRISM.Median.T[Thesis_T$Core.Class == "lassopeptide"])
(mean(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "lassopeptide"])-mean(Thesis_T$PRISM.Median.T[Thesis_T$Core.Class == "lassopeptide"]))/mean(Thesis_T$PRISM.Median.T[Thesis_T$Core.Class == "lassopeptide"])

#PRISM max tanimoto score values
mean(Thesis_T$Prediction.Score)
mean(Thesis_T$PRISM.Max.T)
(mean(Thesis_T$Prediction.Score)-mean(Thesis_T$PRISM.Max.T))/mean(Thesis_T$PRISM.Max.T)

#PRISM max tanimoto score values by class
mean(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "Thiopeptide"])
mean(Thesis_T$PRISM.Max.T[Thesis_T$Core.Class == "Thiopeptide"])
(mean(Thesis_T$Prediction.Score[Thesis_T$Core.Class == "Lanthipeptide"])-mean(Thesis_T$PRISM.Max.T[Thesis_T$Core.Class == "Lanthipeptide"]))/mean(Thesis_T$PRISM.Max.T[Thesis_T$Core.Class == "Lanthipeptide"])


data1 <- data.frame(
  type = c( rep("PRISM Median", 57), rep("ARMRiPP", 57) ),
  value = c( Thesis_T$PRISM.Median.T, Thesis_T$Prediction.Score )
)

data2 <- data.frame(
  type = c( rep("Core Peptide", 57), rep("ARMRiPP", 57) ),
  value = c( Thesis_T$Base.Prediction, Thesis_T$Prediction.Score )
)

data3 <- data.frame(
  type = c( rep("PRISM Max", 57), rep("ARMRiPP", 57) ),
  value = c( Thesis_T$PRISM.Max.T, Thesis_T$Prediction.Score )
)

p <- data2 %>%
  ggplot( aes(x=value, fill=type)) +
  geom_vline(aes(xintercept = 0.875), color="black", size=1)+
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity', binwidth = 0.05, bins=20) +
  scale_fill_manual(values=c("#a64d7a", "#76a5af")) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  xlab("Tanimoto Score") +
  ylab("Count")
show(p)

Thesis_T$Core.Class[Thesis_T$Core.Class == "lassopeptide"] <- "Lasso Peptide" 
Thesis_T$Core.Class[Thesis_T$Core.Class == "thiopeptide"] <- "Thiopeptide" 
Thesis_T$Core.Class[Thesis_T$Core.Class == "lanthipeptide"] <- "Lanthipeptide" 
Thesis_T$Conc.Subclass <- paste(Thesis_T$Core.Class, Thesis_T$Core.Subclass)


long_data <- Thesis_T %>% pivot_longer(cols=c("PRISM.Median.T", "Prediction.Score"),
                                       names_to = "T.Type",
                                       values_to = "T.Value")
long_data$T.Type[long_data$T.Type == "PRISM.Median.T"] <- "PRISM Median"
long_data$T.Type[long_data$T.Type == "Prediction.Score"] <- "ARMRiPP" 

long_data2 <- Thesis_T %>% pivot_longer(cols=c("Base.Prediction", "Prediction.Score"),
                                       names_to = "T.Type",
                                       values_to = "T.Value")
long_data2$T.Type[long_data2$T.Type == "Base.Prediction"] <- "Core Peptide"
long_data2$T.Type[long_data2$T.Type == "Prediction.Score"] <- "ARMRiPP" 

long_data3 <- Thesis_T %>% pivot_longer(cols=c("PRISM.Max.T", "Prediction.Score"),
                                        names_to = "T.Type",
                                        values_to = "T.Value")
long_data3$T.Type[long_data3$T.Type == "PRISM.Max.T"] <- "PRISM Max"
long_data3$T.Type[long_data3$T.Type == "Prediction.Score"] <- "ARMRiPP" 


p <- long_data3 %>%
  ggplot(aes(x=T.Value, fill = T.Type)) + 
  geom_vline(aes(xintercept = 0.875), color="black", size=1)+
  geom_histogram(color="#e9ecef", alpha=0.7, position = 'identity', binwidth = 0.05, bins=20) + 
  scale_fill_manual(values=c("#a64d7a", "#76a5af")) +
  theme(text = element_text(size = 20),legend.position="top", legend.box = "vertical") +
  labs(fill="") +
  xlab("Tanimoto Score") +
  ylab("Count") +
  facet_wrap(~Core.Class)
show(p)

table(Thesis_T$Conc.Subclass)
