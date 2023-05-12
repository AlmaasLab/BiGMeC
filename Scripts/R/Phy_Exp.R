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
install.packages("webr")
library(webr)

src_data <- read.csv("Data\\CarveMe extended model RiPP results.csv")

length(unique(src_data$BGC.Source.Organism))

src_data <- na.omit(src_data)

src_data <- src_data[src_data$Max.NP!=0,]

src_data$Core.Class[src_data$Core.Class == "lassopeptide"] <- "Lasso Peptide" 
src_data$Core.Class[src_data$Core.Class == "thiopeptide"] <- "Thiopeptide" 
src_data$Core.Class[src_data$Core.Class == "lanthipeptide"] <- "Lanthipeptide" 

perf_data <- read.csv("RiPP Pathway Output\\RiPP reconstruction summary.csv")

#ref_data$Core.Number <- perf_data$Core.Number

src_data$Metabolite.list <- perf_data$Metabolite.list[which(perf_data$Name %in% src_data$BGC.ID & perf_data$Core.Number %in% src_data$Core.Number)]
src_data_mature <- src_data[grepl("cyclic_peptide_c",src_data$Metabolite.list,fixed=TRUE)|
                              grepl("macrocyclic_thiopeptide_c",src_data$Metabolite.list,fixed=TRUE)|
                              grepl("macrolactam_peptide_c",src_data$Metabolite.list,fixed=TRUE),]

src_data$Conc.Subclass <- paste(src_data$Core.Class, src_data$Core.Subclass)

#Pie chart of subclasses
src_na <- na.omit(src_data)
p <- src_na %>% ggplot(aes(x="", fill=factor(Clade3))) +
  geom_bar(width=1) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(title="RiPP Subclasses", 
       subtitle="Pie chart for distribution of RiPP Subclass for cores in dataset",
       caption="",
       x="",
       y="") +
  scale_fill_discrete(name="RiPP Subclass") +
  coord_polar(theta = "y", start=0)
show(p)

clades <- read.csv("Data\\clade_dataframe.csv")

src_data$Organism <- perf_data$Organism[match(src_data$BGC.ID, perf_data$Name)]
src_data$Clade1 <- clades$Clade.Level.1[match(src_data$Organism, clades$Organism.Name)]
src_data$Clade2 <- clades$Clade.Level.2[match(src_data$Organism, clades$Organism.Name)]
src_data$Clade3 <- clades$Clade.Level.3[match(src_data$Organism, clades$Organism.Name)]
src_na <- na.omit(src_data)
#Piedonut chart of clade level 2 and 3 of the source GMMs
p <- PieDonut(src_na, aes(Clade2, Clade3), showPieName = FALSE)
show(p)


#Piedonut chart of clade level 2 and 3 of every accessed organism with RiPP pathway
p <- PieDonut(clades, aes(Clade.Level.2, Clade.Level.3), showPieName = FALSE)
show(p)


src_data[is.na(src_data)] = "Streptomyces cattleya NRRL 8057 = DSM 46488"

src_data.summary <- src_data %>%
  group_by(Clade2) %>%
  summarise(
    sd = sd(Gradient, na.rm = TRUE),
    Gradient = mean(Gradient)
  )

p <- src_data.summary %>%
  ggplot(aes(x=Clade2, y=Gradient, fill=Clade2)) +
  geom_bar(stat="summary", fun="mean") +
  ylab("Growth-Production Gradient") +
  xlab("") +
  labs(x="") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  theme(text = element_text(size=15)) +
  scale_fill_discrete(name="Clade Level 2") +
  geom_errorbar(aes(ymin=Gradient-sd, ymax=Gradient+sd), width = 0.4, size=1,color="black")
show(p)


src_data.summary2 <- src_data %>%
  group_by(Clade3) %>%
  summarise(
    sd = sd(Gradient, na.rm = TRUE),
    Gradient = mean(Gradient)
  )

src_data.summary2$Clade3[src_data.summary2$Clade3 == ""] <- "Catenulispora acidiphila DSM 44928"

p <- src_data.summary2 %>%
  ggplot(aes(x=Clade3, y=Gradient, fill=Clade3)) +
  geom_bar(stat="summary", fun="mean") +
  ylab("Growth-Production Gradient") +
  xlab("") +
  labs(x="") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  theme(text = element_text(size=15), axis.title.x = element_blank()) +
  scale_fill_discrete(name="Clade Level 3") +
  geom_errorbar(aes(ymin=Gradient-sd, ymax=Gradient+sd), width = 0.4, size=1,color="black") +
show(p)



