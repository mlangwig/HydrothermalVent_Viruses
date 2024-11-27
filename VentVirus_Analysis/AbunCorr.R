################################### Process abundance data for correlation ####################################
library(dplyr)
library(tidyverse)
library(knitr)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(viridis)
library(pals)
library(RColorBrewer)
library(stringi)
library(ggbreak)

#read input
vir_abun <- read_delim(file = "input/PlumeVentVirus-vs-Reads-CoverM-Count_MinCov.tsv")
mic_abun <- read_delim(file = "input/PlumeVentMicrobialMAGs-vs-Reads-CoverM-Count_MinCov.tsv")

#remove unmapped and unnecessary columns
vir_abun <- vir_abun[-1,]
mic_abun <- mic_abun[-1,]

################################## Grabbing relative abundance ##################################

# #grab relative abundance value
# vir_abun <- vir_abun %>%
#   select(matches("Genome|Relative"))
# 
# mic_abun <- mic_abun %>%
#   select(matches("Genome|Relative"))


#or grab read counts
vir_abun <- vir_abun %>%
  select(matches("Genome|Read Count"))

mic_abun <- mic_abun %>%
  select(matches("Genome|Read Count"))

#add microbe taxonomy to sum abundance by it
mic_abun <- gtdb %>%
  select(c("user_genome", "c")) %>%
  right_join(mic_abun, by = c("user_genome" = "Genome"))
#filter for just Gamma and Camp
mic_abun <- mic_abun %>%
  filter(c %in% c("c__Gammaproteobacteria", "c__Campylobacteria"))

# #group by and sum abundance
# mic_abun_sum <- mic_abun %>%
#   group_by(c) %>%
#   summarize(across(contains("Relative"), sum, na.rm = TRUE))

#group by and sum abundance
mic_abun_sum <- mic_abun %>%
  group_by(c) %>%
  summarize(across(contains("Read Count"), sum, na.rm = TRUE))

#add virus host prediction to sum abundance by it
iphop_gc <- iphop %>%
  separate("Host.genus", c("d", "p", "c", "o", "f", "g"), sep= ";")
#filter for just Gamma and Camp infecting
iphop_gc <- iphop_gc %>%
  filter(c %in% c("c__Gammaproteobacteria", "c__Campylobacteria"))
#remove redundancy
iphop_gc <- iphop_gc %>%
  select(c("Virus", "c")) %>%
  unique()
#print the viruses that have multiple predictions
length(unique(iphop_gc$Virus))
#2567 tells me there are viruses that occur more than once
iphop_gc %>%
  count(Virus) %>%
  filter(n > 1) %>%
  pull(Virus) %>%
  print()
#this shows 2 viruses that have multiple predictions, Brothers_NWCB_S012_vRhyme_bin_152 and Brothers_NWCB_S141_vRhyme_bin_38
# bin 152 has CRISPR to Campylo so I'm keeping that one
# bin 38 has CRISPR to Gamma so I'm keeping that one
iphop_gc <- iphop_gc %>%
  filter(!(Virus == "Brothers_NWCB_S012_vRhyme_bin_152" & c == "c__Gammaproteobacteria")) %>%
  filter(!(Virus == "Brothers_NWCB_S141_vRhyme_bin_38" & c == "c__Campylobacteria"))
#check that viruses are now unique
length(unique(iphop_gc$Virus))
length(iphop_gc$Virus)

#map microbe the virus infects and only keep Gammas and Camps
vir_abun <- iphop_gc %>%
  inner_join(vir_abun, by = c("Virus" = "Genome"))

#erase unique virus ID
vir_abun2 <- vir_abun
vir_abun2$Virus <- "virus"
#join names
vir_abun2$virus_host <- paste(vir_abun2$Virus, vir_abun2$c, sep = "_")
# Move the merged column to be the first column
vir_abun2 <- vir_abun2 %>%
  select(virus_host, everything())

#group by and sum abundance
vir_abun_sum <- vir_abun2 %>%
  group_by(virus_host) %>%
  summarize(across(contains("Read Count"), sum, na.rm = TRUE))

#remove first column name
colnames(vir_abun_sum)[1] <- ""
colnames(mic_abun_sum)[1] <- ""

#fix col names
# #remove rel abun from colnames
# colnames(vir_abun_sum) <- gsub(" Relative Abundance \\(\\%)", "", colnames(vir_abun_sum))
# colnames(mic_abun_sum) <- gsub(" Relative Abundance \\(\\%)", "", colnames(mic_abun_sum))

#remove rel abun from colnames
colnames(vir_abun_sum) <- gsub(" Read Count", "", colnames(vir_abun_sum))
colnames(mic_abun_sum) <- gsub(" Read Count", "", colnames(mic_abun_sum))

#fix abun_names
read_mapping <- abun_names %>%
  select(c("Reads", "Site"))
read_mapping$Site <- gsub("_metaspades_scaffolds.min1000.fasta","", read_mapping$Site)
read_mapping$Site <- gsub(".fasta","", read_mapping$Site)
read_mapping$Site <- gsub("min1000","Plume", read_mapping$Site)
read_mapping$Site <- gsub("_scaffolds_Plume","", read_mapping$Site) 
read_mapping$Site <- gsub("Bowl","Mariner", read_mapping$Site)
read_mapping$Site <- gsub("ELSC","Lau Basin Deposit",read_mapping$Site)
read_mapping$Site <- gsub("Abe","ABE",read_mapping$Site)

#replace . with -
read_mapping$Reads <- gsub("\\.fastq\\.gz$", "", read_mapping$Reads)  # Remove .fastq.gz
read_mapping$Reads <- gsub("\\.", "-", read_mapping$Reads)  # Replace periods with hyphens
read_mapping$Reads <- paste0(read_mapping$Reads, ".fastq.gz")  # Add .fastq.gz back

#create a named vector
name_mapping <- setNames(read_mapping$Site, read_mapping$Reads)

#map
#tst <- mic_abun_sum

# Replace the column names in df2 with the simpler names based on matching

#virus
colnames(vir_abun_sum) <- sapply(colnames(vir_abun_sum), function(x) ifelse(x %in% names(name_mapping), 
                                                          name_mapping[x], 
                                                          x))
#microbe
colnames(mic_abun_sum) <- sapply(colnames(mic_abun_sum), function(x) ifelse(x %in% names(name_mapping), 
                                                                            name_mapping[x], 
                                                                            x))

#write output
write_csv(vir_abun_sum, file = "output/class_VentVirus_abundance_GamCamp_RC.csv")
write_csv(mic_abun_sum, file = "output/class_VentMicrobe_abundance_GamCamp_RC.csv")

########################### Grabbing reads mapped normalized by # reads in sample ###########################


########################### propr corr analysis with example data ###########################

#https://github.com/tpq/propr

#counts <- matrix(rpois(20*50, 100), 20, 50)
#group <- sample(c("A", "B"), size = 20, replace = TRUE)
devtools::install_github("tpq/propr")
library(propr)

# example
# Sample input count data
data <- matrix(c(10, 5, 15, 20, 30, 25), nrow = 2, byrow = TRUE)
# Calculate Propr matrix using correlation coefficient
result_cor <- propr(data, metric = "cor", ivar = "clr")
getMatrix(result_cor)
# Calculate Propr matrix using variance of log-ratio (VLR)
result_vlr <- propr(data, metric = "vlr", ivar = "clr")
getMatrix(result_vlr)
# Calculate Propr matrix using partial correlation coefficient
result_pcor <- propr(data, metric = "pcor", ivar = "clr")
getMatrix(result_pcor)
# Calculate Propr matrix using phi
result_phi <- propr(data, metric = "phi", ivar = "clr")
getMatrix(result_phi)


########################### propr corr analysis with my data ###########################
# my data
# from the github documentation: 
#"counts,  # rows as samples, like it should be" aka input data rows are samples

#virus
data_virus <- vir_abun_sum
#transpose data
data_virus_tr <- as.data.frame(t(data_virus))
# Set the first row as column names
colnames(data_virus_tr) <- data_virus_tr[1, ]
# Remove the first row
data_virus_tr <- data_virus_tr[-1, ]

#microbe
data_microbe <- mic_abun_sum
#transpose data
data_microbe_tr <- as.data.frame(t(data_microbe))
# Set the first row as column names
colnames(data_microbe_tr) <- data_microbe_tr[1, ]
# Remove the first row
data_microbe_tr <- data_microbe_tr[-1, ]

#put together
data_mic_vir <- cbind(data_microbe_tr, data_virus_tr)
#convert to numeric
data_mic_vir <- type.convert(data_mic_vir, as.is = TRUE)

#write output
write_delim(data_mic_vir, file = "output/class_VentMicrobe_abundance_GamCamp_transpose_RC.tsv")


# Calculate Propr matrix using correlation coefficient
result_rho <- propr(data_mic_vir, metric = "rho", ivar = "clr")
getMatrix(result_rho)
# # Calculate Propr matrix using correlation coefficient
# result_cor <- propr(data, metric = "cor", ivar = "clr")
# getMatrix(result_cor)
# # Calculate Propr matrix using variance of log-ratio (VLR)
# result_vlr <- propr(data, metric = "vlr", ivar = "clr")
# getMatrix(result_vlr)
# # Calculate Propr matrix using partial correlation coefficient
# result_pcor <- propr(data, metric = "pcor", ivar = "clr")
# getMatrix(result_pcor)
# Calculate Propr matrix using phi
result_phi <- propr(data_mic_vir, metric = "phi", ivar = "clr")
getMatrix(result_phi)

#determine significance of proportionality?
result_rho <- updateCutoffs.propr(
  result_rho,
  number_of_cutoffs = 100,  # number of cutoffs to estimate FDR
  custom_cutoffs = NULL,  # or specify custom cutoffs
  tails = 'right',  # consider only the positive values ('right') or both sides ('both')
  ncores = 1  # parallelize here
) 
#result_rho <- updateCutoffs(result_rho, cutoff = seq(0, 1, .05))
#get results
getAdjacencyFDR(result_rho)
propr::getCutoffFDR(result_rho)


########################### Descriptive analysis of compositional data ###########################

#Following the guide of https://www.geo.fu-berlin.de/en/v/soga-r/Advances-statistics/Feature-scales/Descriptive_analysis_of_compositional_data/index.html

########################### example #################################

# load data set
g36 <- read.csv("http://userpages.fu-berlin.de/soga/data/raw-data/G36chemical.txt",
                sep='\t',
                row.names = 'Sample')
# exclude columns of no interest 
g36 <- g36[, -c(1:5, ncol(g36)-1, ncol(g36))]
# rename columns for better readability
colnames(g36) <- gsub(".mg.g", "", colnames(g36))

# define a subcomposition
library(compositions)
X = acomp(g36, c("Fe","K","Mg"))
mean(X)
mean(X, robust = TRUE)
mvar(X)
msd(X)
variation(X)


########################### my data #################################

X2 = acomp(data_mic_vir)
mean(X2) #compositional mean, not robust
mean(X2, robust = TRUE) #robust estimation of compositional mean w/ composition package
mvar(X2) #metric variation, total or generalized variance. Average squared distance from center of data.
msd(X2) #metric standard deviation, radial standard deviation on a log scale if the variation was the same in all directions. Averaged spread.

variation(X2) #smaller variation = better proportionality between two components


