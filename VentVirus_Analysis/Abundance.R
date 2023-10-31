########################## Compiling and parsing virus CoverM abundance ############################
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(viridis)
library(pals)
library(RColorBrewer)
########################## inputs ############################

#read input
abun <- read.delim(file = "input/PlumeVentVirus-vs-Reads-CoverM.tsv", header = TRUE)
#read in metadata to change names
abun_names <- read.delim2(file = "input/AssembliesToReads_Mapping.txt", header = FALSE)

########################## name cleaning ############################
#clean names
colnames(abun) = gsub("gz.Relative.Abundance....", "gz", colnames(abun))
#drop SWIR
abun <- abun %>% select(-c("X58P_trim_1.fastq.gz", "SWIR_B_trim_1.fastq.gz"))

#remove .fasta
#from whole df
#abun_names <- as.data.frame(apply(abun_names,2, str_remove_all, "_metaspades_scaffolds.min1000.fasta|.fasta"))
abun_names$V1 <- gsub("_metaspades_scaffolds.min1000.fasta|.fasta","", abun_names$V1) 
#change dashes to dots
abun_names$V2 <- gsub("-","\\.", abun_names$V2) 

########################## add Site name and melt long ############################

#transform from wide to long
abun_long <- melt(abun, id = c("Genome")) 
#rename variable column to Site
abun_long<-rename(abun_long, "Site" = "variable")
#drop unmapped
abun_long<-abun_long[!grepl("unmapped", abun_long$Genome),]

#change the read names to site names
abun_long <- abun_names %>%
  dplyr::select("V1", "V2") %>%
  right_join(abun_long, by = c("V2" = "Site"))
abun_long<-rename(abun_long,"Site" = "V1") %>%
  select(-"V2")

#drop 0s
abun_long<-filter(abun_long, value > 0)

########################## add host metadata ############################

abun_long_iphop <- iphop %>%
  dplyr::select(Virus, Host.genus) %>%
  right_join(abun_long, by = c("Virus" = "Genome"))   
#only keep phylum and class of the tax string
abun_long_iphop <- abun_long_iphop %>% separate(Host.genus, c("d", "p", "c", "o", "f", "g"), 
                                        sep= ";")
#select only phylum and class taxonomy
abun_long_iphop <- abun_long_iphop %>% select(c("Virus","p","c","Site","value"))

#for Proteobacteria keep class, everything else, keep phylum
abun_long_proteo <- abun_long_iphop %>% filter(grepl("p__Proteobacteria", p))
abun_long_proteo <- abun_long_proteo %>% select(-c("p"))
abun_long_proteo <- abun_long_proteo %>% rename("Taxa" = "c")
abun_long_proteo <- abun_long_proteo %>% select(c("Virus","Taxa","Site","value"))


abun_long_iphop <- abun_long_iphop %>% filter(!grepl("p__Proteobacteria", p))
abun_long_iphop <- abun_long_iphop %>% select(-c("c"))
abun_long_iphop <- abun_long_iphop %>% rename("Taxa" = "p")
#put the data frames back together
abun_long_iphop <- rbind(abun_long_iphop, abun_long_proteo)


############################ Set order of sites for plotting #############################
#set order of x axis 
# coverm_long_final$SiteFull <- factor(coverm_long_final$SiteFull, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
#                                                                           "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
#                                                                           "13_July", "13_Oct", "13_Jan", "13_May",
#                                                                           "21_July", "21_Oct", "21_Jan", "21_May",
#                                                                           "24_July", "24_Oct", "24_Jan", "24_May")) 








