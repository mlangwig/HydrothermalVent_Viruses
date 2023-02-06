##################### Integrating all virus output into 1 table #####################

setwd(dir = "~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/VentVirus_Analysis/")

library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)

##################### Read in major inputs #####################

iphop <- read.csv(file = "input/Host_prediction_to_genus_m90.csv", header = TRUE)

vib_amgs <- read.delim2(file = "input/All_AMG_annotations_Vents.tsv", header = TRUE)

checkV <- read.delim2(file = "input/CheckV_quality_vMAGs_vUnbinned.tsv", header = TRUE)

vMAG_mapping <- read.delim2(file = "input/vMAG_scaffold_mapping.txt", header = FALSE)  

vib_type <- read.delim2(file = "input/All_GenQuality_Vents.tsv", header = TRUE)

##################### Create the master table #####################

#select columns from vibrant amgs because has the scaffold and protein IDs
master_table<-select(vib_amgs, c('scaffold', 'protein', 'KO', 'AMG', 'KO.name',
                                 'Pfam', 'Pfam.name', 'VOG', 'VOG.name'))  

#map CheckV results
master_table <- checkV %>%
  dplyr::select(contig_id, contig_length, checkv_quality, provirus, completeness, contamination, warnings) %>%
  right_join(master_table, by = c("contig_id" = "scaffold"))   

#everything that was not mapped is a vMAG so they can now be easily separated

#make Site column
master_table$Site<-master_table$contig_id
master_table <- master_table %>% separate(Site, c("Site", NA), 
                                          sep= "_NODE|_scaffold")

#add vMAG data to the master table
vMAG_mapping<-rename(vMAG_mapping,"vMAG_name" = "V1")
vMAG_mapping<-rename(vMAG_mapping,"vMAG_scaffold" = "V2")

##map vMAG data (like VLOOKUP)
master_table <- vMAG_mapping %>%
  dplyr::select("vMAG_name", "vMAG_scaffold") %>%
  right_join(master_table, by = c("vMAG_scaffold" = "contig_id"))
master_table<-rename(master_table,"vMAG" = "vMAG_name")
master_table<-rename(master_table,"contig_id" = "vMAG_scaffold")

#map VIBRANT type (like VLOOKUP)
master_table <- vib_type %>%
  dplyr::select(scaffold, type) %>%
  right_join(master_table, by = c("scaffold" = "contig_id")) 

##################### Separating the vMAGs and unbinned viruses into 2 master tables #####################

##################### vMAG master table #####################






  