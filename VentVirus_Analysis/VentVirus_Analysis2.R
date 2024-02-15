##################### Integrating all virus output into 1 table #####################

setwd(dir = "~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/VentVirus_Analysis/")

library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)

##################### Read in major inputs #####################

#iphop

#checkV

#lifestyle

gensize <- read.table(file = "input/PlumeVentVirus_Seqkit.txt", sep = " ", header = TRUE)

vib_annos <- read.delim(file = "../../VIBRANT_annos/52samples_VIBRANT_annos.tsv", sep = "\t", fill = TRUE, 
                        header = TRUE, na.strings=c("","NA"))

vMAG_mapping <- read.table(file = "input/PlumeVent_vMAG_scaffold_mapping.txt", sep = "\t")
  vMAG_mapping <- vMAG_mapping %>%
    rename("vMAG" = "V1",
           "scaffold" = "V2") #new = old
  
  
