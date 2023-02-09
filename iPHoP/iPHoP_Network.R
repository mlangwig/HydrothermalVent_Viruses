library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)

setwd("~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/iPHoP/")

######################################### Read the input ##################################################

#iphop
iphop <- read.csv(file = "input/Host_prediction_to_genus_m90.csv", header = TRUE)
#653 unique hosts without any filtering

#CheckV
checkv <- read.delim(file = "input/CheckV_quality_vMAGs_vUnbinned.tsv", header = TRUE, sep = "\t")

#Genome Size
gensize_kb <- read.delim(file = "input/Vent_vUnbinned_vMAG_GenSize_KB.tsv", header = TRUE, sep = "\t")

MAG_tax<-read.delim(file = "input/gtdbtk_v1.5.0_VentMAGs.tsv", header = TRUE)

################################## Map quality data so you can filter ########################################

#map
iphop <- checkv %>%
  dplyr::select(contig_id, contig_length, checkv_quality, provirus, gene_count, viral_genes,
                completeness, contamination, warnings) %>%
  right_join(iphop, by = c("contig_id" = "Virus"))

#map
iphop <- gensize_kb %>%
  dplyr::select(Genome, KB) %>%
  right_join(iphop, by = c("Genome" = "contig_id"))

####################### Remove iphop results that don't match MAG taxonomy ########################################
#I will only keep predictions that match the MAG data that I have

#remove ;s_ in gtdbtk classification for mapping
MAG_tax <- MAG_tax %>% separate(classification, c("classification", NA), sep= "(?=;s__)")
#vlookup mapping
iphop <- MAG_tax %>%
  dplyr::select(classification, user_genome) %>%
  right_join(iphop, by = c("classification" = "Host.genus"))
#drop NAs
iphop <- iphop %>%
  drop_na(user_genome)
#rename classification column
iphop<-rename(iphop,"Host.genus" = "classification")
iphop<-rename(iphop,"Virus" = "Genome")
#258 unique hosts when filter by matching MAG taxonomy

################################### Quality control iphop results ########################################
#I am removing viruses whose checkv quality was Not determined because in my manual inspections,
#this gets rid of a lot of junk

iphop<-iphop[!grepl("Not-determined", iphop$checkv_quality),]

#Removing viruses with the warning "no viral genes detected" because my manual inspections suggest these are not viral
#or are poor enough quality that I don't want to keep

iphop<-iphop[!grepl("no viral genes detected", iphop$warnings),]

#Now removing viruses ≤5kb because I am not sure I trust host predictions to viral fragments
#And I'd like the potential for more genomic context from the virus

#remove , from KB to prevent issues
iphop<-subset(iphop, iphop$KB>=5000)

#Remove viruses with contamination >20% because these don't look great

iphop<-subset(iphop, iphop$contamination<=20)
#194 unique hosts when filtering by all these quality metrics

##########################remove user_genome so I can see the table##########################
iphop<-select(iphop, -c("user_genome"))
iphop<-unique(iphop)



