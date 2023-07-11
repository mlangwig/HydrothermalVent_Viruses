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

checkV <- read.delim2(file = "input/CheckV_quality_vMAGs_vUnbinned_PlumeVent.txt", header = TRUE)

vMAG_mapping <- read.delim2(file = "input/vMAG_scaffold_mapping.txt", header = FALSE)  

vib_type <- read.delim2(file = "input/All_GenQuality_Vents.tsv", header = TRUE)

gensize_kb<-read.delim2(file = "input/GenSize_KB.tsv", header = TRUE)

plume_master_vUnbinned<-read.delim2(file = "~/Google Drive/My Drive/Faith/PlumeViruses/PlumeVirus_Analysis/output/master_table_vUnbinned.tsv",
                                    header = TRUE)

plume_master_vMAGs<-read.delim2(file = "~/Google Drive/My Drive/Faith/PlumeViruses/PlumeVirus_Analysis/output/master_table_vMAGs.tsv",
                                    header = TRUE)

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

master_table_vMAGs <- master_table[is.na(master_table$checkv_quality),]

#because the counts in the vMAG table were off, I noticed 50 scaffolds from VIBRANT disappeared
#in the vMAGs. There are also 15 vMAG scaffolds that are new after binning so they aren't in the 
#original VIBRANT files. I assume these differences are due to Prodigal translations from DNA
#to protein and don't worry about it for the purposes of what I'm doing here. Detailed gene
#analyses will not use or use more than original VIBRANT annotations for vMAGs.

#removing scaffolds that are now gone after binning
scafs_to_remove<-read.delim2(file = "input/vMAG_scaffoldsToRemove.txt", header = FALSE)
scafs_to_remove<-rename(scafs_to_remove,"scaffold" = "V1")

master_table_vMAGs<-master_table_vMAGs %>% anti_join(scafs_to_remove, by = c("protein" = "V1"))
#notice now that the vMAG table is missing 15 protein scaffolds. This doesn't matter for the purposes
#of these plots/analyses but will be included for gene annotation-focused analyses

#remove old checkV
master_table_vMAGs <- select(master_table_vMAGs, -contig_length, -checkv_quality, -provirus,
                             -completeness, -contamination, -warnings)

#map new CheckV (VLOOKUP)
master_table_vMAGs <- checkV %>%
  dplyr::select(contig_id, contig_length, checkv_quality, provirus, completeness, contamination, warnings) %>%
  right_join(master_table_vMAGs, by = c("contig_id" = "vMAG")) 
#rename some columns
master_table_vMAGs<-rename(master_table_vMAGs, "vMAG"  = "contig_id")

##################### vUnbinned master table #####################

master_table_unbinned <- master_table[!is.na(master_table$checkv_quality),]
master_table_unbinned <- select(master_table_unbinned, -vMAG)

################# master tables without proteins for counting ###############
############## plus determing lytic lysogenic type for vMAGs ###############

#################vMAGs
master_table_vMAGs_noProtein<-master_table_vMAGs %>%
  select(vMAG, contig_length, checkv_quality, provirus, completeness, contamination, type, Site)
master_table_vMAGs_noProtein<-unique(master_table_vMAGs_noProtein)

plume_master_table_vMAGs_noProtein<-plume_master_vMAGs %>%
  select(vMAG, contig_length, checkv_quality, provirus, completeness, contamination, type, Site)
plume_master_table_vMAGs_noProtein<-unique(plume_master_table_vMAGs_noProtein)

#subset for vMAGs that have 2 type designations after binning
dups<-master_table_vMAGs_noProtein[duplicated(master_table_vMAGs_noProtein$vMAG),]
two_type_vMAGs<-master_table_vMAGs_noProtein[master_table_vMAGs_noProtein$vMAG %in% dups$vMAG,] #subset table using list of names

##PLUME MASTER TABLE WAS ALREADY CORRECTED FOR LYTIC/LYSOGENIC IN PLUME MASTER TABLE

###REPLACING ALL LYTIC/LYSOGENIC SAME SCAFFOLDS WITH LYSOGENIC ONLY ACCORDING TO VRHYME RULES
#now grab the two types from the master vMAG table so you can change all of them to lysogenic
two_type_to_one_vMAGs<-master_table_vMAGs[master_table_vMAGs$vMAG %in% dups$vMAG,]
two_type_to_one_vMAGs<-two_type_to_one_vMAGs %>% 
  mutate(type = str_replace(type, "lytic", "lysogenic"))

#remove double lytic/lyso type scaffolds from master vMAG df
`%notin%` <- Negate(`%in%`)
master_table_vMAGs<-master_table_vMAGs[master_table_vMAGs$vMAG %notin% dups$vMAG,] #subset table using list of names
#now add the vMAGs back in with their corrected VIBRANT type
master_table_vMAGs<-rbind(master_table_vMAGs, two_type_to_one_vMAGs)

master_table_vMAGs_noProtein<-master_table_vMAGs %>%
  select(vMAG, contig_length, checkv_quality, provirus, completeness, contamination, type, Site)
master_table_vMAGs_noProtein<-unique(master_table_vMAGs_noProtein)

##################vUnbinned
master_table_unbinned_noProtein<-master_table_unbinned %>%
  select(scaffold, contig_length, checkv_quality, provirus, completeness, contamination, type, Site)
master_table_unbinned_noProtein<-unique(master_table_unbinned_noProtein)

plume_master_table_unbinned_noProtein<-plume_master_vUnbinned %>%
  select(scaffold, contig_length, checkv_quality, provirus, completeness, contamination, type, Site)
plume_master_table_unbinned_noProtein<-unique(plume_master_table_unbinned_noProtein)

####################Combine the tables into 1 for plots and counts of all viruses together
#rename vMAG column to combine
master_table_vMAGs_noProtein<-rename(master_table_vMAGs_noProtein, Virus = vMAG)
master_table_unbinned_noProtein<-rename(master_table_unbinned_noProtein, Virus = scaffold)
#combine them
master_table_noProtein<-rbind(master_table_vMAGs_noProtein, master_table_unbinned_noProtein)

#rename vMAG column to combine
plume_master_table_vMAGs_noProtein<-rename(plume_master_table_vMAGs_noProtein, Virus = vMAG)
plume_master_table_unbinned_noProtein<-rename(plume_master_table_unbinned_noProtein, Virus = scaffold)
#combine them
master_table_noProtein<-rbind(master_table_vMAGs_noProtein, master_table_unbinned_noProtein,
                              plume_master_table_vMAGs_noProtein, plume_master_table_unbinned_noProtein)

############################## Replace site names with general name ###################################

#gsub to replace Plume names
master_table_noProtein$Site <- gsub(".*Lau_Basin.*","Lau_Basin",master_table_noProtein$Site) #the placement of the periods is crucial for replacing whole string
master_table_noProtein$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",master_table_noProtein$Site)
master_table_noProtein$Site <- gsub(".*Axial.*","Axial_Seamount",master_table_noProtein$Site)
master_table_noProtein$Site <- gsub(".*ELSC.*","Lau_Basin",master_table_noProtein$Site)
master_table_noProtein$Site <- gsub(".*Brothers.*","Brothers_Volcano",master_table_noProtein$Site)
master_table_noProtein$Site <- gsub(".*Guaymas.*","Guaymas_Basin",master_table_noProtein$Site)
master_table_noProtein$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",master_table_noProtein$Site)
master_table_noProtein$Site <- gsub(".*EPR.*","East_Pacific_Rise",master_table_noProtein$Site)

############################## Useful counts with no protein ###################################

#Count how many lytic and lysogenic predicted by VIBRANT

table(master_table_vMAGs_noProtein['type']) 
#VIBRANT - 823 lysogenic vMAGs, 4835 lytic
table(master_table_vMAGs_noProtein['provirus']) 
#CheckV - 141 lysogenic vMAGs, 5517 lytic

table(master_table_unbinned_noProtein['type']) 
#VIBRANT - 1146 lysogenic vMAGs, 24332 lytic
table(master_table_unbinned_noProtein['provirus']) 
#CheckV - 436 lysogenic vMAGs, 25042 lytic


#Count checkv_quality
table(master_table_vMAGs_noProtein['checkv_quality']) 
#CheckV --> Complete: 24, High-qual: 526, Med-qual: 1041, Low-qual: 3616, Not-determ: 451

table(master_table_unbinned_noProtein['checkv_quality'])
#CheckV --> Complete: 145, High-qual: 185, Med-qual: 418, Low-qual: 16727, Not-determ: 8003


########## Count for all viruses Vent and Plume #######
#LYTIC OR LYSOGENIC
table(master_table_noProtein['type'])
#VIBRANT - 2125 lysogenic, 35,889 lytic
table(master_table_noProtein['provirus'])
#CheckV - 668 lysogenic, 37,346 lytic

#CHECKV QUALITY
table(master_table_noProtein['checkv_quality'])
#CheckV --> Complete: 319, High-qual: 873, Med-qual: 1726, 
#Low-qual: 25549, Not-determ: 9547

################################# Mapping iphop to new tables ##################################
#New tables to avoid duplicates

master_table_unbinned_iphop <- iphop %>%
  dplyr::select(Virus, Host.genus, List.of.methods) %>%
  right_join(master_table_unbinned, by = c("Virus" = "scaffold")) 

master_table_vMAGs_iphop <- iphop %>%
  dplyr::select(Virus, Host.genus, List.of.methods) %>%
  right_join(master_table_vMAGs, by = c("Virus" = "vMAG")) 

##################### writing outputs #####################

write.table(master_table_vMAGs, file = "output/master_table_vMAGs.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(master_table_unbinned, file = "output/master_table_unbinned.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(gensize_kb, file = "output/gensize_VentPlume.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)

######################################## Visualize ##################################################

############################### CheckV quality ##################################

master_table_fig <- master_table_noProtein %>% 
  select(c("Virus", "checkv_quality", "Site")) %>%
  unique() %>%
  group_by(Site) %>%
  count(checkv_quality) #%>%
#filter(!grepl("Low-quality|Not-determined", 
#             checkv_quality)) #drop low quality bc it's telling me nothing

level_order <- c('Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete') 

###plot
dev.off()
p <- ggplot(master_table_fig, aes(x = factor(checkv_quality, levels = level_order), #reorder bc was plotting x axis backwards/upside down
                                  y = n, fill = checkv_quality)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Site) + 
  xlab(element_blank())  +
  ylab("Count") +
  ggtitle("Virus Genome Quality") +
  scale_fill_viridis_d(name = "CheckV Quality", direction = -1) +
  scale_y_continuous(breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000))
#geom_text(aes(label = paste0(n), y = n),
#          hjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 16)) +
  coord_flip()
p

ggsave("output/CheckV_qual_perSite.pdf", p, width = 15, height = 10, dpi = 500)
ggsave("output/CheckV_qual_perSite.png", p, width = 15, height = 10, dpi = 500)

############################### Lytic and Lysogenic ##################################

###subset just virus, site, type
master_table_fig <- master_table_noProtein %>% 
  select(c("Virus", "type", "Site"))

#make a column of 1s for presence
master_table_fig$value <- rep(c(1),each=nrow(master_table_fig)) 

dev.off()
p <- ggplot(master_table_fig, aes(x = type, y = value, fill = type)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Site) + 
  xlab("Type")  +
  ylab("Count") +
  ggtitle("Lytic and Lysogenic Viruses") +
  scale_fill_viridis_d(direction = -1, begin = .25, end = .75)
p <- p + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,  hjust = 1)) +
  coord_flip()
p

ggsave("output/LyticLysogenic.png", p, width = 15, height = 10, dpi = 500)
ggsave("output/LyticLysogenic.pdf", p, width = 15, height = 10, dpi = 500)

############################### Genome Size ##################################
gensize_kb<-read.delim2(file = "input/GenSize_KB.tsv", header = TRUE)

gensize_kb$Site<-gensize_kb$Genome
gensize_kb <- gensize_kb %>% separate(Site, c("Site", NA), 
                                      sep= "_NODE|_scaffold|_vRhyme|_k95")

gensize_kb$Site <- gsub(".*Lau_Basin.*","Lau_Basin",gensize_kb$Site) #the placement of the periods is crucial for replacing whole string
gensize_kb$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",gensize_kb$Site)
gensize_kb$Site <- gsub(".*Axial.*","Axial_Seamount",gensize_kb$Site)
gensize_kb$Site <- gsub(".*ELSC.*","Lau_Basin",gensize_kb$Site)
gensize_kb$Site <- gsub(".*Brothers.*","Brothers_Volcano",gensize_kb$Site)
gensize_kb$Site <- gsub(".*Guaymas.*","Guaymas_Basin",gensize_kb$Site)
gensize_kb$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",gensize_kb$Site)
gensize_kb$Site <- gsub(".*EPR.*","East_Pacific_Rise",gensize_kb$Site)

#map quality so can filter. SKIP THIS TO SEE GENOME SIZES WITHOUT QC FILTERING
gensize_kb <- checkV %>%
  dplyr::select(contig_id, checkv_quality, provirus, completeness, contamination, warnings) %>%
  right_join(gensize_kb, by = c("contig_id" = "Genome")) #%>%
#filter(!grepl("Low-quality|Not-determined", checkv_quality)) #drop low quality bc it's telling me nothing

gensize_kb$checkv_quality <- factor(gensize_kb$checkv_quality, levels=c("Complete", "High-quality", 
                                                                        "Medium-quality", "Low-quality",
                                                                        "Not-determined")) #set order of x axis 

#get the names of Vent viruses >=5kb
gensize_kb_5kb <- gensize_kb %>% filter(KB >= 5)
write.table(gensize_kb_5kb, file = "output/gensize_5kb.tsv", quote = FALSE, sep = "\t",
            col.names = TRUE, row.names = FALSE)

################### Genome Size by CheckV Quality
dev.off()
p <- ggplot(gensize_kb, aes(x = checkv_quality, y = as.numeric(KB), fill = checkv_quality)) + 
  geom_boxplot() + 
  facet_wrap(~Site) + 
  xlab("Site")  +
  ylab("Genome Size (KB)") +
  ggtitle("Viral Genome Size (KB)") +
  guides(fill=guide_legend(title = "CheckV Quality")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_x_discrete(limits = rev(levels(gensize_kb$checkv_quality))) +
  #scale_y_continuous(limits = c(0, 1050), breaks = seq(0, 1050, by = 300)) +
  scale_fill_viridis_d() +
  ylim(0,1050) +
  coord_flip()
p

ggsave("output/GenomeSizeKB.png", p, width = 12, height = 6, dpi = 500)
ggsave("output/GenomeSizeKB.pdf", p, width = 12, dpi = 500)

################### Genome Size Site Totals

master_table_fig <- gensize_kb %>% 
  select(c("contig_id", "KB", "Site")) %>%
  unique() %>%
  group_by(Site) %>%
  count(KB) #%>%

master_table_fig<-filter(master_table_fig, Site == "Brothers_Volcano")

dev.off()
p <- ggplot(master_table_fig, aes(x=n)) + 
  geom_histogram(binwidth = 5) +
  #facet_wrap(~Site) + 
  xlab("Count")  +
  ylab("Genome Size (KB)") +
  ggtitle("Viral Genome Size (KB)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  #guides(fill=guide_legend(title = "CheckV Quality")) +
  theme_bw()
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  #scale_x_discrete(limits = rev(levels(gensize_kb$checkv_quality))) +
  #scale_y_continuous(limits = c(0, 1050), breaks = seq(0, 1050, by = 300)) +
  #scale_fill_viridis_d() +
  #ylim(0,1050) +
  #coord_flip()
p

ggsave("output/GenomeSizeKB_SiteTotals.png", p, width = 12, height = 6, dpi = 500)
ggsave("output/GenomeSizeKB_SiteTotals.pdf", p, width = 12, dpi = 500)
