####### Determine the proportion of Plume:Vent viruses in clusters #######
library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)
library(reshape2)

######################################### Read Input ################################################
dRep_clusters <- read.csv2(file = "Input/Cdb.csv", sep = ",")
dRep_clusters_5kb <- read.csv2(file = "Input/Cdb_VentPlume_5kb_90ani.csv", sep = ",")

mcl_clusters <- read.csv2(file = "Input/vUnbinned_vMAGs_PlumeVents.clusters.csv", sep = ",", header = FALSE)

############################################ dRep ####################################################

#### Transforming names so can get proportions
#separate names by NODE, scaffold, and vRhyme and remove that part of the name
dRep_clusters <- dRep_clusters %>% separate(genome, c("Site", NA), sep= "(?=_NODE|_scaffold|__|_k95)")
unique(dRep_clusters$Site)

dRep_clusters_5kb <- dRep_clusters_5kb %>% separate(genome, c("Site", NA), sep= "(?=_NODE|_scaffold|_vRhyme|_k95)")
unique(dRep_clusters_5kb$Site)
#gsub to replace Plume names
dRep_clusters$Site <- gsub(".*Lau_Basin.*","Plume",dRep_clusters$Site) #the placement of the periods is crucial for replacing whole string
dRep_clusters$Site <- gsub(".*Cayman.*","Plume",dRep_clusters$Site)
dRep_clusters$Site <- gsub(".*Guaymas_Basin.*","Plume",dRep_clusters$Site)
dRep_clusters$Site <- gsub(".*Axial.*","Plume",dRep_clusters$Site)

dRep_clusters_5kb$Site <- gsub(".*Lau_Basin.*","Plume",dRep_clusters_5kb$Site) #the placement of the periods is crucial for replacing whole string
dRep_clusters_5kb$Site <- gsub(".*Cayman.*","Plume",dRep_clusters_5kb$Site)
dRep_clusters_5kb$Site <- gsub(".*Guaymas_Basin.*","Plume",dRep_clusters_5kb$Site)
dRep_clusters_5kb$Site <- gsub(".*Axial.*","Plume",dRep_clusters_5kb$Site)
#separate the plume and vent dataframes to make life easier
dRep_clusters_Plume <- dRep_clusters %>% filter(Site == "Plume")

dRep_clusters_Plume_5kb <- dRep_clusters_5kb %>% filter(Site == "Plume")
#remove Plume from dRep so can make all names left Vent
dRep_clusters <- dRep_clusters %>% filter(Site != "Plume")
dRep_clusters$Site <- "Vent"

dRep_clusters_5kb <- dRep_clusters_5kb %>% filter(Site != "Plume")
dRep_clusters_5kb$Site <- "Vent"
#Now put them back together
dRep_clusters <- rbind(dRep_clusters, dRep_clusters_Plume)

dRep_clusters_5kb <- rbind(dRep_clusters_5kb, dRep_clusters_Plume_5kb)
#Add counts to filter for clusters that occur more than once (not a singleton)
dRep_clusters_noSingle <- dRep_clusters %>%
  add_count(secondary_cluster) %>% 
  filter(n!=1) %>%
  select(-n)

dRep_clusters_noSingle_5kb <- dRep_clusters_5kb %>%
  add_count(secondary_cluster) %>% 
  filter(n!=1) %>%
  select(-n)

#number of clusters with more than 1 rep according to dRep
length(unique(dRep_clusters_noSingle$secondary_cluster))
#641
length(unique(dRep_clusters_noSingle_5kb$secondary_cluster))
#156
#184 with 5kb and 90% ANI

#sort by unique
dRep_clusters_noSingle <- unique(dRep_clusters_noSingle)

dRep_clusters_noSingle_5kb <- unique(dRep_clusters_noSingle_5kb)
#group by secondary cluster, count occurrences of Site
dRep_count <- dRep_clusters_noSingle %>% group_by(secondary_cluster) %>% count(Site)

dRep_count_5kb <- dRep_clusters_noSingle_5kb %>% group_by(secondary_cluster) %>% count(Site)

#see if any cluster now occurs twice
dRep_count <-  dRep_count %>% group_by(secondary_cluster) %>% filter(n()>1)

dRep_count_5kb <-  dRep_count_5kb %>% group_by(secondary_cluster) %>% filter(n()>1)

######################### dRep get counts of clusters across vents #######################
#reget input
dRep_clusters <- read.csv2(file = "Input/Cdb.csv", sep = ",")
dRep_clusters <- dRep_clusters %>% separate(genome, c("Site", NA), sep= "(?=_NODE|_scaffold|_vRhyme|_k95)")

dRep_clusters_5kb <- read.csv2(file = "Input/Cdb_VentPlume_5kb_90ani.csv", sep = ",")
dRep_clusters_5kb <- dRep_clusters_5kb %>% separate(genome, c("Site", NA), sep= "(?=_NODE|_scaffold|_vRhyme|_k95)")

#dRep drop singletons
dRep_clusters_noSingle <- dRep_clusters %>%
  add_count(secondary_cluster) %>%
  filter(n!=1) %>%
  select(-n)

dRep_clusters_noSingle_5kb <- dRep_clusters_5kb %>%
  add_count(secondary_cluster) %>% 
  filter(n!=1) %>%
  select(-n)

#dRep drop Pseudomonas viruses cluster because contamination
dRep_clusters_noSingle <- dRep_clusters_noSingle %>% filter(secondary_cluster!="1047_1")

dRep_clusters_noSingle_5kb <- dRep_clusters_noSingle_5kb %>% filter(secondary_cluster!="347_1")

#Changing names to general site to see clusters from different sites (not same vent field)
dRep_clusters_noSingle$Site <- gsub(".*Lau_Basin.*","Lau_Basin",dRep_clusters_noSingle$Site) #the placement of the periods is crucial for replacing whole string
dRep_clusters_noSingle$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",dRep_clusters_noSingle$Site)
dRep_clusters_noSingle$Site <- gsub(".*ELSC.*","Lau_Basin",dRep_clusters_noSingle$Site)
dRep_clusters_noSingle$Site <- gsub(".*Guaymas.*","Guaymas_Basin",dRep_clusters_noSingle$Site)
dRep_clusters_noSingle$Site <- gsub(".*Brothers.*","Brothers_Volcano",dRep_clusters_noSingle$Site)
dRep_clusters_noSingle$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",dRep_clusters_noSingle$Site)
dRep_clusters_noSingle$Site <- gsub(".*EPR.*","EPR",dRep_clusters_noSingle$Site)
dRep_clusters_noSingle$Site <- gsub(".*Axial.*","Axial_Seamount",dRep_clusters_noSingle$Site) #the placement of the periods is crucial for replacing whole string


dRep_clusters_noSingle_5kb$Site <- gsub(".*Lau_Basin.*","Lau_Basin",dRep_clusters_noSingle_5kb$Site) #the placement of the periods is crucial for replacing whole string
dRep_clusters_noSingle_5kb$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",dRep_clusters_noSingle_5kb$Site)
dRep_clusters_noSingle_5kb$Site <- gsub(".*ELSC.*","Lau_Basin",dRep_clusters_noSingle_5kb$Site)
dRep_clusters_noSingle_5kb$Site <- gsub(".*Guaymas.*","Guaymas_Basin",dRep_clusters_noSingle_5kb$Site)
dRep_clusters_noSingle_5kb$Site <- gsub(".*Brothers.*","Brothers_Volcano",dRep_clusters_noSingle_5kb$Site)
dRep_clusters_noSingle_5kb$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",dRep_clusters_noSingle_5kb$Site)
dRep_clusters_noSingle_5kb$Site <- gsub(".*EPR.*","EPR",dRep_clusters_noSingle_5kb$Site)
dRep_clusters_noSingle_5kb$Site <- gsub(".*Axial.*","Axial_Seamount",dRep_clusters_noSingle_5kb$Site) #the placement of the periods is crucial for replacing whole string

#count reps per cluster
dRep_clusters_noSingle_plot <- dRep_clusters_noSingle %>%
  group_by(secondary_cluster, Site) %>% 
  mutate(n = 1) %>%
  summarise(value=sum(as.numeric(n))) %>%
  group_by(secondary_cluster) %>%
  filter(n()>=2)

dRep_clusters_noSingle_5kb_plot <- dRep_clusters_noSingle_5kb %>%
  group_by(secondary_cluster, Site) %>% 
  mutate(n = 1) %>%
  summarise(value=sum(as.numeric(n))) %>%
  group_by(secondary_cluster) %>%
  filter(n()>=2)

#plot
dev.off()
plot <- dRep_clusters_noSingle_plot %>%
  ggplot(aes(x = secondary_cluster, y = as.numeric(value), fill = Site)) + 
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Cluster", y = "Viral genomes") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 8),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
        #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  scale_fill_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1",
                              "#8c510a", "#000000")) +
  #ggtitle("dRep Clusters 1kb, 95% ANI") +
  coord_flip()
plot

ggsave(plot, filename = "Output/dRep_clusters_1kb_distant_sites.png", dpi = 500, height = 5, width = 6)
ggsave(plot, filename = "Output/dRep_clusters_5kb_90ani_distant_sites.png", dpi = 500)

write.table(dRep_clusters_noSingle_plot, 
            file = "Output/dRep_clusters_1kb_95ani_sharedVents.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

############################################ venn diagram ####################################################

library()


############################################ mcl ####################################################

mcl_clusters <- read.csv2(file = "Input/vUnbinned_vMAGs_PlumeVents.clusters.csv", sep = ",", header = FALSE)

#mcl 1kb clusters with 50AF
mcl_clusters <- read.csv2(file = "Input/out_vUnbinned_vMAGs_PlumeVents_50AF.clusters.csv", sep = ",", header = FALSE)


#add id number to rows
mcl_clusters <- mcl_clusters %>% mutate(id = row_number())
#melt by id
mcl_clusters <- melt(mcl_clusters, id.vars = "id")
#drop variable column
mcl_clusters <- select(mcl_clusters, -variable)
#drop blanks in value column
mcl_clusters <- mcl_clusters %>%
  na_if("") %>%
  na.omit()
#rename column
mcl_clusters <- rename(mcl_clusters, "Site" = "value")

#write the mcl table in this format
# write.table(mcl_clusters,
#             file = "Output/mcl_formatted_table_VentPlumeViruses.tsv", sep = "\t", quote = FALSE,
#             row.names = FALSE, col.names = TRUE)

#replace strings with Vent or Plume
mcl_clusters$Site <- gsub(".*Lau_Basin.*","Plume",mcl_clusters$Site) #the placement of the periods is crucial for replacing whole string
mcl_clusters$Site <- gsub(".*Cayman.*","Plume",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Guaymas_Basin.*","Plume",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Axial.*","Plume",mcl_clusters$Site)
#separate the plume and vent dataframes to make life easier
mcl_clusters_Plume <- mcl_clusters %>% filter(Site == "Plume")
#remove Plume from mcl clusters so can make all names left Vent
mcl_clusters <- mcl_clusters %>% filter(Site != "Plume")
mcl_clusters$Site <- "Vent"
#Now put them back together
mcl_clusters <- rbind(mcl_clusters, mcl_clusters_Plume)

#Add counts to filter for clusters that occur more than once (not a singleton)
mcl_clusters_noSingle <- mcl_clusters %>%
  add_count(id) %>% 
  filter(n!=1) %>%
  select(-n)
#number of clusters with more than 1 rep according to dRep
length(unique(mcl_clusters_noSingle$id))
#2797
#1,020 with 50AF

#number of clusters with 1+ rep from Plume and Vent:
table(mcl_clusters_noSingle$Site)
#699 plume and 1591 vent with 50AF

#group by secondary cluster, count occurrences of Site
mcl_count <- mcl_clusters_noSingle %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
mcl_count <-  mcl_count %>% group_by(id) %>% filter(n()>1)

######################### mcl get counts of clusters across vents #######################

mcl_clusters <- read.csv2(file = "Input/vUnbinned_vMAGs_PlumeVents.clusters.csv", sep = ",", header = FALSE)
#mcl 1kb clusters with 50AF
mcl_clusters <- read.csv2(file = "Input/out_vUnbinned_vMAGs_PlumeVents_50AF.clusters.csv", sep = ",", header = FALSE)

#add id number to rows
mcl_clusters <- mcl_clusters %>% mutate(id = row_number())
#melt by id
mcl_clusters <- melt(mcl_clusters, id.vars = "id")
#drop variable column
mcl_clusters <- select(mcl_clusters, -variable)
#drop blanks in value column
mcl_clusters <- mcl_clusters %>%
  na_if("") %>%
  na.omit()
#rename column
mcl_clusters <- rename(mcl_clusters, "Site" = "value")
#remove Pseudomonas contam cluster
mcl_clusters <- mcl_clusters %>% filter(id!="1")

#replace strings with 7 general site names
mcl_clusters$Site <- gsub(".*Lau_Basin.*","Lau_Basin",mcl_clusters$Site) #the placement of the periods is crucial for replacing whole string
mcl_clusters$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*ELSC.*","Lau_Basin",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Guaymas.*","Guaymas_Basin",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Brothers.*","Brothers_Volcano",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*EPR.*","EPR",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Axial.*","Axial_Seamount",mcl_clusters$Site) #the placement of the periods is crucial for replacing whole string

#Add counts to filter for clusters that occur more than once (not a singleton)
mcl_clusters_noSingle <- mcl_clusters %>%
  add_count(id) %>% 
  filter(n!=1) %>%
  select(-n)

#number of clusters with 1+ rep from all sites:
table(mcl_clusters_noSingle$Site)

#group by secondary cluster, count occurrences of Site
mcl_count <- mcl_clusters_noSingle %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
mcl_count <-  mcl_count %>% group_by(id) %>% filter(n()>1)

############################ plot mcl across sites ###########################

#dRep drop Pseudomonas viruses cluster because contamination
#mcl_count <- mcl_count %>% filter(id!="3") 

#see if any cluster now occurs twice
mcl_count <-  mcl_count %>% group_by(id) %>% filter(n()>1)
mcl_count_3plus <-  mcl_count %>% group_by(id) %>% filter(n()>2)
#Remove BV + Lau overlap bc going to plot separately
mcl_count_test <- mcl_count %>%
  group_by(id) %>% 
  filter(Site!="Brothers_Volcano|Lau_Basin")
  
#plot
dev.off()
plot <- mcl_count %>%
  ggplot(aes(x = as.character(id), y = as.numeric(n), fill = Site)) + 
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Cluster", y = "Viral genomes") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 8),
        axis.text.y = element_text(size = 8),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  scale_fill_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1",
                             "#8c510a", "#000000")) +
  #ggtitle("dRep Clusters 1kb, 95% ANI") +
  coord_flip()
plot

ggsave(plot, filename = "Output/mcl_clusters_1kb_distant_sites.png", dpi = 500, height = 5, width = 6)
ggsave(plot, filename = "Output/mcl_clusters_5kb_90ani_distant_sites.png", dpi = 500)

############################ see sulfur virus relatedness ###########################

sulfur_viruses <- read.csv(file = "../iPHoP/output/iphop_VentPlume_sulfur_50comp.csv")
sulfur_viruses <- sulfur_viruses %>% select(Virus)
sulfur_viruses <- unique(sulfur_viruses)

#mcl 1kb clusters with 50AF
mcl_clusters <- read.csv2(file = "Input/out_vUnbinned_vMAGs_PlumeVents_50AF.clusters.csv", sep = ",", header = FALSE)

#add id number to rows
mcl_clusters <- mcl_clusters %>% mutate(id = row_number())
#melt by id
mcl_clusters <- melt(mcl_clusters, id.vars = "id")
#drop variable column
mcl_clusters <- select(mcl_clusters, -variable)
#drop blanks in value column
mcl_clusters <- mcl_clusters %>%
  na_if("") %>%
  na.omit()
#rename column
mcl_clusters <- rename(mcl_clusters, "Virus" = "value")
#remove Pseudomonas contam cluster
mcl_clusters <- mcl_clusters %>% filter(id!="1")

#fix vRhyme virus genome names
#separate the vrhyme and non-vrhyme to make life easier
mcl_clusters_vrhyme <- mcl_clusters %>% filter(grepl("vRhyme", Virus))
mcl_clusters_vrhyme <- mcl_clusters_vrhyme %>% separate(Virus, c("Virus", NA), 
                                                        sep= "(?=_scaffold|_NODE|_k95)")
mcl_clusters_vrhyme <- mcl_clusters_vrhyme %>% separate(Virus, c("Virus", "Site"), 
                                                        sep= "(__)")
mcl_clusters_vrhyme <- mcl_clusters_vrhyme %>% 
  mutate(Virus = str_replace(Virus, "vRhyme_", "vRhyme_bin_"))
mcl_clusters_vrhyme$Virus <- paste0(mcl_clusters_vrhyme$Site, "_", mcl_clusters_vrhyme$Virus)
mcl_clusters_vrhyme <- mcl_clusters_vrhyme %>% select(-Site)
#put the data frames back together again
mcl_clusters <- mcl_clusters %>% filter(!grepl("vRhyme", Virus))
mcl_clusters <- rbind(mcl_clusters_vrhyme, mcl_clusters)

#filter for only clusters containing sulfur cycling viruses
mcl_clusters_sulfur <- mcl_clusters[mcl_clusters$Virus %in% sulfur_viruses$Virus,] #subset table using list of names
#get list of IDs
mcl_clusters_sulfur <- mcl_clusters_sulfur %>% select(id)
mcl_clusters_sulfur <- unique(mcl_clusters_sulfur)
#27 unique clusters with 50comp sulfur cycling viruses

mcl_clusters_sulfur <- mcl_clusters[mcl_clusters$id %in% mcl_clusters_sulfur$id,]

####################### add average ANI value to mcl sulfur cluster table ###########################
skani_ANI <- read.delim2(file = "Input/skani_ANI_vUnbinned_vMAGs_PlumeVents_50AF_processed.tsv", header = FALSE)
#fix vRhyme names
skani_ANI_vrhyme <- skani_ANI %>% 
  


#Add counts to filter for clusters that occur more than once (not a singleton)
mcl_clusters_noSingle <- mcl_clusters %>%
  add_count(id) %>% 
  filter(n!=1) %>%
  select(-n)

#number of clusters with 1+ rep from all sites:
table(mcl_clusters_noSingle$Site)

#group by secondary cluster, count occurrences of Site
mcl_count <- mcl_clusters_noSingle %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
mcl_count <-  mcl_count %>% group_by(id) %>% filter(n()>1)









