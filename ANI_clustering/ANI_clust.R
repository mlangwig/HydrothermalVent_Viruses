#Script for creating skani processed file, parsing mcl cluster file, calculating average ANI
#per cluster, and mapping virus metadata onto clusters

library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)
library(reshape2)
library(textshape)
library(magrittr)

######################################### major inputs ################################################

skani_ani<-read.delim2(file = "Input/skani_v2/skani2_ANI_VentPlume_200m30cm_3kb.tsv")
mcl_clusters <- read.delim2(file = "Input/skani_v2/dereplicated_virus.clusters", sep = "\t", header = FALSE)

#################################### skani preprocessed ################################################
# After running skani, create the processed file here:

#make number columns numeric
skani_ani <- transform(skani_ani,
                          ANI = as.numeric(ANI),
                          Align_fraction_ref = as.numeric(Align_fraction_ref),
                          Align_fraction_query = as.numeric(Align_fraction_query))
#filter for >50AF
skani_ani<-filter(skani_ani, Align_fraction_ref >= 50.0)
skani_ani<-filter(skani_ani, Align_fraction_query >= 50.0)

#remove .fasta
skani_ani$Ref_file <- gsub(".fasta","",skani_ani$Ref_file)
skani_ani$Query_file <- gsub(".fasta","",skani_ani$Query_file)
#remove 3kb_vMAGs
# skani_ani$Ref_file <- gsub("vMAGs_3kb_50AF/","",skani_ani$Ref_file)
# skani_ani$Query_file <- gsub("vMAGs_3kb_50AF/","",skani_ani$Query_file)
skani_ani$Ref_file <- gsub("3kb_vMAGs/","",skani_ani$Ref_file)
skani_ani$Query_file <- gsub("3kb_vMAGs/","",skani_ani$Query_file)

###################### create file for mcl clustering ######################
ani <- skani_ani %>% 
  select(-Ref_file, -Query_file)

#only keep lowest AF
ani$align_frac = ifelse(ani$Align_fraction_ref < ani$Align_fraction_query,
                        ani$Align_fraction_ref, ani$Align_fraction_query)

#normalize ANI by lowest AF
ani$ANI_norm = ani$ANI*ani$align_frac/100^2

#remove extraneous columns
ani <- ani %>%
  select(-ANI, -Align_fraction_ref, -Align_fraction_query, -align_frac)

# remove  ANI below 70
ani<-filter(ani, ANI_norm >= 0.70)

# File "ani" and the following table output are what was used as input for 
#mcl clustering - these are the viruses that are clustered

# write.table(ani, file = "Input/skani_v2/skani2_ANI_VentPlume_200m30cm_3kb_preprocessed.tsv",
#             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


######################################### mcl ################################################

#add id number to rows
mcl_clusters <- mcl_clusters %>% mutate(id = row_number())
#melt by id
mcl_clusters <- melt(mcl_clusters, id.vars = "id")
#drop variable column
mcl_clusters <- select(mcl_clusters, -variable)
#remove NAs introduced after melting
mcl_clusters <- mcl_clusters %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  na.omit()
#rename column
mcl_clusters <- rename(mcl_clusters, "Site" = "value")

#remove singleton clusters
mcl_clusters <- mcl_clusters %>%
  group_by(id) %>%
  mutate(count=n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  filter(id!="1") %>% #remove Pseudomonas virus cluster
  ungroup()

#write the mcl table in this format
# write.table(mcl_clusters,
#             file = "Output/mcl_formatted_table_VentPlumeViruses.tsv", sep = "\t", quote = FALSE,
#             row.names = FALSE, col.names = TRUE)

################################ calculate average ANI per cluster ##########################################

#remove rows where comparing self to self
ani = subset(ani, ani$Ref_name != ani$Query_name)

#remove duplicate rows when comparing same thing but in different order
ani <- ani %>% 
  mutate(nv1 = paste0(Ref_name, ANI_norm),
         nv2 = paste0(Query_name, ANI_norm)) %>% 
  unique_pairs("nv1", "nv2") %>% 
  select(-nv1, -nv2)

#map cluster number to ANI table
#CREATE ani_perClust_filter_map and map ID to Ref_name
ani <- mcl_clusters %>%
  dplyr::select("id", "Site") %>%
  left_join(ani, by = c("Site" = "Ref_name")) %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  na.omit() %>%
  rename("Ref_name" = "Site")

#temporarily change col name for mapping
mcl_clusters <- mcl_clusters %>%
  rename("id_Q" = "id")

#map onto Query now instead of Ref
ani <- mcl_clusters %>%
  dplyr::select("id_Q", "Site") %>%
  left_join(ani, by = c("Site" = "Query_name")) %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  na.omit() %>%
  rename("Query_name" = "Site")
#change col name back
mcl_clusters <- mcl_clusters %>%
  rename("id" = "id_Q")

#filter out IDs that aren't the same
ani <- subset(ani, id_Q == id)

#get average ANI per cluster
ani <- ani %>% 
  group_by(id) %>% 
  mutate(ANI_mean = mean(ANI_norm))

################################ map virus metadata to clusters ##########################################

#melt the ani file so can do vlookups
ani_long <- ani %>% select(c('id', 'ANI_mean', 'Ref_name', 'Query_name'))
ani_long <- melt(ani_long, id = c('id', 'ANI_mean')) %>%
  select(-'variable') %>%
  unique()

##create needed input to map metadata to clusters
#change one column name so they match
vUnbinned_VP_master <- vUnbinned_VP_master %>% rename(vMAG = scaffold)
#get columns of interest
vMAG <- vMAG_VP_master %>% select(c('vMAG', 'type', 'contig_length',
                                    'checkv_quality', 'provirus',
                                    'completeness', 'contamination'))

vUnbinned <- vUnbinned_VP_master %>% select(c('vMAG', 'type', 'contig_length',
                                              'checkv_quality', 'provirus',
                                              'completeness', 'contamination'))
#bind/concatenate
allVirus_master <- allVirus_master %>% rbind(vMAG, vUnbinned) %>%
  unique()

#####change vMAG names from first scaffold name to file name for mapping

#separate the vrhyme and non-vrhyme to make life easier
ani_long_vrhyme <- ani_long %>% 
  filter(grepl("vRhyme", value)) %>%
  rename(value, "Virus"= "value" )

ani_long_vrhyme <- ani_long_vrhyme %>% separate(Virus, c("Virus", NA), 
                                                        sep= "(?=_scaffold|_NODE|_k95)")
ani_long_vrhyme <- ani_long_vrhyme %>% separate(Virus, c("Virus", "Site"), 
                                                        sep= "(__)")
ani_long_vrhyme <- ani_long_vrhyme %>% 
  mutate(Virus = str_replace(Virus, "vRhyme_", "vRhyme_bin_"))
ani_long_vrhyme$Virus <- paste0(ani_long_vrhyme$Site, "_", ani_long_vrhyme$Virus)
ani_long_vrhyme <- ani_long_vrhyme %>% select(-Site)

#put the data frames back together again
ani_long <- ani_long %>% filter(!grepl("vRhyme", value)) %>%
  rename("Virus" = "value")
ani_long <- rbind(ani_long_vrhyme, ani_long)

#vlookup
ani_long_metadata <- allVirus_master %>%
  dplyr::select('vMAG', 'type', 'contig_length',
                'checkv_quality', 'provirus',
                'completeness', 'contamination') %>%
  right_join(ani_long, by = c("vMAG" = "Virus"))


################################ see clusters that have vent and plume ##########################################

#replace strings with Vent or Plume
ani_long_metadata$Site <- gsub(".*Lau_Basin.*","Plume",ani_long_metadata$vMAG) #the placement of the periods is crucial for replacing whole string
ani_long_metadata$Site <- gsub(".*Cayman.*","Plume",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*Guaymas_Basin.*","Plume",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*Axial.*","Plume",ani_long_metadata$Site)
#separate the plume and vent dataframes to make life easier
temp_plume <- ani_long_metadata %>% filter(Site == "Plume")
#remove Plume from mcl clusters so can make all names left Vent
ani_long_metadata <- ani_long_metadata %>% filter(Site != "Plume")
ani_long_metadata$Site <- "Vent"
#Now put them back together
ani_long_metadata <- rbind(ani_long_metadata, temp_plume)

#Add counts to filter for clusters that occur more than once (not a singleton)
# temp_count <- ani_long_metadata %>%
#   add_count(id) %>% 
#   filter(n!=1) %>%
#   select(-n)

#number of viruses from Plume and Vent:
table(ani_long_metadata$Site)
#406 plume and 1030 vent with skani parameters, 3kb, and 50AF

#group by cluster, count occurrences of Site
temp_count <- ani_long_metadata %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
temp_count <-  temp_count %>% group_by(id) %>% filter(n()>1)
#no cluster contains plume + vent

################################ see clusters that contain distinct sites ##########################################

#replace strings with 7 general site names
ani_long_metadata$Site <- gsub(".*Lau_Basin.*","Lau_Basin",ani_long_metadata$vMAG)
ani_long_metadata$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*ELSC.*","Lau_Basin",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*Guaymas.*","Guaymas_Basin",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*Brothers.*","Brothers_Volcano",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*EPR.*","EPR",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub(".*Axial.*","Axial_Seamount",ani_long_metadata$Site) 

#number of clusters with 1+ rep from all sites:
table(ani_long_metadata$Site)

#count occurrences of Site
temp_count <- ani_long_metadata %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
temp_count <-  temp_count %>% group_by(id) %>% filter(n()>1)

############################ visualize counts across sites ###########################
temp_count$id <- as.character(temp_count$id)

#plot
dev.off()
plot <- temp_count %>%
  ggplot(aes(x = reorder(id, rev(sort(as.numeric(id)))), y = as.numeric(n), fill = Site)) + #reorder is a fun new trick! to sort the order for plotting without changing the str type
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

#ggsave(plot, filename = "Output/mcl_GeoDistinct_clusters.png", dpi = 500, height = 6, width = 6)







