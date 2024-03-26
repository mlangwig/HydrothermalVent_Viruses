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
library(stringi)

######################################### major inputs ################################################

ani<-read.delim2(file = "Input/skani_49962/skani_v2_ANI_VentVirus_200m30cm_3kb_49962_50AF_processed.tsv",
                 header = FALSE)
ani <- ani %>%
  rename("Ref_name" = "V1",
         "Query_name" = "V2",
         "ANI_norm" = "V3")
mcl_clusters <- read.delim2(file = "Input/skani_49962/dereplicated_virus_49962.clusters", sep = "\t", header = FALSE)
genomad <- read.delim2(file = "../VentVirus_Analysis/input/49962_final_VentViruses_Nlinked_taxonomy_parsed.tsv")
iphop <- read.csv(file = "../VentVirus_Analysis/input/Host_prediction_to_genus_m90_49962.csv")
metadata <- read.delim2(file = "../VentVirus_Analysis/output/master_table_VentPlumeViruses_simple.tsv")
metadata <- metadata %>%
  select(c(vMAG, type:f)) %>%
  unique()

#################################### skani preprocessed ################################################
# # After running skani, create the processed file here:
# 
# #make number columns numeric
# skani_ani <- transform(skani_ani,
#                           ANI = as.numeric(ANI),
#                           Align_fraction_ref = as.numeric(Align_fraction_ref),
#                           Align_fraction_query = as.numeric(Align_fraction_query))
# #filter for >50AF
# skani_ani<-filter(skani_ani, Align_fraction_ref >= 50.0)
# skani_ani<-filter(skani_ani, Align_fraction_query >= 50.0)
# 
# #remove .fasta
# skani_ani$Ref_file <- gsub(".fasta","",skani_ani$Ref_file)
# skani_ani$Query_file <- gsub(".fasta","",skani_ani$Query_file)
# #remove 3kb_vMAGs
# # skani_ani$Ref_file <- gsub("vMAGs_3kb_50AF/","",skani_ani$Ref_file)
# # skani_ani$Query_file <- gsub("vMAGs_3kb_50AF/","",skani_ani$Query_file)
# skani_ani$Ref_file <- gsub("3kb_vMAGs/","",skani_ani$Ref_file)
# skani_ani$Query_file <- gsub("3kb_vMAGs/","",skani_ani$Query_file)
# 
# ###################### create file for mcl clustering ######################

#THIS WAS COMPLETED ON THE SERVER USING AN R SCRIPT THAT DOES THE SAME EXACT AS THE FOLLOWING

# ani <- skani_ani %>%
#   select(-Ref_file, -Query_file)
# 
# #only keep lowest AF
# ani$align_frac = ifelse(ani$Align_fraction_ref < ani$Align_fraction_query,
#                         ani$Align_fraction_ref, ani$Align_fraction_query)
# 
# #normalize ANI by lowest AF
# ani$ANI_norm = ani$ANI*ani$align_frac/100^2
# 
# #remove extraneous columns
# ani <- ani %>%
#   select(-ANI, -Align_fraction_ref, -Align_fraction_query, -align_frac)
# 
# # remove  ANI below 70
# ani<-filter(ani, ANI_norm >= 0.70)


# File "ani" and the following table output are what was used as input for
#mcl clustering - these are the viruses that are clustered
# 
# # write.table(ani, file = "Input/skani_v2/skani2_ANI_VentPlume_200m30cm_3kb_12_31_50AF_preprocessed.tsv",
# #             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


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
  #filter(id!="1") %>% #remove Pseudomonas virus cluster
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
ani$ANI_norm <- as.numeric(ani$ANI_norm)
ani <- ani %>% 
  group_by(id) %>% 
  mutate(ANI_mean = mean(ANI_norm)) %>%
  ungroup()

################################ map virus metadata to clusters ##########################################

#melt the ani file so can do vlookups
ani_long <- ani %>% select(c('id', 'ANI_mean', 'Ref_name', 'Query_name'))
ani_long <- melt(ani_long, id = c('id', 'ANI_mean')) %>%
  select(-'variable') %>%
  unique()

#check how many clusters are composed of just 2 viruses
ani_2clusts <- ani_long %>%
  group_by(id) %>%
  filter(n() == 2)
length(unique(ani_2clusts$id))
#687 are pairwise

#####change vMAG names from first scaffold name to file name for mapping
#separate the vrhyme and non-vrhyme to make life easier
ani_long_vrhyme <- ani_long %>% 
  filter(grepl("vRhyme", value)) %>%
  rename(value, "Virus"= "value" )
#151 vMAGs in ani clusts

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
ani_long_metadata <- metadata %>%
  right_join(ani_long, by = c("vMAG" = "Virus"))

#NOTE THAT MAPPING IPHOP DUPLICATES SOME VIRUSES BC SOME VIRUSES HAVE 1+ HOST PREDICTION
#add iphop-predicted hosts
ani_long_metadata <- iphop %>%
  dplyr::select('Virus', 'Host.genus', 'List.of.methods') %>%
  right_join(ani_long_metadata, by = c("Virus" = "vMAG")) %>%
  rename("vMAG" = "Virus")

#add site column for ani_long_metadata
ani_long_metadata$Site <- ani_long_metadata$vMAG 
ani_long_metadata <- ani_long_metadata %>%
  separate(Site, c("Site", NA), sep= "_NODE|_scaffold|_vRhyme|_k95")

#clean up site names
#faster replace all the naming patterns
ani_long_metadata$Site <- stri_replace_all_regex(ani_long_metadata$Site,
                                                pattern=c("_A[0-9]",
                                                          "_T[0-9][0-9]", "_T[0-9]", "_S0[0-9][0-9]",
                                                          "_S1[0-9][0-9]", "_[0-9][0-9][0-9]-[0-9][0-9][0-9]",
                                                          "-38[0-9]"),
                                                replacement='',
                                                vectorize=FALSE)
ani_long_metadata$Site <- gsub("*_M1[0-9]","",ani_long_metadata$Site)
ani_long_metadata$Site <- gsub("*_M[0-9]","",ani_long_metadata$Site)

########################## see clusters that have both deposit and plume viruses #######################################

#replace strings with Deposit or Plume
ani_long_metadata$Site_type <- gsub(".*Lau_Basin.*","Plume",ani_long_metadata$Site) #the placement of the periods is crucial for replacing whole string
ani_long_metadata$Site_type <- gsub(".*Cayman.*","Plume",ani_long_metadata$Site_type)
ani_long_metadata$Site_type <- gsub(".*Guaymas_Basin.*","Plume",ani_long_metadata$Site_type)
ani_long_metadata$Site_type <- gsub(".*Axial.*","Plume",ani_long_metadata$Site_type)
#separate the plume and vent dataframes to make life easier
temp_plume <- ani_long_metadata %>% filter(Site_type == "Plume")
#remove Plume from mcl clusters so can make all names left Vent
ani_long_metadata <- ani_long_metadata %>% filter(Site_type != "Plume")
ani_long_metadata$Site_type <- "Deposit"
#Now put them back together
ani_long_metadata <- rbind(ani_long_metadata, temp_plume)

#number of viruses from Plume and Deposit:
table(ani_long_metadata$Site_type)
#744 plume and 1259 deposit

#group by cluster, count occurrences of Site
temp_count <- ani_long_metadata %>% group_by(id) %>% count(Site_type)
#see if any cluster now occurs twice
temp_count <-  temp_count %>% group_by(id) %>% filter(n()>1)
# --> no cluster contains plume + deposit <--

########################## see clusters that contain geographically distinct sites ##########################################

#replace strings with 7 general site names
ani_long_metadata$Site_gen <- gsub(".*Lau_Basin.*","Lau_Basin",ani_long_metadata$Site)
ani_long_metadata$Site_gen <- gsub(".*Cayman.*","Mid_Cayman_Rise",ani_long_metadata$Site_gen)
ani_long_metadata$Site_gen <- gsub(".*ELSC.*","Lau_Basin",ani_long_metadata$Site_gen)
ani_long_metadata$Site_gen <- gsub(".*Guaymas.*","Guaymas_Basin",ani_long_metadata$Site_gen)
ani_long_metadata$Site_gen <- gsub(".*Brothers.*","Brothers_Volcano",ani_long_metadata$Site_gen)
ani_long_metadata$Site_gen <- gsub(".*MAR.*","Mid_Atlantic_Ridge",ani_long_metadata$Site_gen)
ani_long_metadata$Site_gen <- gsub(".*EPR.*","EPR",ani_long_metadata$Site_gen)
ani_long_metadata$Site_gen <- gsub(".*Axial.*","Axial_Seamount",ani_long_metadata$Site_gen) 

#count occurrences of Site
temp_count <- ani_long_metadata %>% group_by(id) %>% count(Site_gen)
#see if any cluster now occurs twice
temp_count <-  temp_count %>% group_by(id) %>% filter(n()>1)
length(unique(temp_count$id))
#65 clusters contain cross site viruses

#see these clusters from ani_long_metadata
ids <- as.integer(unique(temp_count$id))
ani_long_meta_gd <- ani_long_metadata %>% filter(ani_long_metadata$id %in% ids)
length(unique(ani_long_meta_gd$vMAG))
#152 viruses in clusters from geo distinct locations

#remove undetermined from ani_long_meta_gd
# test <- ani_long_meta_gd %>% filter(checkv_quality != "Not-determined")
# ani_long_meta_gd <- ani_long_meta_gd %>% filter(taxonomy != "NA") %>%
#   filter(taxonomy != "Unclassified")

write.table(ani_long_meta_gd, file = "Output/ani_metadata_GeoDistinct.tsv", quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

############################ visualize counts across sites ###########################
temp_count$id <- as.character(temp_count$id)

#plot
dev.off()
plot <- temp_count %>%
  ggplot(aes(x = reorder(id, rev(sort(as.numeric(id)))), y = as.numeric(n), fill = Site_gen)) + #reorder is a fun new trick! to sort the order for plotting without changing the str type
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
  scale_y_continuous(expand = c(0, 0)) #+ #turn this on to make it look aligned with ticks
  # scale_fill_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1",
  #                            "#8c510a", "#000000")) +
  #ggtitle("dRep Clusters 1kb, 95% ANI") +
  #coord_flip()
plot

#ggsave(plot, filename = "Output/mcl_GeoDistinct_clusters.png", dpi = 500, height = 6, width = 6)


################################ intra-vent viral clusters ##########################################

#remove clusters that are geo distinct
id_rem <- unique(ani_long_meta_gd$id)
"%ni%" <- Negate("%in%")
ani_long_meta_iv <- ani_long_metadata %>% filter(ani_long_metadata$id %ni% id_rem)

## Remove clusters that only have ANI within a site e.g. clust 7 all NWCB
#count occurrences of Site
temp_count <- ani_long_meta_iv %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
temp_count <-  temp_count %>% group_by(id) %>% filter(n()>1)
iv_ids <- unique(temp_count$id)

#create iv file of just the ones of interest (actually IV)
#THIS IF FINAL IV FILE
#note will have issues if you reordered it in CircosPlot.R
ani_long_meta_iv <- ani_long_meta_iv %>% filter(ani_long_meta_iv$id %in% iv_ids)
length(unique(ani_long_meta_iv$id))
# 462 clusters with intra vent related viruses


#HERE FOR FILT VERSION to do counts
# #remove viruses that are Unclassified or NA
ani_long_meta_iv_filt <- ani_long_meta_iv %>%
  #filter(taxonomy != "NA") %>%
  #filter(taxonomy != "Unclassified") %>%
  #filter(checkv_quality != "Low-quality") %>%
  #filter(checkv_quality != "Not-determined") %>%
  filter(grepl("Axial|EPR", Site)) #%>%
  #filter(ANI_mean > 0.85)
  # filter(grepl("Lau_Basin_Kilo_Moana", Site)) %>%
  # filter(grepl("Lau_Basin_Abe", Site))
  

#write the table
# write.table(ani_long_meta_iv, file = "Output/ani_metadata_IntraVent.tsv", quote = FALSE, row.names = FALSE,
#             col.names = TRUE, sep = "\t")

################################## calculate num of times combos of sites occur ################################################

plot_iv <- ani_long_meta_iv_filt %>%
  group_by(id) %>% 
  count(Site)

#count the occurrences of dif combos of sites grouped by cluster ID
plot_iv <- ani_long_meta_iv_filt %>%
  group_by(id) %>%
  mutate(Site = toString(Site)) %>% #toString to get them in the comma separated list 
  select(id, Site) %>%
  unique()
test <- table(plot_iv$Site)
plot_iv <- as.data.frame(test)

length(unique(plot_iv$Var1)) # can unique_pairs be used? --> table(unique_pairs(test$Site))

################################## plot id and site composition ################################################

#modify input for faceting
plot_iv$Site <- plot_iv$Var1 
plot_iv$Site <- as.character(plot_iv$Site)
plot_iv <- plot_iv %>% separate(Site, c("Site", NA), sep= "_")
#remove long names in Var1
plot_iv$Var1 <- gsub("Lau_Basin_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("Brothers_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("Axial_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("Cayman_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("ELSC_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("Guaymas_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("EPR_","", plot_iv$Var1)
plot_iv$Var1 <- gsub("_"," ", plot_iv$Var1)
#fix up site names for outer label
plot_iv$Site <- gsub("Lau","Lau Basin Plume", plot_iv$Site)
plot_iv$Site <- gsub("ELSC","Lau Basin Vent", plot_iv$Site)
plot_iv$Site <- gsub("Brothers","Brothers Volcano", plot_iv$Site)
plot_iv$Site <- gsub("Guaymas","Guaymas Basin", plot_iv$Site)
plot_iv$Site <- gsub("Axial","Axial Seamount", plot_iv$Site)
plot_iv$Site <- gsub("EPR","East Pacific Rise", plot_iv$Site)
plot_iv$Site <- gsub("Cayman","Mid-Cayman Rise", plot_iv$Site)

#factor outer label sites
plot_iv$Site_f = factor(plot_iv$Site, levels=c('Axial Seamount','Mid-Cayman Rise','Lau Basin Plume', 'Lau Basin Vent',
                                               'East Pacific Rise', 'Guaymas Basin', 'Brothers Volcano'))
#add ordering column for plotting in order within each facet
# plot_iv <- plot_iv %>%
#   group_by(Site) %>%
#   arrange(desc(Freq), .by_group = T)
# 
# test <- plot_iv %>%
#   mutate(Var1 = reorder_within(as.character(Var1), desc(as.numeric(Freq)), Site_f)) %>%
#   separate(Var1, c("Var1", NA), sep= "__")

#install.packages("tidytext")
library(tidytext)

#thank you https://juliasilge.com/blog/reorder-within/ for your helpful example

colors <- c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478",
                      "Guaymas Basin" = "#28827A", "Lau Basin Plume" = "#3F78C1",
                      "Lau Basin Vent" = "#72a0db", "Mid-Cayman Rise" = "#000000", "East Pacific Rise" = "#CE9A28")

dev.off()
plot <- plot_iv %>%
  ggplot(aes(x = as.numeric(Freq), y = reorder_within(as.character(Var1), desc(as.numeric(Freq)), Site_f), fill = Site_f), alpha = Var1) + #reorder within to get desc nums in faceted plot
  geom_bar(stat = "identity") +
  scale_alpha_continuous(range= c(0.1,1)) +
  scale_fill_manual(values = colors, name = "Site") +
  labs(y = "Cluster", x = "Viral clusters") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 1, vjust = .5, size = 8),
        axis.text.x = element_text(size = 8),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.text.y = element_text(angle = 0)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_x_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  scale_y_reordered(limits=rev) +
  # scale_fill_brewer(palette = "Set1",
  #                   name = "Site") +
  # scale_alpha_discrete(range= c(0.4,1)) +
  facet_grid(Site_f ~ ., scales = "free_y", space = "free_y") #+
  #ggtitle("dRep Clusters 1kb, 95% ANI") +
  #coord_flip()
plot

################################## plot ani vs virus completeness ################################################

skani_ani<-read.delim2(file = "Input/skani_v2/skani2_ANI_VentPlume_200m30cm_3kb.tsv")

#make number columns numeric
skani_ani <- transform(skani_ani,
                       ANI = as.numeric(ANI),
                       Align_fraction_ref = as.numeric(Align_fraction_ref),
                       Align_fraction_query = as.numeric(Align_fraction_query))

#remove .fasta
skani_ani$Ref_file <- gsub(".fasta","",skani_ani$Ref_file)
skani_ani$Query_file <- gsub(".fasta","",skani_ani$Query_file)
#remove 3kb_vMAGs
skani_ani$Ref_file <- gsub("3kb_vMAGs/","",skani_ani$Ref_file)
skani_ani$Query_file <- gsub("3kb_vMAGs/","",skani_ani$Query_file)

#create ani file
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

#remove rows where comparing self to self
ani = subset(ani, ani$Ref_name != ani$Query_name)

#remove duplicate rows when comparing same thing but in different order
ani <- ani %>% 
  mutate(nv1 = paste0(Ref_name, ANI_norm),
         nv2 = paste0(Query_name, ANI_norm)) %>% 
  unique_pairs("nv1", "nv2") %>% 
  select(-nv1, -nv2)

#melt the data
ani_long <- ani %>% pivot_longer(cols = c('Ref_name', 'Query_name')) %>%
  select(-name) %>%
  group_by(value) %>% 
  mutate(ANI_mean = mean(ANI_norm)) %>%
  select(-ANI_norm) %>%
  unique()

###change vMAG names from first scaffold name to file name for mapping

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
ani_long_metadata_all <- allVirus_master %>%
  dplyr::select('vMAG', 'type', 'contig_length',
                'checkv_quality', 'provirus',
                'completeness', 'contamination') %>%
  right_join(ani_long, by = c("vMAG" = "Virus"))

##################### remove e coli contam viruses from allVirus_master ###########################
#remove names from the list
allVirus_master <- allVirus_master %>%
  filter(!str_detect(vMAG, remove.list))
#replace equals sign with underscore to get complete removal
allVirus_master_nu <- allVirus_master
allVirus_master_nu$vMAG <- gsub('=','_',allVirus_master_nu$vMAG)

#now replace again and should work for all 
allVirus_master_nu <- allVirus_master_nu %>%
  filter(!str_detect(vMAG, remove.list))

#only med-quality better virus master
allVirus_master_nu_hq <- allVirus_master_nu %>%
  filter(!grepl("Low-quality|Not-determined", checkv_quality))

##################################################

##################### remove e coli contam viruses from gensize ###########################

#remove names from the list
gensize_ne <- gensize %>%
  filter(!str_detect(file, remove.list))

sum(gensize_ne$sum_len >= 1000 & gensize_ne$sum_len <= 5000)

##################################################

##################### remove e coli contam viruses from genomad ###########################

#replace equals sign with underscore to get complete removal
genomad_ne <- genomad
genomad_ne$seq_name <- gsub('=','_',genomad_ne$seq_name)

#remove names from the list
genomad_ne <- genomad_ne %>%
  filter(!str_detect(seq_name, remove.list))

#get only Caudos
genomad_ne_caud <- genomad_ne %>%
  filter(grepl("Duplodnaviria", taxonomy)) #%>%
  filter(grepl("Caudoviricetes", taxonomy)) #Caud = 17,964

##################################################

#drop NAs to see if that fixes plotting probs
ani_long_metadata_all <- ani_long_metadata_all %>%
  filter_at(vars(completeness), all_vars(!is.na(.)))
#6755 viruses remaining from 8364 that had ANI averages


library(scales)
ani_long_metadata_all$completeness <- as.numeric(ani_long_metadata_all$completeness)

#### plot ANI and completeness
#plot
dev.off()
plot <- ani_long_metadata_all %>%
  ggplot(aes(x = as.numeric(completeness), y = as.numeric(ANI_mean))) + #reorder is a fun new trick! to sort the order for plotting without changing the str type 
  geom_point() +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Percent Completeness", y = "ANI") +
  #guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 8),
        #axis.text.y = element_text(size = 2),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  #scale_x_continuous(breaks = c(0, 100, by = 20)) #+ #turn this on to make it look aligned with ticks
  # scale_fill_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1",
  #                            "#8c510a", "#000000")) +
  #ggtitle("dRep Clusters 1kb, 95% ANI") +
  coord_flip() +
  geom_hline(yintercept = 0.70, linetype = 3)
plot


########################### unused
# library(RColorBrewer)
# n <- 61
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))


##create needed input to map metadata to clusters
#change one column name so they match
# 
# vUnbinned_VP_master <- vUnbinned_VP_master %>% rename(vMAG = scaffold)
# #get columns of interest
# vMAG <- vMAG_VP_master %>% select(c('vMAG', 'type', 'contig_length',
#                                     'checkv_quality', 'provirus',
#                                     'completeness', 'contamination'))
# 
# vUnbinned <- vUnbinned_VP_master %>% select(c('vMAG', 'type', 'contig_length',
#                                               'checkv_quality', 'provirus',
#                                               'completeness', 'contamination'))
# #bind/concatenate
# allVirus_master <- rbind(vMAG, vUnbinned) %>%
#   unique()
# #allVirus_master %>% 

# #add genomad taxonomy
# ani_long_metadata <- genomad %>%
#   dplyr::select('seq_name', 'taxonomy', 'n_genes', 'n_hallmarks') %>%
#   right_join(ani_long_metadata, by = c("seq_name" = "vMAG")) %>%
#   rename("vMAG" = "seq_name")

#Add counts to filter for clusters that occur more than once (not a singleton)
# temp_count <- ani_long_metadata %>%
#   add_count(id) %>% 
#   filter(n!=1) %>%
#   select(-n)

# #Create Site column and create site names
# ani_long_metadata$Site <- ani_long_metadata$vMAG
# ani_long_metadata <- ani_long_metadata %>% separate(Site, c("Site", NA), sep= "(?=_scaffold|_NODE|_k95|_vRhyme)")
# 
# #Remove some parts of site names
# ani_long_metadata$Site <- gsub("_A[0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("_M10","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("_T[0-9][0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("_T[0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("_M10","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*_S0[0-9][0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*_S1[0-9][0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*_M1[0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*_M[0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*_[0-9][0-9][0-9]-[0-9][0-9][0-9]","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*-380","",ani_long_metadata$Site)
# ani_long_metadata$Site <- gsub("*-384","",ani_long_metadata$Site)
# 
# #see these clusters from ani_long_metadata
# #ids <- as.integer(unique(temp_count$id))
