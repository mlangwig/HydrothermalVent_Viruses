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

#change vMAG names from first scaffold to file name for mapping


#vlookup
test <- allVirus_master %>%
  dplyr::select('vMAG', 'type', 'contig_length',
                'checkv_quality', 'provirus',
                'completeness', 'contamination') %>%
  right_join(ani_long, by = c("vMAG" = "value"))

