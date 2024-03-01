library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)
library(reshape2)
library(textshape)
library(magrittr)
library(readr)
library(data.table)

################################# get protein annotations  ####################################

########################### inputs
vibrant <- read_tsv(file = "Input/52_VIBRANT_annos_long_viruses.tsv", col_names = FALSE)
vibrant <- vibrant %>%
  rename("protein" = "X1") %>%
  rename("id" = "X2") %>%
  rename("evalue" = "X3") %>%
  rename("score" = "X4")
#this is the VIBRANT annotation data from all of the hmmtbl_parse.tsv files VIBRANT produces

#master table for annotation mapping
vib_annos <- read_tsv(file = "../VentVirus_Analysis/output/master_table_VentPlumeViruses_simple.tsv")
vib_annos <- vib_annos %>%
  select(c(KO, KO.name, Pfam, Pfam.name, VOG, VOG.name)) %>%
  unique()

kos <- vib_annos %>%
  select(KO, KO.name) %>%
  rename("id" = "KO",
         "anno" = "KO.name") %>%
  drop_na()
pfams <- vib_annos %>%
  select(Pfam, Pfam.name) %>%
  rename("id" = "Pfam",
         "anno" = "Pfam.name") %>%
  drop_na()
vogs <- vib_annos %>%
  select(VOG, VOG.name) %>%
  rename("id" = "VOG",
         "anno" = "VOG.name") %>%
  drop_na()
vib_annos <- rbind(kos,pfams,vogs)
vib_annos <- unique(vib_annos)

###########################

#filter for best supported annotation based on bit score and evalue
#check if something has same bit score and e value
test <- vibrant %>%
  group_by(protein) %>%
  filter(n() > 1, all(score == first(score), evalue == first(evalue)))
#one example of this, Lau_Basin_Abe_k95_1492896_flag_1_multi_24.9933_len_5352_6
vibrant_best <- vibrant %>%
  group_by(protein) %>%
  slice(which.max(score)) %>% #take highest bit score
  slice(which.min(evalue)) #%>% #take lowest e value
#this code kept the first anno for the Lau virus with two hits same scores, I am keeping it because
#this hit has a higher v-score (Pfam v-score =1.0, VOG v-score = 0.79)

#add annotation info to vibrant_best
vibrant_best_anno <- vib_annos %>%
  right_join(vibrant_best, by = c("id" = "id"))

######################################### mmseqs Input ################################################
mmseqs <- read.csv2(file = "Input/PlumeVent_mmseqs_clusters_49962viruses.tsv", sep = "\t", header = FALSE)
#595,416 proteins clustered
#note that output from mmseqs2 is in long format so repeat proteins in first column
#mmseqs does not output singletons

mmseqs <- mmseqs %>%
  rename("cluster.representative" = "V1") %>%
  rename("cluster.member" = "V2")
############################################ mmseqs ####################################################

#Filter for non-self match in cluster
mmseqs <- mmseqs %>%
  filter(cluster.representative != cluster.member)
#this should now be number of clusters
length(unique(mmseqs_dif$cluster.representative)) #74,940

#change format to long
mmseqs_long <- mmseqs %>%
  group_by(cluster.representative) %>%
  mutate(id = cur_group_id()) %>%
  pivot_longer(cols = c('cluster.representative', 'cluster.member')) %>%
  select(-name) %>%
  rename("genome" = "value") %>%
  ungroup() %>%
  group_by(id) %>%
  unique() %>%
  ungroup()

########### mmseqs file where you map metadata in additional cols
#column without the vRhyme name for mapping
#create new mmseqs file
mmseqs_long_meta <- mmseqs_long
#new col
mmseqs_long_meta$genome_vRhyme <- mmseqs_long_meta$genome
#separate col
mmseqs_long_meta <- mmseqs_long_meta %>% separate(genome_vRhyme, c(NA, "genome_vRhyme"), sep = "__")
#put non vRhyme names back in the column
mmseqs_long_meta <- mmseqs_long_meta %>%
  mutate(genome_vRhyme = if_else(is.na(genome_vRhyme), genome, genome_vRhyme))

#column with Plume/Deposit for parsing
mmseqs_long_meta$PD <- mmseqs_long_meta$genome_vRhyme
mmseqs_long_meta <- mmseqs_long_meta %>%
  mutate(PD = gsub(".*Lau_Basin.*","Plume", PD), # a with acute
         PD = gsub(".*Cayman.*","Plume", PD),
         PD = gsub(".*Guaymas_Basin.*","Plume", PD), # a with acute
         PD = gsub(".*Axial.*","Plume", PD) # a with acute
  )

#replace everything not plume with deposit
mmseqs_long_meta <- mmseqs_long_meta %>%
  mutate(PD = ifelse(grepl("Plume", PD), PD, "Deposit"))

#column with general site name for parsing
mmseqs_long_meta$Site <- mmseqs_long_meta$genome_vRhyme
mmseqs_long_meta <- mmseqs_long_meta %>%
  mutate(Site = gsub(".*Lau_Basin.*","Lau_Basin", Site),
         Site = gsub(".*Cayman.*","Cayman", Site),
         Site = gsub(".*ELSC.*","Lau_Basin", Site), #not distinguishing Laus here bc PD will get at that
         Site = gsub(".*Guaymas.*","Guaymas", Site),
         Site = gsub(".*Brothers.*","Brothers", Site),
         Site = gsub(".*MAR.*","MAR", Site), 
         Site = gsub(".*EPR.*","EPR", Site),
         Site = gsub(".*Axial.*","Axial", Site)
  )

#add VIBRANT annotations to mmseqs meta
vibrant_best_anno <- vibrant_best_anno %>%
  rename("anno_id" = "id")
mmseqs_long_meta_annos <- vibrant_best_anno %>%
  right_join(mmseqs_long_meta, by = c("protein" = "genome_vRhyme")) %>%
  drop_na()

# Find which clusters have proteins from Plume vs Deposit
mmseqs_PD_ids <- mmseqs_long_meta %>%
  select(id, PD) %>%
  unique() %>%
  group_by(id) %>% 
  filter(n()>1) %>%
  select(id) %>%
  unique()
length(unique(mmseqs_PD_ids$id))
# 152 clusters that have proteins from plumes and deposits
# 138 of these are also have proteins from geographically distinct vents

mmseqs_PD <- mmseqs_long_meta %>%
  filter(id %in% mmseqs_PD_ids$id)
# 773 proteins from the 152 clusters = proteins from plumes and deposits

# Find which clusters have proteins from Geographically Distinct locations
mmseqs_GD_ids <- mmseqs_long_meta %>%
  select(id, Site) %>%
  unique() %>%
  group_by(id) %>% 
  filter(n()>1) %>%
  select(id) %>%
  unique()
length(unique(mmseqs_GD_ids$id))
# 23,351 clusters that have proteins from plumes and deposits

mmseqs_GD <- mmseqs_long_meta %>%
  filter(id %in% mmseqs_GD_ids$id)
# 84,223 proteins from the 23k clusters = proteins from geographically distinct locations
#Not removing the duplicates in PD from GD because making different points

# Make a file of PD and GD together
mmseqs_PD_GD_ids <- rbind(mmseqs_PD_ids, mmseqs_GD_ids) #23,503 clusts
mmseqs_PD_GD_ids <- unique(mmseqs_PD_GD_ids) #23,365 clusts = 138 clusts same between PD and GD
mmseqs_PD_GD <- mmseqs_long_meta %>%
  filter(id %in% mmseqs_PD_GD_ids$id)

# write_delim(mmseqs_PD_GD, file = "Output/virus_proteins_mmseqsClusts_PD_GD.tsv", 
#             col_names = TRUE, delim = "\t")
# write_delim(mmseqs_GD, file = "Output/virus_proteins_mmseqsClusts_GD.tsv", 
#             col_names = TRUE, delim = "\t")
# write_delim(mmseqs_PD, file = "Output/virus_proteins_mmseqsClusts_PD.tsv", 
#             col_names = TRUE, delim = "\t")

#add the mmseqs annotation info to the PD and GD cluster file to see how many have annotations/what they are

mmseqs_PD_GD_annos <- mmseqs_long_meta_annos %>%
  dplyr::select(c("anno_id", "anno", "protein")) %>%
  right_join(mmseqs_PD_GD, by = c("protein" = "genome_vRhyme")) %>%
  drop_na() %>%
  rename("protein_no_vRhyme" = "protein",
         "protein" = "genome") %>%
  select(c(id, protein, protein_no_vRhyme,PD, Site, anno_id, anno))






################### unused

# test <- read.delim2("../../VirusGenomes/49962_faas/virus_prots_list_renamed.txt", header = FALSE)
# test <- test %>%
#   rename("protein" = "V1")
# test$protein <- gsub("=","_",test$protein)
# test2 <- read.delim2("../../mmseqs/master_table_prots.tsv", header = TRUE)
# test3 <- setdiff(test, test2) 
# test4 <- setdiff(test2, test) 

#test3 <- as.data.frame(setdiff(vibrant_best$id, vib_annos$id))


