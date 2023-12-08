######################################### Read Input ################################################
mmseqs <- read.csv2(file = "../../mmseqs/PlumeVent_mmseqs_clusters.tsv", sep = "\t")
#595,465 proteins clustered
#note that output from mmseqs2 is in long format so repeat proteins in first column
#mmseqs does not output singletons

############################################ mmseqs ####################################################

#Filter for non-self match in cluster
mmseqs_dif <- mmseqs %>%
  filter(cluster.representative != cluster.member)
#this should now be number of clusters
length(unique(mmseqs_dif$cluster.representative)) #74,963

# mmseqs_dif <- mmseqs_dif %>%
#   add_count(cluster.representative)

#this doesn't do what I thought:
# #Filter for duplicates in cluster representative
# mmseqs_noSingle <- mmseqs %>%
#   add_count(cluster.representative) %>% 
#   filter(n!=1) %>%
#   select(-n)
#

#change format to longer
mmseqs_dif_long <- mmseqs_dif %>%
  group_by(cluster.representative) %>%
  mutate(id = cur_group_id()) %>%
  pivot_longer(cols = c('cluster.representative', 'cluster.member')) %>%
  select(-name) %>%
  rename("genome" = "value") %>%
  ungroup() %>%
  group_by(id) %>%
  unique() %>%
  ungroup()

##### Find which clusters have proteins from Plume vs Vent
#substitute names in second column
#gsub to replace Plume names
mmseqs_PV <- mmseqs_dif_long
mmseqs_PV$genome <- gsub(".*Lau_Basin.*","Plume",mmseqs_PV$genome) #the placement of the periods is crucial for replacing whole string
mmseqs_PV$genome <- gsub(".*Cayman.*","Plume",mmseqs_PV$genome)
mmseqs_PV$genome <- gsub(".*Guaymas_Basin.*","Plume",mmseqs_PV$genome)
mmseqs_PV$genome <- gsub(".*Axial.*","Plume",mmseqs_PV$genome)

#Changing names to general site to see clusters from different sites (not same vent field)
mmseqs_GD <- mmseqs_dif_long
mmseqs_GD$genome <- gsub(".*Lau_Basin.*","Lau_Basin",mmseqs_GD$genome) #the placement of the periods is crucial for replacing whole string
mmseqs_GD$genome <- gsub(".*Cayman.*","Cayman",mmseqs_GD$genome)
mmseqs_GD$genome <- gsub(".*ELSC.*","Lau_Basin",mmseqs_GD$genome)
mmseqs_GD$genome <- gsub(".*Guaymas.*","Guaymas",mmseqs_GD$genome)
mmseqs_GD$genome <- gsub(".*Brothers.*","Brothers",mmseqs_GD$genome)
mmseqs_GD$genome <- gsub(".*MAR.*","MAR",mmseqs_GD$genome)
mmseqs_GD$genome <- gsub(".*EPR.*","EPR",mmseqs_GD$genome)
mmseqs_GD$genome <- gsub(".*Axial.*","Axial",mmseqs_GD$genome)

#separate the plume and vent dataframes 
mmseqs_P <- mmseqs_PV %>% filter(genome == "Plume") #22,656 plume
#remove Plume so can make all names left Vent
mmseqs_PV <- mmseqs_PV %>% filter(genome != "Plume")
mmseqs_PV$genome <- "Vent"
#Now put them back together
mmseqs_PV <- rbind(mmseqs_P, mmseqs_PV)

#sort by unique
mmseqs_PV_uniq <- unique(mmseqs_PV)
mmseqs_GD <- unique(mmseqs_GD)
#do any of the unique clusters now occur twice? (in cluster.rep column)
test <- mmseqs_PV_uniq %>% group_by(id) %>% filter(n()>1)
length(unique(test$id))
test2 <- mmseqs_GD %>% group_by(id) %>% filter(n()>1)
length(unique(test2$id))
#152 clusters that have a protein shared between vent and plume, made up of 763 proteins
#23,354-152= 23,202 clusters that have proteins that occur in geographically distinct vents
#84,225-763= 83,462 proteins in geo distinct vents

#make a list of the shared proteins
prots_PV_list <- as.data.frame(unique(test$id)) %>%
  rename("PV_clusters" = "unique(test$id)")
#grab the protein names shared between vent and plume
prots_PV <- mmseqs_dif_long[mmseqs_dif_long$id %in% prots_PV_list$PV_clusters,] #subset table using list of names

prots_GD_list <- as.data.frame(unique(test2$id)) %>%
  rename("GD_clusters" = "unique(test2$id)")
#grab the protein names shared between geo distinct vents
prots_GD <- mmseqs_dif_long[mmseqs_dif_long$id %in% prots_GD_list$GD_clusters,]
#remove the PV clusters from this
ids <- as.data.frame(unique(prots_PV$id)) %>%
  rename("id" = "unique(prots_PV$id)")
prots_GD <- prots_GD[!(prots_GD$id %in% ids$id),]

#remove the clusters in PV from GD
prots_GD_PV <- rbind(prots_GD, prots_PV)
prots_GD_PV <- prots_GD_PV %>%
  group_by(id,genome) %>%
  unique() %>%
  ungroup()

write_delim(prots_GD_PV, file = "Output/virus_proteins_mmseqsClusts_GD_PV.tsv", 
            col_names = TRUE, delim = "\t")
write_delim(prots_GD, file = "Output/virus_proteins_mmseqsClusts_GD.tsv", 
            col_names = TRUE, delim = "\t")
write_delim(prots_PV, file = "Output/virus_proteins_mmseqsClusts_PV.tsv", 
            col_names = TRUE, delim = "\t")

iris %>% 
  group_by(Species) %>% 
  arrange(Petal.Length) %>% 
  slice(1)

################################# get protein annotations  ####################################
#read metadata inputs to map annotation info
vent_vMAG <- read.delim2(file = "../VentVirus_Analysis/output/master_table_vMAGs.tsv")
plume_vMAG <- read.delim2(file = "~/Google Drive/My Drive/Faith/PlumeViruses/PlumeVirus_Analysis/output/master_table_vMAGs.tsv")
vent_vUnbinned <- read.delim2(file = "../VentVirus_Analysis/output/master_table_unbinned.tsv")
plume_vUnbinned <- read.delim2(file = "~/Google Drive/My Drive/Faith/PlumeViruses/PlumeVirus_Analysis/output/master_table_vUnbinned.tsv")
#cat together
vMAG_VP_master <- rbind(vent_vMAG, plume_vMAG)
vUnbinned_VP_master <- rbind(vent_vUnbinned, plume_vUnbinned)

#select master tables by column name - I want them in this order just because
vMAG_VP_master_sub <- select(vMAG_VP_master, c('protein', 'KO', 'AMG', 'KO.name',
                                               'Pfam', 'Pfam.name', 'VOG', 'VOG.name',
                                               'contig_length', 'checkv_quality', 'provirus',
                                               'completeness', 'contamination', 'warnings'))
vUnbinned_VP_master_sub <- select(vUnbinned_VP_master, c('protein', 'KO', 'AMG', 'KO.name',
                                                    'Pfam', 'Pfam.name', 'VOG', 'VOG.name',
                                                    'contig_length', 'checkv_quality', 'provirus',
                                                    'completeness', 'contamination', 'warnings'))
#combine
PV_master <- rbind(vMAG_VP_master_sub, vUnbinned_VP_master_sub)

#now get the proteins file to map the metadata onto
proteins_ventplume <- 
  mmseqs[mmseqs$cluster.representative %in% mmseqs_noSingle_uniq$cluster.representative,] #subset table using list of names

proteins_geodistinct <- 
  mmseqs[mmseqs$cluster.representative %in% mmseqs_noSingle_uniq$cluster.representative,] #subset table using list of names

#split vRhyme apart to change names
proteins_ventplume_vmags <- proteins_ventplume %>% filter(grepl('vRhyme', cluster.member))
proteins_geodistinct_vmags <- proteins_geodistinct %>% filter(grepl('vRhyme', cluster.member))

#get rid of vRhyme in mmseqs names for right join
proteins_ventplume_vmags <- proteins_ventplume_vmags %>% separate(cluster.member, c(NA, "cluster.member"), 
                                                      sep="(__)")
proteins_geodistinct_vmags <- proteins_geodistinct_vmags %>% separate(cluster.member, c(NA, "cluster.member"), 
                                                                  sep="(__)")
#remove Plume so can make all names left Vent
proteins_ventplume <- proteins_ventplume %>% filter(!grepl('vRhyme', cluster.member))
proteins_geodistinct <- proteins_geodistinct %>% filter(!grepl('vRhyme', cluster.member))
#Now put them back together
proteins_ventplume <- rbind(proteins_ventplume, proteins_ventplume_vmags)
proteins_geodistinct <- rbind(proteins_geodistinct, proteins_geodistinct_vmags)

#Finally, vlookup/right_join to merge
PV_master <- proteins_ventplume %>%
  dplyr::select("cluster.representative", "cluster.member") %>%
  left_join(PV_master, by = c("cluster.member" = "protein"))

PV_master_geodistinct <- proteins_geodistinct %>%
  dplyr::select("cluster.representative", "cluster.member") %>%
  left_join(PV_master, by = c("cluster.member" = "protein"))

################################# protein annotations with DRAMv  ####################################

#input
dramv <- read.delim2(file = "../../DRAMv/annotations_renamed.tsv", header = TRUE, sep = "\t")

#filter input for prots GD and PV



################################# write tables  ####################################
write.table(PV_master, file = "Output/mmseqs_proteinsVentPlume.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t")
write.table(PV_master_geodistinct, file = "Output/mmseqs_proteinsGeoDistinct.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t")
