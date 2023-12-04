######################################### Read Input ################################################
mmseqs <- read.csv2(file = "Input/PlumeVent_mmseqs_clusters.tsv", sep = "\t")
#595,541 proteins clustered
#note that output from mmseqs2 is in long format so repeat proteins in first column
#mmseqs does not output singletons

############################################ mmseqs ####################################################

#Filter for non-self match in cluster
mmseqs_dif <- mmseqs %>%
  filter(cluster.representative != cluster.member)
#this should now be number of clusters
length(unique(mmseqs_dif$cluster.representative))

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
mmseqs_dif$cluster.member <- gsub(".*Lau_Basin.*","Lau_Basin",mmseqs_dif$cluster.member) #the placement of the periods is crucial for replacing whole string
mmseqs_dif$cluster.member <- gsub(".*Cayman.*","Cayman",mmseqs_dif$cluster.member)
mmseqs_dif$cluster.member <- gsub(".*ELSC.*","ELSC",mmseqs_dif$cluster.member)
mmseqs_dif$cluster.member <- gsub(".*Guaymas.*","Guaymas",mmseqs_dif$cluster.member)
mmseqs_dif$cluster.member <- gsub(".*Brothers.*","Brothers",mmseqs_dif$cluster.member)
mmseqs_dif$cluster.member <- gsub(".*MAR.*","MAR",mmseqs_dif$cluster.member)
mmseqs_dif$cluster.member <- gsub(".*EPR.*","EPR",mmseqs_dif$cluster.member)
mmseqs_dif$cluster.member <- gsub(".*Axial.*","Axial",mmseqs_dif$cluster.member)

#separate the plume and vent dataframes 
mmseqs_P <- mmseqs_PV %>% filter(genome == "Plume") #22,656 plume
#remove Plume so can make all names left Vent
mmseqs_PV <- mmseqs_PV %>% filter(genome != "Plume")
mmseqs_PV$genome <- "Vent"
#Now put them back together
mmseqs_PV <- rbind(mmseqs_P, mmseqs_PV)

#sort by unique
mmseqs_PV_uniq <- unique(mmseqs_PV)
#do any of the unique clusters now occur twice? (in cluster.rep column)
test <- mmseqs_PV_uniq %>% group_by(id) %>% filter(n()>1)
#163 clusters that have a protein shared between vent and plume, X 326 proteins that occur in both vent and plume
#49,362 proteins that occur in geographically distinct vents

#make a list of the shared proteins
prots_PV_list <- as.data.frame(unique(test$id)) %>%
  rename("PV_clusters" = "unique(test$id)")
#grab the protein names shared between vent and plume
prots_PV <- mmseqs_dif_long[mmseqs_dif_long$id %in% prots_PV_list$PV_clusters,] #subset table using list of names


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

################################# write tables  ####################################
write.table(PV_master, file = "Output/mmseqs_proteinsVentPlume.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t")
write.table(PV_master_geodistinct, file = "Output/mmseqs_proteinsGeoDistinct.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t")
