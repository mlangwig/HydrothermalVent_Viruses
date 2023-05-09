######################################### Read Input ################################################
mmseqs <- read.csv2(file = "Input/PlumeVent_mmseqs_clusters.tsv", sep = "\t")
#595,541 clusters, 460,233 without singletons

############################################ mmseqs ####################################################

#Filter for duplicates in cluster representative
mmseqs_noSingle <- mmseqs %>%
  add_count(cluster.representative) %>% 
  filter(n!=1) %>%
  select(-n)

#substitute names in second column
#gsub to replace Plume names
mmseqs_noSingle$cluster.member <- gsub(".*Lau_Basin.*","Plume",mmseqs_noSingle$cluster.member) #the placement of the periods is crucial for replacing whole string
mmseqs_noSingle$cluster.member <- gsub(".*Cayman.*","Plume",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*Guaymas_Basin.*","Plume",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*Axial.*","Plume",mmseqs_noSingle$cluster.member)

#Changing names to general site to see clusters from different sites (not same vent field)
mmseqs_noSingle$cluster.member <- gsub(".*Lau_Basin.*","Lau_Basin",mmseqs_noSingle$cluster.member) #the placement of the periods is crucial for replacing whole string
mmseqs_noSingle$cluster.member <- gsub(".*Cayman.*","Cayman",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*ELSC.*","ELSC",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*Guaymas.*","Guaymas",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*Brothers.*","Brothers",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*MAR.*","MAR",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*EPR.*","EPR",mmseqs_noSingle$cluster.member)
mmseqs_noSingle$cluster.member <- gsub(".*Axial.*","Axial",mmseqs_noSingle$cluster.member)

#separate the plume and vent dataframes 
mmseqs_noSingle_Plume <- mmseqs_noSingle %>% filter(cluster.member == "Plume")
#remove Plume so can make all names left Vent
mmseqs_noSingle <- mmseqs_noSingle %>% filter(cluster.member != "Plume")
mmseqs_noSingle$cluster.member <- "Vent"
#Now put them back together
mmseqs_noSingle <- rbind(mmseqs_noSingle, mmseqs_noSingle_Plume)

#sort by unique
mmseqs_noSingle_uniq <- unique(mmseqs_noSingle)
#do any of the unique clusters now occur twice? (in cluster.rep column)
mmseqs_noSingle_uniq <- mmseqs_noSingle_uniq %>% group_by(cluster.representative) %>% filter(n()>1)
#326 proteins that occur in both vent and plume
#49,362 proteins that occur in geographically distinct vents

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
