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

tst2 <- prots_GD_PV %>%
  group_by(id) %>%
  mutate(yes_in_group = any(str_detect(genome, 'Axial'))) %>%
  filter(yes_in_group == TRUE)

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

#use prots_PV to get list of Plume Vent proteins and subset from PV_master
PV_prot_list <- unique(prots_PV$genome)
PV_annos <- PV_master %>% filter(PV_master$protein %in% PV_prot_list)

################################# protein annotations with DRAMv  ####################################

#input
dramv <- read.delim2(file = "../../DRAMv/annotations_renamed.tsv", header = TRUE, sep = "\t")

#filter input for prots GD and PV
dramv_pdgd <- dramv %>% filter(dramv$protein %in% prots_GD_PV$genome)

#zoom in on KEGG annotations
dramv_pdgd_kegg <- dramv_pdgd
dramv_pdgd_kegg[dramv_pdgd_kegg == ''] <- NA
dramv_pdgd_kegg <- dramv_pdgd_kegg %>% 
  drop_na(ko_id)

#dramv annos for PV proteins
PV_prot_list_ne <- gsub("=","_",PV_prot_list)
PV_annos_dram <- dramv %>% filter(dramv$protein %in% PV_prot_list_ne)

################################# protein annotations with PHROGs  ####################################

#PHROGs annos for PV proteins
#phrogs results
phrogs <- read.delim2(file = "../../PHROGs/convert_alis_80cov_75pi.tsv", header = TRUE, sep = "\t")
#mapping file
phrogs_map <- read.delim2(file = "../../PHROGs/phrog_annot_v4.tsv", header = TRUE, sep = "\t")

#remove phrogs string
phrogs$query <- gsub("phrog_", "", phrogs$query)
phrogs$query <- as.integer(phrogs$query)
phrogs$evalue <- as.numeric(phrogs$evalue)
#choose best phrogs annotation based on e value so 1 anno per scaffold
phrogs <- phrogs %>%
  group_by(target) %>%
  filter(evalue == min(evalue))
#vlookup
phrogs_map <- phrogs_map %>%
  dplyr::select("phrog", "color", "annot", "category") %>%
  right_join(phrogs, by = c("phrog" = "query"))

#drop excess cols
phrogs_map <- phrogs_map %>%
  select(c("phrog","color","annot","category","target", "evalue"))

#get only PV proteins
phrogs_map_PV <- phrogs_map %>%
  filter(phrogs_map$target %in% PV_prot_list)

################################### protein annotations with VIBRANT  #####################################################

#read in the input
#vibrant <- read.delim2(file = "../../VIBRANT_annos/52samples_VIBRANT_annos.tsv", header = TRUE, na.strings=c("","NA"))
vibrant <- read_tsv(file = "../../VIBRANT_annos/52_VIBRANT_annos.txt")
#change several columns to numeric
# vibrant %>% mutate_at(c('KO.evalue', 'KO.score', 'Pfam.evalue', 'Pfam.score',
#                         'VOG.evalue', 'VOG.score'), as.numeric) %>% str()

#remove proteins that have no annotation
# vibrant <- vibrant %>%
#   group_by(protein) %>%
#   filter(!(all(is.na(KO) & is.na(Pfam) & is.na(VOG))))
#192,991 of the original 595,578 proteins have an annotation

#add col of best bit score between KEGG, Pfam, and VOG
# vibrant_filt <- vibrant %>% 
#   mutate(highest_score = pmax(`Pfam score`, `VOG score`, `KO score`, na.rm = T))

library(data.table)

vibrant <- data.table(vibrant)

vibrant_filter <- vibrant[order(-score, evalue), .SD[1], by = protein]
#filter
# vibrant_filt <- vibrant %>%
#   group_by(protein) %>%
#   arrange(desc(score), evalue) %>%
#   slice(1) %>%
#   ungroup()

tst <- vibrant_filter[grepl("Axial_Plume", protein), ]

#read in Spencer-created file
vibrant <- read.csv(file = "../../VIBRANT_annos/Maggie_Proteins.csv")

#fixing the one that was dropped - tie breaker was higher viral score
vibrant <- vibrant %>%
  mutate(
    best = ifelse(dat.protein == "Lau_Basin_Abe_k95_1492896_flag=1_multi=24.9933_len=5352_6", "Pfam score", best),
    best_name = ifelse(dat.protein == "Lau_Basin_Abe_k95_1492896_flag=1_multi=24.9933_len=5352_6", "Phage Tail Collar Domain", best_name)
  )

#check how many of each
table(vibrant$best)
#KOs  31,382 | 30,192 Pfams | 131,417 VOGs
#192,991 proteins annotated with these 3 databases compared to starting 595,578 proteins

################################# write tables  ####################################
write.table(PV_master, file = "Output/mmseqs_proteinsVentPlume.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t")
write.table(PV_master_geodistinct, file = "Output/mmseqs_proteinsGeoDistinct.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t")

########################## adjust mmseqs file so have version of protein name without vRhyme ####################################
#the KEGG annotations have no vRhyme names because annotations were done before vRhyme so creating column
#without vRhyme name to map to that

#create new mmseqs file
mmseqs_dif_long_names <- mmseqs_dif_long
#new col
mmseqs_dif_long_names$genome_vRhyme <- mmseqs_dif_long_names$genome
#separate col
mmseqs_dif_long_names <- mmseqs_dif_long_names %>% separate(genome_vRhyme, c(NA, "genome_vRhyme"), sep = "__")
#put non vRhyme names back in the column
mmseqs_dif_long_names <- mmseqs_dif_long_names %>%
  mutate(genome_vRhyme = if_else(is.na(genome_vRhyme), genome, genome_vRhyme))

#get mmseqs names version of just GD PV proteins
mmseqs_GD_PV <- mmseqs_dif_long_names[mmseqs_dif_long_names$id %in% prots_GD_PV$id,]
  

# #write the output to get list of protein names without vRhyme...and just to have the file
# write.table(mmseqs_dif_long_names, file = "Output/mmseqs_long.tsv", col.names = TRUE,
#             row.names = FALSE, sep = "\t", quote = FALSE)

########### map VIBRANT output onto mmseqs clusters ##################

#map vibrant to mmseqs clusters with vlookup right_join
mmseqs_vibrant <- vibrant %>%
  dplyr::select(c(dat.protein, best_name)) %>%
  right_join(mmseqs_dif_long_names, by = c("dat.protein" = "genome_vRhyme")) %>%
  rename("protein_noVrhyme" = "dat.protein",
         "prot_name" = "best_name") %>%
  filter(!is.na(prot_name))
#81,994 mmseqs proteins with vibrant annotation

#create a Site column
mmseqs_vibrant <- mmseqs_vibrant %>%
  mutate(Site = protein_noVrhyme) %>%
  separate(Site, c("Site", NA), sep = "_scaffold|_k95|_NODE|-38[0-9]|_S[0-9][0-9][0-9]|_M1[0-9]|_35[0-9]-[0-9][0-9][0-9]|_1[0-9][0-9]-[0-9][0-9][0-9]")

#filter that file for just the geo distinct proteins
#grab just proteins that are shared between deposit and plume and between distinct vent fields
mmseqs_vibrant <- mmseqs_vibrant[mmseqs_vibrant$id %in% prots_GD_PV$id,] #subset table using list of names
length(unique(mmseqs_vibrant$genome)) #40,261 proteins GD and PV have annotation of 84,261 GD and PV proteins
tail(names(sort(table(mmseqs_vibrant$id))), 1) #to view largest protein cluster

#see proteins that only occur 1-2 times
# Group by the values in the specified column
mmseqs_vibrant_rare <- mmseqs_vibrant %>%
  group_by(prot_name) %>%
  mutate(count = n()) %>%
  filter(count >= 1 & count <= 2)

########### use the file of mmseqs KEGG and phrogs annotations you made on the server to ##################
#134,975 have annotations compared to 210,191 in list...but I guess phrogs annotations could be unknown

#input of mmseqs proteins that have kegg and/or phrogs annotation - kegg done with VIBRANT
mmseqs_kegg_phrog <- read.delim2(file = "../../mmseqs/mmseqs_KEGG_phrogs_annos.txt",
                                 header = FALSE)
#input kegg mapping file
kegg <- read.delim2(file = "~/Google Drive/My Drive/databases/kegg_mappingfile.txt")


#filter for better hit based on bitscore
mmseqs_kegg_phrog$V3 <- as.numeric(mmseqs_kegg_phrog$V3)
#make filtered file
mmseqs_kegg_phrog_filt <- mmseqs_kegg_phrog %>%
  group_by(V2) %>%
  filter(V3 == max(V3))
#confirm now no duplicate names
length(unique(mmseqs_kegg_phrog_filt$V2)) #134,975 unique proteins, 61,768 clusters with kegg/phrog anno

#map on the mmseqs clusters
mmseqs_kegg_phrog_filt <- mmseqs_dif_long_names %>%
  dplyr::select(id, genome_vRhyme) %>%
  right_join(mmseqs_kegg_phrog_filt, by = c("genome_vRhyme" = "V2"))
##add annotation description
#kegg
kegg <- kegg %>%
  select(c("KO","Description", "pathway")) %>%
  unique() %>%
  #filter_all(all_vars(. != "")) %>%
  rename("db_id" = "KO") %>%
  mutate(pathway = ifelse(db_id == "K11180", "Sulfur metabolism", pathway)) %>%
  mutate(pathway = ifelse(db_id == "K11181", "Sulfur metabolism", pathway)) %>%
  mutate(pathway = ifelse(db_id == "K23077", "Sulfur relay system", pathway))

#phrog
phrog <- phrogs_map %>%
  select(c("phrog", "annot", "category")) %>%
  unique() %>%
  rename("db_id" = "phrog") %>%
  rename("Description" = "annot") %>%
  rename("pathway" = "category")

#concatenate
kegg_phrog<-rbind(kegg,phrog)

#map more descriptive annos onto table
mmseqs_kegg_phrog_filt$V1 <- gsub("phrog_","",mmseqs_kegg_phrog_filt$V1)
#vlookup to join
mmseqs_kegg_phrog_filt <- kegg_phrog %>%
  dplyr::select(db_id, Description, pathway) %>%
  right_join(mmseqs_kegg_phrog_filt, by = c("db_id" = "V1"))

#create a Site column
mmseqs_kegg_phrog_filt <- mmseqs_kegg_phrog_filt %>%
  mutate(Site = genome_vRhyme) %>%
  separate(Site, c("Site", NA), sep = "_scaffold|_k95|_NODE|-38[0-9]|_S[0-9][0-9][0-9]|_M1[0-9]|_35[0-9]-[0-9][0-9][0-9]|_1[0-9][0-9]-[0-9][0-9][0-9]")

write.table(mmseqs_kegg_phrog_filt, file = "Output/mmseqs_kegg_phrog_filtered.tsv", col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)

#subset kegg and phrog file to determine how many clusters have examples where the annotations do not match
mmseqs_prot_mismatches <- mmseqs_kegg_phrog_filt %>%
  group_by(id) %>%
  filter(!all(pathway == first(pathway)))
length(unique(mmseqs_prot_mismatches$id))
#27,164 clusters of 61,768 clusters have annotations that are different - a lot of them are something annotated and 
#then "unknown" so maybe the unknowns are being assigned a function? 0:-)

#grab just proteins that are shared between deposit and plume and between distinct vent fields
mmseqs_kegg_phrog_filt <- mmseqs_kegg_phrog_filt[mmseqs_kegg_phrog_filt$id %in% prots_GD_PV$id,] #subset table using list of names
length(unique(prots_GD_PV$genome)) #84,261 proteins GD and PV
length(unique(mmseqs_kegg_phrog_filt$genome_vRhyme)) #55,777 proteins annotated
tail(names(sort(table(mmseqs_dif_long_names$id))), 1) #to view largest protein cluster

#get potential aux annots
mmseqs_kegg_phrog_aux <- mmseqs_kegg_phrog_filt %>%
  filter(pathway == "moron, auxiliary metabolic gene and host takeover")
  
#choose best rep of the cluster based on bitscore
mmseqs_kegg_phrog_filt2 <- mmseqs_kegg_phrog_filt %>%
  group_by(id) %>%
  filter(all(pathway == first(pathway)) | V3 == max(V3))
  ungroup()
  
###################################### plot the cluster annos #####################################################

#Plotting with all annos, even if multiple per cluster

plot <- mmseqs_vibrant 

#more general site name
plot$Site <- gsub("_Plume|_Seawater"," Seamount", plot$Site)
plot$Site <- gsub("_Diffuse|_LC|_NWCA|_NWCB|_UC"," Volcano", plot$Site)
plot$Site <- gsub("Cayman_Deep|Cayman_Shallow","Mid-Cayman Rise", plot$Site)
plot$Site <- gsub(".*EPR.*","East Pacific Rise", plot$Site)
plot$Site <- gsub(".*ELSC.*","Lau Basin Deposit", plot$Site)
plot$Site <- gsub("Guaymas_Basin","Guaymas Basin Plume", plot$Site)
plot$Site <- gsub(".*Guaymas_[0-9].*","Guaymas Basin Deposit", plot$Site)
plot$Site <- gsub(".*MAR.*","Mid-Atlantic Ridge", plot$Site)
plot$Site <- gsub(".*Lau_Basin.*","Lau Basin Plume", plot$Site)  
  
plot <- plot %>%
  select(prot_name, Site) %>%
  group_by(prot_name, Site) %>%
  count(prot_name) %>%
  ungroup()

#add percentage column so can plot where all sum to 100%
plot <- plot %>%
  group_by(Site) %>%
  mutate(percentage = n / sum(n) * 100)

#only grabbing top portion of data because so many low count categories
plot <- plot %>%
  group_by(Site) %>% #now group by site
  arrange(desc(n)) %>%
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 5) #top 13 rank is also acceptable - 15 unique pathways

############ Plotting with all annos, even if multiple per cluster ####################
  
#distinct colors
library(RColorBrewer)
n <- 24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

#wrap text for ggplot titles
# Helper function for string wrapping. 
# Default 20 character target width.
swr = function(string, nwrap=14) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Year
plot$Site = swr(plot$Site)
  
####The following produces Figure X, which was modified in Biorender
dev.off()
p <- plot %>%
  ggplot(aes(x = Site, y = as.numeric(percentage), fill = prot_name)) + #y = reorder(Site, value, sum)
  geom_bar(stat = "identity") +
  # scale_fill_viridis_d(begin = .5,
  #                      end = 0) +
  scale_fill_manual(values = col_vector) +
  labs(y = "Annotation", x = "") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 1, vjust = .5)
        axis.ticks.x = element_blank(),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        #legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_x_discrete(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  #geom_hline(yintercept = 16.5) +
  ggtitle("Selected Annotations of Proteins Shared Between Geographically Separated Vents") + #Change for top X grabbed
  facet_grid(.~Site, scales = "free") 
  #scale_y_discrete(limits=rev)
  #coord_flip()
p
  
#ggsave("Output/barplot_GD_PV_protein_annos.png", p, width = 13, dpi = 500,
#       bg = "transparent")
  
################### create donut plots showing % proteins shared between dif vents ######################

################### create input ######################

#need total proteins at each site
#then number of proteins from site shared with another --> convert into %

# library(stringi)
# #getting total number of proteins from each site
# #Add site column
# mmseqs <- mmseqs %>%
#   mutate(Site = cluster.member) %>%
#   separate(Site, c("Site", NA), sep = "_scaffold|_k95|_NODE|-38[0-9]|_S[0-9][0-9][0-9]|_M1[0-9]|_35[0-9]-[0-9][0-9][0-9]|_1[0-9][0-9]-[0-9][0-9][0-9]|_M[0-9]|_A[0-9]|_T[0-9]")
# mmseqs$Site <- stri_replace_all_regex(mmseqs$Site, pattern=c("vRhyme_[0-9]__", "vRhyme_[0-9][0-9]__",
#                                                        "vRhyme_[0-9][0-9][0-9]__"),
#                                                 replacement='',
#                                                 vectorize=FALSE)


#create Site file in mmseqs GD nad PV file
mmseqs_GD_PV <- mmseqs_GD_PV %>%
  mutate(Site = genome_vRhyme) %>%
  separate(Site, c("Site", NA), sep = "_scaffold|_k95|_NODE|-38[0-9]|_S[0-9][0-9][0-9]|_M1[0-9]|_35[0-9]-[0-9][0-9][0-9]|_1[0-9][0-9]-[0-9][0-9][0-9]|_M[0-9]|_A[0-9]|_T[0-9]")

#sort to avoid headaches
mmseqs_GD_PV <- mmseqs_GD_PV %>%
  group_by(id) %>%
  arrange(Site, .by_group = T) %>%
  ungroup()

#change names to simpler ones
mmseqs_GD_PV <- mmseqs_GD_PV %>%
  mutate(GenSite = Site) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*Lau_Basin.*", "Lau Basin Plume")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*Brothers.*", "Brothers Volcano")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*ELSC.*", "Lau Basin Deposit")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*Axial.*", "Axial Seamount")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*Cayman.*", "Mid-Cayman Rise")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*EPR.*", "East Pacific Rise")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., "Guaymas_Basin", "Guaymas Basin Plume")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*Guaymas_[0-9].*", "Guaymas Basin Deposit")) %>%
  mutate_at(vars(GenSite), ~ str_replace(., ".*MAR.*", "Mid-Atlantic Ridge"))

#get total protein count
prot_totals <- mmseqs_GD_PV %>%
  group_by(GenSite) %>%
  count()

library(widyr)
#pairwise_count to get counts of pairs within a group! 
mmseqs_GD_PV_counts <- pairwise_count(mmseqs_GD_PV, GenSite, id, sort = TRUE) 

#group by and count - for GD PV proteins, not all proteins that exist in the dataset
# #remove reverse duplicates if needed
# tst <- tst %>%
#   mutate(nv1 = paste0(item1, item2),
#          nv2 = paste0(item2, item1)) %>% 
#   unique_pairs("nv1", "nv2") %>% 
#   select(-nv1, -nv2)
  
#map total proteins for that sample
mmseqs_GD_PV_counts <- prot_totals %>%
  dplyr::select(GenSite, n) %>%
  right_join(mmseqs_GD_PV_counts, by = c("GenSite" = "item1"))
mmseqs_GD_PV_counts <- mmseqs_GD_PV_counts %>%
  rename("total" = "n.x",
         "count" = "n.y",
         "item1" = "GenSite") %>%
  mutate(prop = count/total*100)
#proportion isn't making sense with 7 site cateogry...because protein hits within group aren't showing up?
#i.e. Axial Seawater to Axial Plume would just come out as Axial Seamount

#add general site name
# mmseqs_GD_PV_counts <- mmseqs_GD_PV_counts %>%
#   mutate(GenSite = item1) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*Lau_Basin.*", "Lau_Basin Plume")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*Brothers.*", "Brothers Volcano")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*ELSC.*", "Lau Basin Deposit")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*Axial.*", "Axial Seamount")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*Cayman.*", "Mid-Cayman Rise")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*EPR.*", "East Pacific Rise")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., "Guaymas_Basin", "Guaymas Basin Plume")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*Guaymas_[0-9].*", "Guaymas Basin Deposit")) %>%
#   mutate_at(vars(GenSite), ~ str_replace(., ".*MAR.*", "Mid-Atlantic Ridge"))
  
################### plot ######################

# #calculate quantiles of interest
# library(data.table)
# mm_counts <- mmseqs_GD_PV_counts; setDT(mm_counts)
# mm_counts <- mm_counts[, .N, by = .(item1, item2, prop, GenSite)]
# 
# diam <- diamonds; setDT(diam)
# tabulated <- diam[, .N, by = .(cut, color, clarity)]

colors <- c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478", "East Pacific Rise" = "#CE9A28",
            "Guaymas Basin Deposit" = "#28827A", "Guaymas Basin Plume" = "#499e97", 
            "Lau Basin Deposit" = "#72a0db","Lau Basin Plume" = "#3F78C1","Mid-Atlantic Ridge" = "#8c510a",
            "Mid-Cayman Rise" = "#000000")

# Make the plot
dev.off()
p <- ggplot(mmseqs_GD_PV_counts, aes(x=2, y=count, fill=item2)) +
  geom_bar(position = 'fill', stat = 'identity')  +
  facet_wrap(~item1) + 
  theme_void() +
  xlim(0.5, 2.5) +
  coord_polar(theta = 'y') + 
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = colors, name = "Site")
p


data <- data.frame(
  category=c("A", "B", "C"),
  count=c(10, 60, 30)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)
# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))
# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

tst<-mmseqs_GD_PV_counts
tst <- tst %>%
  group_by(item1) %>%
  mutate(fraction = count/sum(count)) %>%
  ungroup()
# Compute the cumulative percentages (top of each rectangle)
tst <- tst %>%
  group_by(item1) %>%
  mutate(ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1)),
         labelPosition = (ymax + ymin) / 2,
         label = paste0(count)) %>%
  ungroup()
  
# Make the plot
dev.off()
p <- ggplot(tst, aes(x=2, y=count, fill=item2, ymax=ymax, ymin=ymin)) +
  geom_bar(position = 'fill', stat = 'identity')  +
  facet_wrap(~item1) + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=2.5) + 
  theme_void() +
  theme(panel.spacing = unit(1.5, "lines")) +
  xlim(0.5, 2.5) +
  coord_polar(theta = 'y') + 
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = colors, name = "Site")
p
  
ggsave(p, filename = "Output/protein_clust_donuts.svg", width = 11, height = 6)
ggsave(p, filename = "Output/protein_clust_donuts.png")

#and of course the text getting cut off in coord polar is not fixable with code: https://github.com/tidyverse/ggplot2/issues/2698
  
  
  

  
  