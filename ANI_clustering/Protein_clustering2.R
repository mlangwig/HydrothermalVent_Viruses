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
#595,416 viral proteins were clustered
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
length(unique(mmseqs$cluster.representative)) #74,960
#number of proteins in all clusters
length(unique(mmseqs$cluster.member)) #135,225

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
# 138 of these also have proteins from geographically distinct vents

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

################################ mmseqs PD GD add annotations ###################################

#add the mmseqs annotation info to the PD and GD cluster file to see how many have annotations/what they are
mmseqs_PD_GD_annos <- mmseqs_long_meta_annos %>%
  dplyr::select(c("anno_id", "anno", "protein")) %>%
  right_join(mmseqs_PD_GD, by = c("protein" = "genome_vRhyme")) %>%
  drop_na() %>%
  rename("protein_no_vRhyme" = "protein",
         "protein" = "genome") %>%
  select(c(id, protein, protein_no_vRhyme,PD, Site, anno_id, anno))

#column with specific site name for counting
mmseqs_PD_GD_annos$SiteDetail <- mmseqs_PD_GD_annos$protein
mmseqs_PD_GD_annos <- mmseqs_PD_GD_annos %>%
  mutate(SiteDetail = gsub(".*Lau_Basin.*", "Lau Basin Plume", SiteDetail),
         SiteDetail = gsub(".*ELSC.*", "Lau Basin Deposit", SiteDetail), #not distinguishing Laus here bc PD will get at that
         SiteDetail = gsub(".*Guaymas_Basin.*", "Guaymas Basin Plume", SiteDetail),
         SiteDetail = gsub(".*Guaymas_[0-9].*", "Guaymas Basin Deposit", SiteDetail),
         SiteDetail = gsub(".*Brothers.*", "Brothers Volcano", SiteDetail),
         SiteDetail = gsub(".*Cayman.*", "Mid-Cayman Rise", SiteDetail), 
         SiteDetail = gsub(".*EPR.*", "East Pacific Rise", SiteDetail),
         SiteDetail = gsub(".*Axial.*", "Axial Seamount", SiteDetail),
         SiteDetail = gsub(".*MAR.*", "Mid-Atlantic Ridge", SiteDetail)
  )

tst <- mmseqs_PD_GD_annos %>%
  group_by(id) %>%
  filter(all(SiteDetail %in% c("Brothers Volcano", "Lau Basin Deposit"))) %>% #keep only those rows that have ONLY Brothers and Lau Deposit
  ungroup()
length(unique(tst$id))

tst2 <- tst %>%
  group_by()

#write_delim(mmseqs_PD_GD_annos, file = "Output/mmseqs_proteins_PD_GD_annos.tsv",
#            col_names = TRUE, delim = "\t")

tst <- mmseqs_PD_GD_annos %>%
  count(anno)

#see just the plume deposit cluster annotations
mmseqs_PD_annos <- mmseqs_PD_GD_annos %>%
  filter(id %in% mmseqs_PD_ids$id)

mmseqs_PD_GD_annos_tmp <- mmseqs_PD_GD_annos %>%
  filter(PD == "Plume") %>%
  group_by(id) %>%
  filter(all(c("Guaymas", "Lau_Basin") %in% Site))

#################### mmseqs calculate and plot protein overlap among sites ################################

#column with specific site name for counting
mmseqs_PD_GD$SiteDetail <- mmseqs_PD_GD$genome_vRhyme
mmseqs_PD_GD <- mmseqs_PD_GD %>%
  mutate(SiteDetail = gsub(".*Lau_Basin.*", "Lau Basin Plume", SiteDetail),
         SiteDetail = gsub(".*ELSC.*", "Lau Basin Deposit", SiteDetail), #not distinguishing Laus here bc PD will get at that
         SiteDetail = gsub(".*Guaymas_Basin.*", "Guaymas Basin Plume", SiteDetail),
         SiteDetail = gsub(".*Guaymas_[0-9].*", "Guaymas Basin Deposit", SiteDetail),
         SiteDetail = gsub(".*Brothers.*", "Brothers Volcano", SiteDetail),
         SiteDetail = gsub(".*Cayman.*", "Mid-Cayman Rise", SiteDetail), 
         SiteDetail = gsub(".*EPR.*", "East Pacific Rise", SiteDetail),
         SiteDetail = gsub(".*Axial.*", "Axial Seamount", SiteDetail),
         SiteDetail = gsub(".*MAR.*", "Mid-Atlantic Ridge", SiteDetail)
  )

#get total protein count
prot_totals <- mmseqs_PD_GD %>%
  group_by(SiteDetail) %>%
  count()

library(widyr)
#pairwise_count to get counts of pairs within a group! 
mmseqs_PD_GD_counts <- pairwise_count(mmseqs_PD_GD, SiteDetail, id, sort = TRUE) 

#checking with easy examples to see if counts are correct
# test <- mmseqs_PD_GD %>%
#   group_by(id) %>%
#   filter(any(SiteDetail == "East Pacific Rise") & any(SiteDetail == "Axial Seamount"))

#CHECK TO CONFIRM THAT COUNTS ARE MADE CORRECTLY
check <- mmseqs_PD_GD %>%
  group_by(id) %>%
  filter(all(SiteDetail %in% c("Axial Seamount", "Mid-Cayman Rise"))) %>%
  ungroup()
length(unique(check$id)) #19

check <- mmseqs_PD_GD %>%
  group_by(id) %>%
  filter(all(SiteDetail %in% c("Brothers Volcano", "Lau Basin Deposit"))) %>%
  ungroup()
length(unique(check$id)) #15,874

#check for 3 sites occurring together
# Define the required sites
required_sites <- c("Brothers Volcano", "Lau Basin Deposit", "Mid-Atlantic Ridge")
# Filter out rows based on SiteDetail and group by id
check <- mmseqs_PD_GD %>%
  group_by(id) %>%
  filter(setequal(unique(SiteDetail), required_sites)) %>%
  ungroup()
length(unique(check$id)) #1,007

######################################## UpSet plot ################################################

#based on tutorial here: https://github.com/const-ae/ggupset
#install.packages("ComplexUpset")
library(ComplexUpset)

dev.off()
mmseqs_PD_GD_plot <- mmseqs_PD_GD
mmseqs_PD_GD_plot <- mmseqs_PD_GD_plot %>% #test2 <- 
  group_by(id, SiteDetail) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = SiteDetail, values_from = count, values_fill = 0) %>%
  as.data.frame() %>%
  #mutate_if(is.numeric, as.numeric) %>%
  mutate(id = as.integer(id)) %>%
  mutate_at(vars(-id), ~ifelse(. > 1, 1, .)) %>%
  as.data.frame() %>%
  select(`Axial Seamount`:`East Pacific Rise`)

p <- ComplexUpset::upset(mmseqs_PD_GD_plot,
      rev(c("Brothers Volcano", "Lau Basin Deposit", "Mid-Atlantic Ridge",
        "Guaymas Basin Deposit", "East Pacific Rise", "Lau Basin Plume",
        "Guaymas Basin Plume", "Axial Seamount", "Mid-Cayman Rise")),
      sort_sets=FALSE,
      name = 'Site',
      min_size = 15,
      width_ratio=0.12,
      set_sizes = 
        upset_set_size() + ylab('Number of proteins'),
      base_annotations = list(
        'Intersection size'=intersection_size(
          text=list(size = 2)
        )
      ),
      themes = upset_modify_themes(
        list(
          #'intersections_matrix'=theme(text=element_text(angle=90)),
          'overall_sizes'=theme(axis.text.x = element_text(size = 8, angle = 90))
        )
      ),
      matrix = (
        intersection_matrix(geom=geom_point(shape='circle filled', size=3)) +
          scale_color_manual(values=c("Axial Seamount" = "#4F508C", "Guaymas Basin Plume" = "#63c2ba",
                                                            "Lau Basin Plume" = "#72a0db", "Guaymas Basin Deposit" = "#28827A",
                                                            "Brothers Volcano" = "#B56478", "Mid-Cayman Rise" = "#000000",
                                                            "Lau Basin Deposit" = "#3F78C1", "Mid-Atlantic Ridge" = "#8c510a",
                                                            "East Pacific Rise" = "#CE9A28"),
                                                   guide=guide_legend(override.aes=list(shape='circle'))
                             )),
                    queries = list(upset_query(set="Axial Seamount", fill="#4F508C"),
                                   upset_query(set="Guaymas Basin Plume", fill="#63c2ba"),
                                   upset_query(set="Lau Basin Plume", fill="#72a0db"),
                                   upset_query(set="Guaymas Basin Deposit", fill="#28827A"),
                                   upset_query(set="Brothers Volcano", fill="#B56478"),
                                   upset_query(set="Mid-Cayman Rise", fill="#000000"),
                                   upset_query(set="Lau Basin Deposit", fill="#3F78C1"),
                                   upset_query(set="Mid-Atlantic Ridge", fill="#8c510a"),
                                   upset_query(set="East Pacific Rise", fill="#CE9A28"))
      )
p

#ggsave(p, filename = "Output/protein_clust_UpSet.svg", width = 10)
#ggsave(p, filename = "Output/protein_clust_UpSet.png", width = 10)

################################# Creating a modified bar plot above UpSet ################################################

#values taken straight from UpSet barplot (don't judge me) + col of smallest total number of proteins for that cluster
df <- data.frame(value=c(15874,1007,900,779,639,610,539,481,463,444,184,154,141,127,119,116,113,91,81,60,
                         52,40,36,34,32,30,26,22,19,15,15),
                 total_prots=c(19122,3111,3111,1561,1919,3111,1232,1919,1232,1919,699,397,1919,699,699,
                               1232,397,699,1919,699,1756,699,1919,1919,397,699,397,699,397,397,397))

df <- df %>%
  mutate(values_norm = (value/total_prots*100)) %>%
  mutate(across(where(is.numeric), round, 0)) %>%
  mutate(order=row_number())

dev.off()
p2<-ggplot(data=df, aes(x=factor(order),y=values_norm)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=values_norm), vjust=-0.25, size = 6)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_discrete(expand = c(0,0))+
  ylab("Percent of normalized shared protein clusters")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 20),
        axis.ticks.y = element_line(colour = "lightgrey"),
        panel.background = element_blank(),
        panel.grid.minor = element_line(color = "lightgrey"), 
        panel.grid.major = element_line(color = "lightgrey"))
p2

ggsave(p2, filename = "Output/protein_clust_UpSet_barplot_norm.svg", width = 14, height = 5)
ggsave(p2, filename = "Output/protein_clust_UpSet_barplot_norm.png", width = 14, height = 5)

################### unused

################################## donut plot ############################################


# colors <- c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478", "East Pacific Rise" = "#CE9A28",
#             "Guaymas Basin Deposit" = "#28827A", "Guaymas Basin Plume" = "#63c2ba", ##499e97
#             "Lau Basin Deposit" = "#3F78C1","Lau Basin Plume" = "#72a0db","Mid-Atlantic Ridge" = "#8c510a",
#             "Mid-Cayman Rise" = "#000000")
# 
# #new col for calculating fraction of total
# mmseqs_PD_GD_counts <- mmseqs_PD_GD_counts %>%
#   group_by(item1) %>%
#   mutate(fraction = n/sum(n)) %>%
#   ungroup()
# # Compute the cumulative percentages (top of each rectangle)
# mmseqs_PD_GD_counts <- mmseqs_PD_GD_counts %>%
#   group_by(item1) %>%
#   mutate(ymax = cumsum(fraction),
#          ymin = c(0, head(ymax, n=-1)),
#          labelPosition = (ymax + ymin) / 2,
#          label = paste0(n)) %>%
#   ungroup()
# 
# dev.off()
# p <- ggplot(mmseqs_PD_GD_counts, aes(x=2, y=n, fill=item2, ymax=ymax, ymin=ymin)) +
#   geom_bar(position = 'fill', stat = 'identity')  +
#   facet_wrap(~item1) + 
#   geom_text( x=2, aes(y=labelPosition, label=label), size=2.5, color = "white") + 
#   theme_void() +
#   theme(panel.spacing = unit(1.5, "lines")) +
#   xlim(0.5, 2.5) +
#   coord_polar(theta = 'y') + 
#   labs(x=NULL, y=NULL) +
#   scale_fill_manual(values = colors, name = "Site")
# p
# #the spacing of the number text in the donuts is off but good enough, I'm going to fix
# #in Inkscape
# 
# #ggsave(p, filename = "Output/protein_clust_donuts.svg", width = 11, height = 6)
# #ggsave(p, filename = "Output/protein_clust_donuts.png")



# #movies example
# movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )

# #install.packages("ggupset")
# library(ggupset)
# 
# #for ggupset
# test <- mmseqs_PD_GD %>%
#   select(c(id, SiteDetail)) %>%
#   unique() %>%
#   group_by(id) %>%
#   #count(id)
#   summarize(Sites = list(SiteDetail))
# # test <- mmseqs_PD_GD %>%
# #   group_by(id) %>%
# #   summarize(Sites = list(SiteDetail))
# 
# prot_totals <- prot_totals %>%
#   rename("prot_totals" = "n")
# #for UpSetR  
# # test <- prot_totals %>%
# #   right_join(mmseqs_PD_GD, by = c("SiteDetail" = "SiteDetail"))
# test1 <- mmseqs_PD_GD
# test1 <- test1 %>% filter(id == "74916")

# #creating metadata for the UpSetR plot
# sets <- as.data.frame(colnames(test2))
# sets <- sets %>%
#   rename("Sites" = "colnames(test2)") %>%
#   slice(-1) %>%
#   mutate(Sites_Map = Sites)
# PD <- data.frame(Type = c("Plume", "Plume", "Plume", "Deposit", "Deposit", "Plume",
#                "Deposit", "Deposit", "Deposit"))
# metadata <-cbind(sets, PD)

# set.metadata = list(data = metadata, 
#                     plots = list(list(type = "text", column = "Type", assign = 9,
#                                       colors = c(Plume = "navy", Deposit = "black")),
#                                  list(type = "matrix_rows", column = "Sites_Map", 
#                                       colors = c("Axial Seamount" = "#4F508C", "Guaymas Basin Plume" = "#63c2ba",
#                                                  "Lau Basin Plume" = "#72a0db", "Guaymas Basin Deposit" = "#28827A",
#                                                  "Brothers Volcano" = "#B56478", "Mid-Cayman Rise" = "#000000",
#                                                  "Lau Basin Deposit" = "#3F78C1", "Mid-Atlantic Ridge" = "#8c510a",
#                                                  "East Pacific Rise" = "#CE9A28"), alpha = 0.8)))


# #install.packages("UpSetR")
# library(UpSetR)
# test <- mmseqs_PD_GD
# test <- test %>% #test <- 
#   distinct(id, SiteDetail) %>% #distinct(id, SiteDetail)
#   unnest(cols = SiteDetail) %>%
#   mutate(SiteMember=1) %>%
#   pivot_wider(names_from = SiteDetail, values_from = SiteMember, values_fill = list(SiteMember = 0)) %>%
#   as.data.frame() #%>%
#   UpSetR::upset(sets = c("Axial Seamount", "Brothers Volcano", "East Pacific Rise", 
#                          "Guaymas Basin Deposit", "Guaymas Basin Plume", 
#                          "Lau Basin Deposit", "Lau Basin Plume", "Mid-Atlantic Ridge",
#                          "Mid-Cayman Rise"), keep.order = TRUE)
# dev.off()
# test2 <- mmseqs_PD_GD
# test2 %>% #test2 <- 
#   group_by(id, SiteDetail) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   pivot_wider(names_from = SiteDetail, values_from = count, values_fill = 0) %>%
#   as.data.frame() %>%
#   mutate_if(is.numeric, as.numeric) %>%
#   mutate(id = as.integer(id)) %>%
#   mutate_at(vars(-id), ~ifelse(. > 1, 1, .)) %>%
#   as.data.frame() %>%
#   UpSetR::upset(sets = c("Brothers Volcano", "Lau Basin Deposit", "Mid-Atlantic Ridge",
#                          "Guaymas Basin Deposit", "East Pacific Rise", 
#                          "Lau Basin Plume", "Guaymas Basin Plume", "Axial Seamount", 
#                          "Mid-Cayman Rise"),
#                 keep.order = TRUE,
#                 nintersects = NA,
#                 order.by = "freq",
#                 sets.x.label = "Number of proteins",
#                 queries = list(
#                   upset_query(set="Brothers Volcano", fill="#B56478")
#                 ))

#very simple upset plot
# dev.off()
# p <- ggplot(test, aes(x = Sites)) +
#   geom_bar() +
#   #geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
#   scale_x_upset(n_intersections = 20)
# p

# test <- read.delim2("../../VirusGenomes/49962_faas/virus_prots_list_renamed.txt", header = FALSE)
# test <- test %>%
#   rename("protein" = "V1")
# test$protein <- gsub("=","_",test$protein)
# test2 <- read.delim2("../../mmseqs/master_table_prots.tsv", header = TRUE)
# test3 <- setdiff(test, test2) 
# test4 <- setdiff(test2, test) 

#test3 <- as.data.frame(setdiff(vibrant_best$id, vib_annos$id))