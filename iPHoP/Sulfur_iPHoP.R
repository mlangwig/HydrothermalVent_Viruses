library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)
library(ggplot2)

setwd("~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/iPHoP/")

######################################### Read the input ##################################################

#iphop Vent
iphop_genome <- read.csv(file = "../../iPHoP/PlumeVent_Host_prediction_to_genome_m90.csv", header = TRUE)
iphop_genus <- read.csv(file = "../../iPHoP/PlumeVent_Host_prediction_to_genus_m90.csv", header = TRUE)

sulfur_mags <- read.csv(file = "../../sulfur_cyclers/subset_sulfurHMM_MAGs_final_uniq_namesFixed.txt", header = FALSE)

mag_sulfur_hits <- read.delim(file = "input/dissim_sulfur_MAGs_GeneContent.txt", header = TRUE)
mapping_acc <- read.delim(file = "input/hmm_mapping.txt", header = FALSE)

################### vlookup to match iphop genus and genome output ###############################

#remove ;s_ in gtdbtk classification for mapping
iphop_genome <- iphop_genome %>% separate(Host.taxonomy, c("Host.taxonomy", NA), sep= "(?=;s__)")

#map MAG name from genome file to genus file, matching by virus name and host taxonomy
iphop_genus_genome <- iphop_genome %>%
  dplyr::select(Virus, Host.genome, Host.taxonomy) %>%
  right_join(iphop_genus, iphop_genome, by = c("Virus" = "Virus", "Host.taxonomy" = "Host.genus"))

#drop matches to GCA genomes aka not MAGs
iphop_genus_genome <- iphop_genus_genome %>% 
  filter(!grepl("GCA_|GCF",Host.genome)) %>%
  mutate(across('Host.genome', str_replace, '.fna.noVirusContam', '')) %>%
  drop_na(Host.genome)  #drop the NAs in host.genome column bc appear tax doesnt match MAG

############################# subset for sulfur MAGs only ##################################################

iphop_genus_genome_sulfur <- 
  iphop_genus_genome[iphop_genus_genome$Host.genome %in% sulfur_mags$V1,]

################################### Quality control iphop results ########################################

iphop_genus_genome_sulfur <- sites_iphop %>%
  dplyr::select(contig_id, KB, checkv_quality, provirus,
                completeness, contamination, warnings) %>%
  right_join(iphop_genus_genome_sulfur, by = c("contig_id" = "Virus"))

#I am removing viruses whose checkv quality was Not determined because in my manual inspections,
#this gets rid of a lot of junk

iphop_genus_genome_sulfur<-iphop_genus_genome_sulfur[!grepl("Not-determined", iphop_genus_genome_sulfur$checkv_quality),]

#Removing viruses with the warning "no viral genes detected" because my manual inspections suggest these are not viral
#or are poor enough quality that I don't want to keep

iphop_genus_genome_sulfur<-iphop_genus_genome_sulfur[!grepl("no viral genes detected", iphop_genus_genome_sulfur$warnings),]

#Now removing viruses ≤5kb because I am not sure I trust host predictions to viral fragments
#And I'd like the potential for more genomic context from the virus

#filter for viruses with genome >5 KB
iphop_genus_genome_sulfur<-subset(iphop_genus_genome_sulfur, iphop_genus_genome_sulfur$KB>=5)

#Remove viruses with contamination >20% because these don't look great

iphop_genus_genome_sulfur<-subset(iphop_genus_genome_sulfur, iphop_genus_genome_sulfur$contamination<=20)
#768 unique viruses infecting 704 unique sulfur cycling MAGs when filtering by all these quality metrics
#114 unique sulfur cycling microbial taxa

#map viral taxonomy
iphop_genus_genome_sulfur <- virus_tax %>%
  dplyr::select(genome, lineage) %>%
  right_join(iphop_genus_genome_sulfur, by = c("genome" = "contig_id"))
iphop_genus_genome_sulfur <- iphop_genus_genome_sulfur %>% rename("Virus" = "genome")

#filter for >50% complete
iphop_genus_genome_sulfur_50comp<-subset(iphop_genus_genome_sulfur, iphop_genus_genome_sulfur$completeness>=50)
#make sure taxonomy is correct because I noticed errors further down (this seems to be because iphop sometimes doesn't get as specific with taxonomy compared to gtdb results on MAG)
mag_gtdb_VentPlume$user_genome <- sub(".fasta.noVirusContam", "", mag_gtdb_VentPlume$user_genome)

iphop_genus_genome_sulfur_50comp <- mag_gtdb_VentPlume %>%
  dplyr::select(user_genome, classification) %>%
  right_join(iphop_genus_genome_sulfur_50comp, by = c("user_genome" = "Host.genome"))
iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>% rename("Host.genome" = "user_genome")

################################### Formatting to write table ########################################

#separate taxonomy
iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>% separate(classification, c("d", "p", "c", "o", "f", "g"), 
                                                              sep= ";")

iphop_genus_genome_sulfur_50comp$VirusSite <- iphop_genus_genome_sulfur_50comp$Virus # copy column
iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>% separate(VirusSite, c("VirusSite", NA),
                                                              sep= "(?=_NODE|_k95|_scaffold|_vRhyme)") #separate by NODE and k95
iphop_genus_genome_sulfur_50comp$MAGSite <- iphop_genus_genome_sulfur_50comp$Host.genome # copy column
iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>% separate(MAGSite, c("MAGSite", NA),
                                                                                  sep= "(?=_maxbin|_metabat|_UWMA|_HVA)") #separate by NODE and k95

#replace specific site to get general sites - do this for 2 columns at once
sulf_cols<-c("VirusSite", "MAGSite")
iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*ELSC.*","Lau_Basin")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*Brothers.*","Brothers_Volcano")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*Guaymas.*","Guaymas_Basin")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*MAR.*","Mid_Atlantic_Ridge")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*EPR.*","East_Pacific_Rise")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*Cayman.*","Mid_Cayman_Rise")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*Lau.*","Lau_Basin")) %>%
  mutate_at(vars(sulf_cols), ~ str_replace(., ".*Axial.*","Axial_Seamount"))

iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>% separate(List.of.methods, c("Method", NA),
                                                                                  sep= ";")

iphop_genus_genome_sulfur_50comp <- iphop_genus_genome_sulfur_50comp %>% separate(lineage, 
                                                              c("rV", "kV", "pV", "cV", "oV", "fV", "gV", "sV"),
                                                              sep= ";")

write.csv(file = "output/iphop_VentPlume_sulfur_50comp.csv", iphop_genus_genome_sulfur_50comp,
          row.names = FALSE, quote = FALSE)

########### Making input data for a bubble plot of sulfur cycling microbes and their viruses ########################################
mag_sulfur_hits <- read.delim(file = "input/dissim_sulfur_MAGs_GeneContent.txt", header = TRUE)

mag_sulfur_hits <- mapping_acc %>%
  dplyr::select(V1, V2) %>%
  right_join(mag_sulfur_hits, by = c("V1" = "accession")) %>%
  separate(protein, c("MAG", NA),
           sep= "(?=_NODE|_scaffold)")

#fix MAG names so they'll match the iphop table
mag_sulfur_Axial <- mag_sulfur_hits %>% filter(str_detect(MAG,"Axial"))
mag_sulfur_CaymanLauGuay <- mag_sulfur_hits %>% filter(str_detect(MAG, "Cayman|Lau|Guaymas_Basin"))
#remove strings after last 2 underscores for Cayman and Lau
mag_sulfur_CaymanLauGuay$MAG <- sub("_[^_]+$", "", mag_sulfur_CaymanLauGuay$MAG)
mag_sulfur_CaymanLauGuay$MAG <- sub("_[^_]+$", "", mag_sulfur_CaymanLauGuay$MAG)

mag_sulfur_NotAxialCaymanLauGuay <- mag_sulfur_hits %>% filter(!str_detect(MAG,"Axial|Cayman|Lau|Guaymas_Basin"))
#remove anything after last underscore
mag_sulfur_NotAxialCaymanLauGuay$MAG <- sub("_[^_]+$", "", mag_sulfur_NotAxialCaymanLauGuay$MAG) 

#combine them all together again with corrected names
mag_sulfur_hits <- rbind(mag_sulfur_Axial,mag_sulfur_NotAxialCaymanLauGuay,mag_sulfur_CaymanLauGuay)
#get rid of e value and bit score so can filter for unique number of hits to a gene
mag_sulfur_hits <- select(mag_sulfur_hits,-c(evalue,score))
mag_sulfur_hits <- unique(mag_sulfur_hits)

#map sulfur cycling genes of MAGs onto table for reformatting and plotting

mag_sulfur_plotting <- mag_sulfur_hits %>%
  dplyr::select(V1, V2, MAG) %>%
  right_join(iphop_genus_genome_sulfur_50comp, by = c("MAG" = "Host.genome")) %>%
  rename("accession" = "V1") %>%
  rename("gene" = "V2")

mag_sulfur_plotting <- mag_sulfur_plotting %>%
  select(gene, MAG, Virus, d, c, Method, VirusSite, MAGSite) %>%
  group_by(d, c, gene) %>%
  count(gene) %>%
  #arrange() %>%
  mutate(c = str_replace(c, "c__", "")) %>%
  mutate(d = str_replace(d, "d__", "")) %>%
  arrange() %>%
  na_if('')

mag_sulfur_plotting$c <- mag_sulfur_plotting$c %>% replace_na("Bacteria") #replace blank cell with Bac

#count number of occurrences in iphop_genus file of unique virus infecting MAG
sulf_mag_virus_count <- iphop_genus_genome_sulfur_50comp %>%
  group_by(c, Virus) %>%
  summarise(count = n_distinct(Virus)) %>%
  count(c, name = "num_infecting_viruses") %>%
  mutate(c = str_replace(c, "c__", "")) %>%
  na_if('') %>%
  mutate(num_infecting_viruses = paste0("(", num_infecting_viruses, ")"))
sulf_mag_virus_count$c <- sulf_mag_virus_count$c %>% replace_na("Bacteria") #replace blank cell with Bac

mag_sulfur_plotting <- sulf_mag_virus_count %>%
  dplyr::select(c, num_infecting_viruses) %>%
  right_join(mag_sulfur_plotting, by = c("c" = "c"))

mag_sulfur_plotting <- mag_sulfur_plotting %>% 
  unite(class, c(c, num_infecting_viruses), sep = " ", remove = FALSE)

################################### plotting ########################################

mag_sulfur_plotting$gene <- factor(mag_sulfur_plotting$gene, 
                                   levels = c("aprA", "dsrA", "dsrB", "dsrD",
                                           "phsA", "soxB", "soxC", "soxY",
                                           "fccB", "sor"))

# mag_sulfur_plotting$d <- factor(mag_sulfur_plotting$d)
mag_sulfur_plotting$c <- factor(mag_sulfur_plotting$c)

#log transform values for better visualization
mag_sulfur_plotting <- mag_sulfur_plotting %>% 
  mutate(log_n = log(n)) %>%
  mutate(log_n_plus1 = log_n+1)

# mag_sulfur_plotting_a <- mag_sulfur_plotting %>% 
#   filter(str_detect(d,"Archaea"))
# mag_sulfur_plotting_b <- mag_sulfur_plotting %>% 
#   filter(str_detect(d,"Bacteria"))

dev.off()
p <- ggplot(mag_sulfur_plotting, aes(y=class, x=gene))+
  geom_point(aes(size=log_n_plus1, color = log_n_plus1))+ 
  #scale_size_continuous(range = c(1, 7))+
  scale_color_continuous(guide="legend", 
                         type = "viridis")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=-90),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-.10),
        text = element_text(color="black"),
        legend.position="right")+
  scale_x_discrete(position="top")+
  xlab("")+
  ylab("Microbial Class (Number of Viruses)")+
  labs(size = "Number of genes",
       color = "Number of genes") +
  facet_grid(d ~ ., scales = "free_y", space = "free") + #wow this took me awhile - remember you can't factor the y axis and then get it to facet properly
  scale_y_discrete(limits=rev) # has to be this instead of factor
p  

ggsave("output/iphop_sulfur_50comp_genes.png", width = 7, height = 6)
  
  
# ########################## old ##########################
# #map
# iphop_VentPlume_sulfur <- sulfur_cyclers %>%
#   dplyr::select(Virus, MAG_w_sulfur, Method, GTDBtk_v1_5_0) %>%
#   right_join(iphop_VentPlume, by = c("Virus" = "Virus")) %>%
#   drop_na(MAG_w_sulfur)
# 
# iphop_VentPlume_sulfur_high_qual <- iphop_VentPlume_sulfur %>% filter(!str_detect(checkv_quality,
#                                                                                   "Low-quality"))
# 
# iphop_VentPlume_sulfur <- iphop_VentPlume_sulfur %>% rename("VirusSite" = "Site")
# iphop_VentPlume_sulfur$MAGSite <- iphop_VentPlume_sulfur$MAG_w_sulfur # copy column
# iphop_VentPlume_sulfur <- iphop_VentPlume_sulfur %>% separate(MAGSite, c("Site", NA),
#                                                               sep= "(?=_maxbin|_metabat)") #separate by NODE and k95
# iphop_VentPlume_sulfur <- iphop_VentPlume_sulfur %>% rename("MAGSite" = "Site")
# 
# #replace specific site to get general sites - do this for 2 columns at once
# sulf_cols<-c("VirusSite", "MAGSite")
# iphop_VentPlume_sulfur <- iphop_VentPlume_sulfur %>%
#   mutate_at(vars(sulf_cols), ~ str_replace(., ".*ELSC.*","Lau_Basin")) %>%
#   mutate_at(vars(sulf_cols), ~ str_replace(., ".*Brothers.*","Brothers_Volcano")) %>%
#   mutate_at(vars(sulf_cols), ~ str_replace(., ".*Guaymas.*","Guaymas_Basin")) %>%
#   mutate_at(vars(sulf_cols), ~ str_replace(., ".*MAR.*","Mid_Atlantic_Ridge")) %>%
#   mutate_at(vars(sulf_cols), ~ str_replace(., ".*EPR.*","East_Pacific_Rise"))
# 
# #separate taxonomy
# iphop_VentPlume_sulfur <- iphop_VentPlume_sulfur %>% separate(Host.genus, c("d", "p", "c", "o", "f", "g"), 
#                                                               sep= ";")
# # iphop_VentPlume_sulfur <- iphop_VentPlume_sulfur %>% separate(lineage, c("r", "k", "p", "c", "o", "f", "g", "s"), 
# #                                                               sep= ";")
