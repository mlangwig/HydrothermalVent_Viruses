##################### Integrating all virus output into 1 table #####################

setwd(dir = "~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/VentVirus_Analysis/")

library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)

##################################### Read in major inputs ##########################################

#checkV
checkv <- read.delim(file = "input/quality_summary.tsv", header = TRUE)
#remove cols
checkv <- checkv %>%
  select(-c(contig_length, provirus, proviral_length, completeness_method, kmer_freq))
table(checkv$checkv_quality)

#lifestyle, VIBRANT predicted
vib_type <- read.delim(file = "input/52samples_VIBRANT_type.tsv", header = TRUE)
vib_type$scaffold <- gsub("=","_",vib_type$scaffold)

#genome size, Seqkit
gensize <- read.table(file = "input/PlumeVentVirus_Seqkit_final.txt", sep = " ", header = TRUE)
#count range of 1-5kb genome sizes
count <- sum(gensize$sum_len >= 1000 & gensize$sum_len <= 5000)

#VIBRANT annotations of KEGG, Pfam, and VOGs
vib_annos <- read.delim(file = "../../VIBRANT_annos/52samples_VIBRANT_annos.tsv", sep = "\t", fill = TRUE, 
                        header = TRUE, na.strings=c("","NA"))
#replace = with _
vib_annos$protein <- gsub("=","_",vib_annos$protein)
vib_annos$scaffold <- gsub("=","_",vib_annos$scaffold)

#vMAG mapping vMAG name to scaffolds
vMAG_mapping <- read.table(file = "input/vMAG_scaffolds_5708_final.txt", sep = "\t")
#rename unnamed cols
vMAG_mapping <- vMAG_mapping %>%
  rename("vMAG" = "V1",
         "scaffold" = "V2") #new = old

#geNomad taxonomy
genomad_tax <- read.table(file = "input/49962_final_VentViruses_Nlinked_taxonomy_parsed.tsv", header = TRUE)
#remove cols
genomad_tax <- genomad_tax %>%
  select(-c(n_genes_with_taxonomy, agreement, taxid)) %>%
  separate(lineage, c("r", "k", "p", "c", "o", "f", "g", "s"), 
           sep= ";")
table(genomad_tax$r)
table(genomad_tax$c)

#iPHoP host predictions
iphop <- read.csv(file = "input/Host_prediction_to_genus_m90_49962.csv",
                  header = TRUE)

iphop_filt <- iphop %>%
  separate(Host.genus, c("d", "p", "c", "o", "f", "g"), sep= ";") %>%
  group_by(Virus) %>%
  filter(Confidence.score == max(Confidence.score, na.rm = TRUE)) %>% #note doesn't get rid of all duplicates because some have same confidence score
  filter(c == "c__Gammaproteobacteria")

virus_hallmarks <- read.delim2(file = "input/virus_hallmarks_sum.tsv")
  
####################################### Determine lytic vs lysogenic #############################################
################# for vMAGs that may have had dif types binned together

#remove 18 E coli virus contaminants
remove.list <- paste(c("ELSC_Bowl_M2_NODE_26941_length_5513_cov_204.284070",
                       "Cayman_Deep_k95_329100_flag_3_multi_740.0000_len_5481",
                       "ELSC_Abe_A3_NODE_22123_length_5513_cov_37.688452",
                       "Guaymas_Basin_k95_702414_flag_3_multi_156.0000_len_5481",
                       "Lau_Basin_Tahi_Moana_k95_522185_flag_3_multi_299.0000_len_5481",
                       "ELSC_Vai_Lili_V2_NODE_3240_length_5513_cov_24.324359",
                       "ELSC_Mariner_M17_NODE_20378_length_5513_cov_309.771259",
                       "ELSC_Bowl_M1_NODE_4552_length_5513_cov_398.155589",
                       "Cayman_Shallow_k95_556392_flag_3_multi_763.1417_len_5481",
                       "ELSC_Tui_Malila_T10_NODE_9861_length_5513_cov_130.186966",
                       "ELSC_Abe_A1_NODE_14649_length_5513_cov_85.977163",
                       "ELSC_Tui_Malila_T11_NODE_11702_length_5513_cov_423.693093",
                       "ELSC_Mariner_M10_NODE_9821_length_5513_cov_36.731526",
                       "Lau_Basin_Mariner_k95_379953_flag_3_multi_286.0000_len_5481",
                       "Lau_Basin_Abe_k95_1566522_flag_3_multi_441.0000_len_5481",
                       "ELSC_Tui_Malila_T2_NODE_29080_length_5513_cov_220.766803",
                       "Lau_Basin_Tui_Malila_k95_411308_flag_3_multi_152.0000_len_5481",
                       "Lau_Basin_Kilo_Moana_k95_205532_flag_3_multi_561.0000_len_5481"), collapse = '|')
vib_type <- vib_type %>%
  filter(!str_detect(scaffold, remove.list))

#add VIBRANT lytic lysogenic prediction to VIBRANT annotations table
vib_type_final <- vib_annos %>%
  right_join(vib_type, by = c("scaffold" = "scaffold")) 

#add vRhyme MAG names to VIBRANT scaffold IDs
vib_type_final <- vib_type_final %>%
  right_join(vMAG_mapping, by = c("scaffold" = "scaffold"))

#subset just the scaffold, vMAG, and type cols
vib_type_final <- vib_type_final %>%
  select(c(scaffold, type, vMAG)) %>%
  unique()
#count number of lytic and lysogenic scaffolds in vMAGs
count <- vib_type_final %>%
  group_by(vMAG) %>%
  summarize(lytic_count = sum(type == "lytic"),
            lysogenic_count = sum(type == "lysogenic"))
#which have ≥2 - get a list and split those vMAGs
count <- count %>%
  filter(lysogenic_count >= 2)
# write.table(count, file = "../../vRhyme/vMAG_114_lysogenic_list.txt", quote = FALSE,
#             row.names = FALSE, sep = "\t")

################# 

#now determine whether vMAGs are lytic or lysogenic based on presence of 1 lysogenic scaffold
vib_type_final <- vib_type_final %>%
  select(c('type','vMAG')) %>%
  unique()
#if when grouping by vMAG there is lytic+lysogenic, change all types for that vMAG to lysogenic
vib_type_vMAG <- vib_type_final %>%
  group_by(vMAG) %>%
  mutate(type = ifelse("lytic" %in% type & "lysogenic" %in% type, "lysogenic", type)) %>%
  ungroup() %>%
  unique()
table(vib_type_vMAG$type)
#569 lysogenic vMAGs, 5,139 lytic

#get VIBRANT types of non-vMAGs
vib_type_vUnbinned <- vib_type %>%
  select(c('scaffold', 'type')) %>%
  unique()
`%notin%` <- Negate(`%in%`)
vib_type_vUnbinned<-vib_type_vUnbinned[vib_type_vUnbinned$scaffold %notin% vMAG_mapping$scaffold,]
#should be 44,254 or the number of unbinned viruses
vib_type_vUnbinned <- vib_type_vUnbinned %>%
  select(c('type','scaffold')) %>%
  rename("vMAG" = "scaffold")
#make final VIBRANT type
vib_type <- rbind(vib_type_vMAG, vib_type_vUnbinned)
table(vib_type$type)
#2,391 lysogenic viruses, 47,571 lytic

####################################### Create the master tables #############################################

#get vMAG names in vib_annos
vib_annos <- vib_annos %>%
  left_join(vMAG_mapping, by = c("scaffold" = "scaffold"))
#replace NAs for unbinned with their scaffold name
vib_annos$vMAG <- ifelse(is.na(vib_annos$vMAG), vib_annos$scaffold, vib_annos$vMAG)
#remove E coli contaminant viruses
vib_annos <- vib_annos %>%
  filter(!str_detect(scaffold, remove.list))
length(unique(vib_annos$vMAG))
#49,962 viruses, same as total

#create master table and add CheckV output
master_table <- vib_annos %>%
  select(vMAG, protein:VOG.v.score) %>%
  right_join(checkv, by = c("vMAG" = "contig_id"))

#add seqkit genome size output
master_table <- gensize %>%
  select(file, num_seqs, min_len, avg_len, max_len, sum_len) %>%
  right_join(master_table, by = c("file" = "vMAG")) %>%
  rename("vMAG" = "file")

#add genomad taxonomy
master_table <- master_table %>%
  select(vMAG, protein:warnings, num_seqs:sum_len)
master_table <- genomad_tax %>%
  right_join(master_table, by = c("genome" = "vMAG")) %>%
  rename("vMAG" = "genome")
master_table <- master_table %>%
  select(vMAG, protein:warnings, num_seqs:sum_len, r:f) #only to family because no genus or species names

#add lifestyle type
master_table <- vib_type %>%
  right_join(master_table, by = c("vMAG" = "vMAG"))

#add number of virus hallmarks
master_table <- virus_hallmarks %>%
  right_join(master_table, by = c("gene" = "vMAG")) %>%
  rename("vMAG" = "gene")

# write.table(master_table, file = "output/master_table_VentPlumeViruses.tsv", col.names = TRUE,
#             quote = FALSE, row.names = FALSE, sep = "\t")

#make master table without proteins for accurate counting
master_table_noProtein <- master_table %>%
  select(vMAG:type, gene_count:f) %>%
  unique()

# write.table(master_table_noProtein, file = "output/master_table_VentPlumeViruses_noProtein.tsv", col.names = TRUE,
#             quote = FALSE, row.names = FALSE, sep = "\t")

#master table no proteins of med-quality and better only
master_table_noProtein_hq <- master_table %>%
  select(vMAG, type, gene_count:f) %>%
  unique() %>%
  filter(!(checkv_quality %in% c("Low-quality", "Not-determined")))
table(master_table_noProtein_hq$type)
#328 lysogenic, 1,505 lytic for med quality and better viruses

#create a master table reduced info for easier viewing
master_table_simple <- master_table %>%
  select(vMAG, scaffold, protein, KO, AMG, KO.name, Pfam, 
         Pfam.name, VOG, VOG.name, type, checkv_quality, completeness, contamination, sum_len:f)
# write.table(master_table_simple, file = "output/master_table_VentPlumeViruses_simple.tsv", col.names = TRUE,
#             quote = FALSE, row.names = FALSE, sep = "\t")

#add iphop host prediction
master_table_iphop <- iphop %>%
  right_join(master_table_noProtein, by = c("Virus" = "vMAG")) %>%
  rename("vMAG" = "Virus")

#add number of virus hallmarks
master_table_iphop <- virus_hallmarks %>%
  right_join(master_table_iphop, by = c("gene" = "vMAG")) %>%
  rename("vMAG" = "gene")

#filter iphop results for just those viruses with hallmarks
# master_table_iphop_filt <- master_table_iphop %>%
#   filter(total_hallmarks > 0) 

master_table_iphop_filt <- master_table_iphop
master_table_iphop_filt <- master_table_iphop_filt[!is.na(master_table_iphop_filt$Host.genus),]

#only keep phylum and class of the tax string
master_table_iphop_filt <- master_table_iphop_filt %>% separate("Host.genus", c("d", "p", "c", "o", "f", "g"), 
                                                            sep= ";")

#add Site labels
master_table_iphop_filt$Site <- master_table_iphop_filt$vMAG
master_table_iphop_filt <- master_table_iphop_filt %>% separate(Site, c("Site", NA), 
                                      sep= "_NODE|_scaffold|_vRhyme|_k95")
#add deposit and plume labels
master_table_iphop_filt$Locat <- master_table_iphop_filt$Site
master_table_iphop_filt$Locat <- gsub(".*Axial.*","Plume", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*Cayman.*","Plume", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub("Guaymas_Basin","Plume", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*Lau.*","Plume", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*Brothers.*","Deposit", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*ELSC.*","Deposit", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*EPR.*","Deposit", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*Guaymas.*","Deposit", master_table_iphop_filt$Locat)
master_table_iphop_filt$Locat <- gsub(".*MAR.*","Deposit", master_table_iphop_filt$Locat)

master_table_iphop_filt3 <- master_table_iphop_filt %>%
  filter(str_detect(Locat, "Plume")) #%>%
  #filter(str_detect(Locat, "Plume")) #%>%
  #filter(str_detect(d, "Archaea")) %>%
  #filter(str_detect(p, "p__Pseudomonadota")) #%>%
  #filter(str_detect(p, "p__Methanobacteriota_B")) #%>%
  #filter(str_detect(p, "p__Nanoarchaeota"))
  #filter(str_detect(p, "p__Bacteroidota"))
  #filter(str_detect(c, "c__Gammaproteobacteria"))
  #filter(str_detect(c, "c__Alphaproteobacteria"))
  #filter(str_detect(p, "p__Campylobacterota"))
  #filter(str_detect(g, "g__Colwellia"))

#check how many viruses have >1 host prediction
master_table_iphop_filt3 <- iphop %>%
  group_by(Virus) %>% 
  filter(n()>1) #%>%
length(unique(master_table_iphop_filt3$Virus)) #614 viruses

####################################### vMAGs >4 seqs #############################################

# vMAG_5plus_prot <- master_table %>%
#   filter(str_detect(vMAG, "vRhyme")) %>%
#   filter(num_seqs > 4)
# #1,204 vMAGs have >4 seqs
# # write.table(vMAG_5plus_prot, file = "output/vMAGs_5plus_prots.tsv", col.names = TRUE,
# #             row.names = FALSE, quote = FALSE, sep = "\t")
# 
# vMAG_5plus <- master_table_noProtein %>%
#   filter(str_detect(vMAG, "vRhyme")) %>%
#   filter(num_seqs > 4)
# 
# #getting list of vMAGs with >10 scaffolds to return to unbinned fraction
# vMAG_11plus <- master_table_noProtein %>%
#   filter(num_seqs > 10) %>%
#   select(vMAG)
# # write.table(vMAG_11plus, file = "output/vMAG_11scaffold_plus.txt", col.names = FALSE,
# #             row.names = FALSE, sep = "\t", quote = FALSE)

############################# Generating supplementary figures #################################

############################### CheckV Quality ##################################

master_table_fig <- master_table_noProtein %>% 
  select(c("vMAG", "checkv_quality")) %>%
  #group_by(Site) %>%
  count(checkv_quality)

master_table_fig$checkv_quality <- factor(master_table_fig$checkv_quality, 
                                          levels = c('Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete'))
#level_order <- c('Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete') 

library(scales) #for label_comma()
###plot
dev.off()
p <- ggplot(master_table_fig, aes(x = factor(checkv_quality), #reorder bc was plotting x axis backwards/upside down, , levels = level_order
                                  y = n, fill = checkv_quality)) + 
  geom_bar(stat = "identity") + 
  #facet_wrap(~Site) + 
  xlab(element_blank())  +
  ylab("Count") +
  ggtitle("Virus Genome Quality") +
  scale_fill_viridis_d(name = "CheckV Quality", 
                       guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40000), 
                     breaks = seq(0, 40000, by = 10000), labels = label_comma()) +
  scale_x_discrete(expand = c(0, 0))
#geom_text(aes(label = paste0(n), y = n),
#          hjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 16)) +
  coord_flip()
p

ggsave("output/CheckV_quality.pdf", p, width = 15, height = 10, dpi = 500)
ggsave("output/CheckV_quality.png", p, dpi = 500, width = 10, height = 6) #width = 15, height = 10,

############################### Genome Size ##################################

master_table_noProtein$checkv_quality <- factor(master_table_noProtein$checkv_quality, 
                                                levels = c('Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'))

dev.off()
p <- ggplot(master_table_noProtein, aes(x = factor(checkv_quality), y = sum_len, fill = checkv_quality)) + 
  geom_boxplot() + 
  #facet_wrap(~Site) + 
  xlab("CheckV Quality")  +
  ylab("Genome Size") +
  ggtitle("Viral Genome Sizes") +
  guides(fill=guide_legend(title = element_blank())) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_x_discrete(limits = rev(levels(master_table_noProtein$checkv_quality)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 600000), breaks = seq(0, 600000, by = 100000), expand = c(0, 0),
                     labels = label_comma()) + #, expand = c(0, 0)
  scale_fill_viridis_d(guide = guide_legend(reverse = TRUE), alpha = 0.8) +
  #ylim(1000,565000) +
  coord_flip()
p

ggsave("output/GenomeSize.png", p, dpi = 500) #, width = 12, height = 6,
ggsave("output/GenomeSize.pdf", p, width = 12, dpi = 500)

############################### Virus taxonomy ##################################

master_table_fig <- master_table_noProtein %>% 
  select(c("vMAG", "r", "c")) %>%
  group_by(r) %>%
  count(c) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(r = gsub("\\br__\\b", "Unknown Realm", r)) %>%
  mutate(c = gsub("\\bc__\\b", "Unknown Class", c)) %>%
  mutate(r = gsub("r__", "", r)) %>%
  mutate(c = gsub("c__", "", c)) %>%
  #filter(!(c %in% c("Caudoviricetes"))) %>% #remove Caudos?
  #filter(!(c %in% c("Unknown Class"))) %>% #remove unknown class predictions
  arrange(r) %>%
  mutate(log_n = log(n) + 1)

# #add category for Caudo to plot separately because such an outlier
# master_table_fig$caudo <- ifelse(master_table_fig$n > 1000, "Caudo", "Other")
# master_table_fig <- master_table_fig %>%
#   mutate(caudo = ifelse(c == "Unknown Class", "Unknown", caudo))

#factor for plotting class in order of realm group
master_table_fig$r <- factor(master_table_fig$r, levels = unique(master_table_fig$r))
master_table_fig$c <- factor(master_table_fig$c, levels = unique(master_table_fig$c))
#aster_table_fig$caudo <- factor(master_table_fig$caudo, levels = unique(master_table_fig$caudo))

###plot
dev.off()
p <- ggplot(master_table_fig, aes(x = log_n, y = c, fill = r)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.2) + 
  #geom_bar(data = subset(master_table_fig, caudo == "Unknown"), stat = "identity", position = "dodge", width = 1) +  # Custom thickness for the facet where caudo == "Unknown"
  xlab(expression(paste("Number of genomes ", (log[10]+1))))  +
  ylab("Class") +
  ggtitle("Virus geNomad Taxonomy") +
  scale_fill_viridis_d(name = "Viral Realm", direction = -1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 12), 
                     breaks = seq(0, 12, by = 3)) +
  scale_y_discrete(expand = c(0, 0), limits = rev(levels(master_table_fig$c))) +  # Reverse the order of y-axis labels
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed"),
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), strip.background = element_blank()) #+
  #facet_wrap(~caudo, scales = "free", ncol = 1, strip.position = "left")
p

ggsave("output/genomad_tax.png", p, dpi = 500, width = 6, height = 4) #, width = 12, height = 6,
ggsave("output/genomad_tax.pdf", p, dpi = 500, width = 6, height = 4)

############################### Host predictions ##################################
#add site info of vMAG to iphop table
master_table_iphop$vMAG_Site <- master_table_iphop$vMAG
master_table_iphop <- master_table_iphop %>%
  separate(vMAG_Site, c("vMAG_Site", NA), sep = "_vRhyme|_NODE|_k95|_scaffold")
master_table_iphop$vMAG_Site <- stri_replace_all_regex(master_table_iphop$vMAG_Site,
                                                       pattern=c("_A[0-9]",
                                                                 "_T[0-9][0-9]", "_T[0-9]", "_S0[0-9][0-9]",
                                                                 "_S1[0-9][0-9]", "_[0-9][0-9][0-9]-[0-9][0-9][0-9]",
                                                                 "-38[0-9]"),
                                                       replacement='',
                                                       vectorize=FALSE)

master_table_iphop$vMAG_Site <- gsub("*_M1[0-9]","",master_table_iphop$vMAG_Site)
master_table_iphop$vMAG_Site <- gsub("*_M[0-9]","",master_table_iphop$vMAG_Site)

write.table(master_table_iphop, file = "output/master_table_iphop.tsv", quote = FALSE,
            row.names = FALSE, sep = "\t")

#get just Proteobacteria to highlight their breakdown
master_table_iphop_proteos <- master_table_iphop %>%
  select(vMAG, Host.genus) %>%
  separate(Host.genus, c('d', 'p', 'c', 'o', 'f', 'g'), sep= ";") %>%
  select(vMAG, d, p, c) %>%
  filter(grepl("Pseudomonadota", p)) %>%
  group_by(d, c) %>%
  count(c) %>%
  ungroup() %>%
  mutate(c = str_replace(c, "c__", "")) %>%
  mutate(d = str_replace(d, "d__", "")) %>%
  mutate(c = str_replace(c, "$", "*")) %>%
  rename("p" = "c") #arbitrary rename for merging
  
#everything else minus Proteos
master_table_fig <- master_table_iphop %>%
  select(vMAG, Host.genus) %>%
  separate(Host.genus, c('d', 'p', 'c', 'o', 'f', 'g'), sep= ";") %>%
  select(vMAG, d, p, c) %>% #, c
  group_by(d, p) %>%
  count(p) %>%
  ungroup() %>%
  drop_na() %>%
  #mutate(c = str_replace(c, "c__", "")) %>%
  mutate(p = str_replace(p, "p__", "")) %>%
  mutate(d = str_replace(d, "d__", "")) %>%
  mutate_at(vars(p), ~na_if(., "")) %>%
  drop_na() %>%
  filter(!grepl("Pseudomonadota", p))

#bind the proteos and whole count together
master_table_fig <- rbind(master_table_fig, master_table_iphop_proteos)

#remove host predictions less than 10 for clarity of plot
master_table_fig <- master_table_fig %>%
  filter(n > 10)

###plot
dev.off()
p <- ggplot(master_table_fig, aes(x = n, 
                                  y = factor(p, levels = rev(levels(factor(p)))), 
                                  fill = p)) + #, fill = p
  geom_bar(position = "dodge", stat = "identity", width = 0.2) + 
  #geom_bar(data = subset(master_table_fig, caudo == "Unknown"), stat = "identity", position = "dodge", width = 1) +  # Custom thickness for the facet where caudo == "Unknown"
  xlab("Number of host predictions")  +
  ylab("Microbial Host Phyla") +
  ggtitle("") +
  scale_fill_viridis_d(guide="none") + #, direction = -1, name = "Microbial Phyla"
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +  # Reverse the order of y-axis labels
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ d, scales = "free") #ncol = 1, strip.position = "left"
p

#ggsave("output/iphop_hosts.png", p, dpi = 500, width = 10, height = 5) #, width = 12, height = 6,
#ggsave("output/iphop_hosts.pdf", p, dpi = 500, width = 6, height = 4)


##############################    Version of plot where faceting by site

#get just Proteobacteria to highlight their breakdown
master_table_iphop_proteos <- master_table_iphop %>%
  select(vMAG, Host.genus, vMAG_Site) %>%
  separate(Host.genus, c('d', 'p', 'c', 'o', 'f', 'g'), sep= ";") %>%
  select(vMAG, d, p, c, vMAG_Site) %>%
  filter(grepl("Pseudomonadota", p)) %>%
  group_by(d, c, vMAG_Site) %>%
  count(c) %>%
  ungroup() %>%
  mutate(c = str_replace(c, "c__", "")) %>%
  mutate(d = str_replace(d, "d__", "")) %>%
  mutate(c = str_replace(c, "$", "*")) %>%
  rename("p" = "c") #arbitrary rename for merging

#everything else minus Proteos
master_table_fig <- master_table_iphop %>%
  select(vMAG, Host.genus, vMAG_Site) %>%
  separate(Host.genus, c('d', 'p', 'c', 'o', 'f', 'g'), sep= ";") %>%
  select(vMAG, d, p, c, vMAG_Site) %>% #, c
  group_by(d, p, vMAG_Site) %>%
  count(p) %>%
  ungroup() %>%
  drop_na() %>%
  #mutate(c = str_replace(c, "c__", "")) %>%
  mutate(p = str_replace(p, "p__", "")) %>%
  mutate(d = str_replace(d, "d__", "")) %>%
  mutate_at(vars(p), ~na_if(., "")) %>%
  drop_na() %>%
  filter(!grepl("Pseudomonadota", p))

#bind the proteos and whole count together
master_table_fig <- rbind(master_table_fig, master_table_iphop_proteos)

#remove host predictions less than 10 for clarity of plot
master_table_fig <- master_table_fig %>%
  filter(n > 10)

###plot
dev.off()
p <- ggplot(master_table_fig, aes(x = n, 
                                  y = factor(p, levels = rev(levels(factor(p)))), 
                                  fill = p)) + #, fill = p
  geom_bar(position = "dodge", stat = "identity", width = 0.2) + 
  #geom_bar(data = subset(master_table_fig, caudo == "Unknown"), stat = "identity", position = "dodge", width = 1) +  # Custom thickness for the facet where caudo == "Unknown"
  xlab("Number of host predictions")  +
  ylab("Microbial Host Phyla") +
  ggtitle("") +
  scale_fill_viridis_d(guide="none") + #, direction = -1, name = "Microbial Phyla"
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +  # Reverse the order of y-axis labels
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ vMAG_Site, scales = "free") #ncol = 1, strip.position = "left"
p

#ggsave("output/iphop_hosts.png", p, dpi = 500, width = 10, height = 5) #, width = 12, height = 6,
#ggsave("output/iphop_hosts.pdf", p, dpi = 500, width = 6, height = 4)


