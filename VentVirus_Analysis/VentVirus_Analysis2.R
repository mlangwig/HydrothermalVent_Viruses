##################### Integrating all virus output into 1 table #####################

setwd(dir = "~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/VentVirus_Analysis/")

library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)

########################################## Read in major inputs ##########################################

#iphop

checkv <- read.delim(file = "input/quality_summary.tsv", header = TRUE)
#remove cols
checkv <- checkv %>%
  select(-c(contig_length, provirus, proviral_length, completeness_method, kmer_freq))
table(checkv$checkv_quality)

vib_type <- read.delim(file = "input/52samples_VIBRANT_type.tsv", header = TRUE)
vib_type$scaffold <- gsub("=","_",vib_type$scaffold)

gensize <- read.table(file = "input/PlumeVentVirus_Seqkit_final.txt", sep = " ", header = TRUE)
#count range of 1-5kb genome sizes
count <- sum(gensize$sum_len >= 1000 & gensize$sum_len <= 5000)

vib_annos <- read.delim(file = "../../VIBRANT_annos/52samples_VIBRANT_annos.tsv", sep = "\t", fill = TRUE, 
                        header = TRUE, na.strings=c("","NA"))
#replace = with _
vib_annos$protein <- gsub("=","_",vib_annos$protein)
vib_annos$scaffold <- gsub("=","_",vib_annos$scaffold)

vMAG_mapping <- read.table(file = "input/vMAG_scaffolds_6088_final.txt", sep = "\t")
#rename unnamed cols
vMAG_mapping <- vMAG_mapping %>%
  rename("vMAG" = "V1",
         "scaffold" = "V2") #new = old

genomad_tax <- read.table(file = "input/43995_final_VentViruses_Nlinked_taxonomy_parsed.tsv", header = TRUE)
#remove cols
genomad_tax <- genomad_tax %>%
  select(-c(n_genes_with_taxonomy, agreement, taxid)) %>%
  separate(lineage, c("d", "p", "c", "o", "f", "g"), 
           sep= ";")
table(genomad_tax$d)
table(genomad_tax$o)
  
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
#640 lysogenic vMAGs, 5,448 lytic

#get VIBRANT types of non-vMAGs
vib_type_vUnbinned <- vib_type %>%
  select(c('scaffold', 'type')) %>%
  unique()
`%notin%` <- Negate(`%in%`)
vib_type_vUnbinned<-vib_type_vUnbinned[vib_type_vUnbinned$scaffold %notin% vMAG_mapping$scaffold,]
#should be 37,907 or the number of unbinned viruses
vib_type_vUnbinned <- vib_type_vUnbinned %>%
  select(c('type','scaffold')) %>%
  rename("vMAG" = "scaffold")
#make final VIBRANT type
vib_type <- rbind(vib_type_vMAG, vib_type_vUnbinned)
table(vib_type$type)
#2,391 lysogenic viruses, 41,604 lytic

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
#43,995 viruses, same as total

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
  select(vMAG, protein:warnings, num_seqs:sum_len, d:g)

#add lifestyle type
master_table <- vib_type %>%
  right_join(master_table, by = c("vMAG" = "vMAG"))

#add iphop host prediction
#master_table_iphop <- iphop %>%
# right_join(master_table, by = c("" = ""))


#make master table without proteins for accurate counting
master_table_noProtein <- master_table %>%
  select(vMAG, type, gene_count:g) %>%
  unique()

#master table no proteins of med-quality and better only
master_table_noProtein_hq <- master_table %>%
  select(vMAG, type, gene_count:g) %>%
  unique() %>%
  filter(!(checkv_quality %in% c("Low-quality", "Not-determined")))
table(master_table_noProtein_hq$type)
#395 lysogenic, 1686 lytic for med quality and better viruses

############################# Generating supplementary figures #################################

######### CheckV quality

master_table_fig <- master_table_noProtein %>% 
  select(c("vMAG", "checkv_quality")) %>%
  #group_by(Site) %>%
  count(checkv_quality)

master_table_fig$checkv_quality <- factor(master_table_fig$checkv_quality, 
                                          levels = c('Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete'))
#level_order <- c('Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete') 

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
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))
  #scale_y_continuous(breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000))
#geom_text(aes(label = paste0(n), y = n),
#          hjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 16)) +
  coord_flip()
p

ggsave("output/CheckV_quality.pdf", p, width = 15, height = 10, dpi = 500)
ggsave("output/CheckV_quality.png", p, width = 15, height = 10, dpi = 500)

######### Genome Size

master_table_noProtein$checkv_quality <- factor(master_table_noProtein$checkv_quality, 
                                                levels = c('Complete', 'High-quality', 'Medium-quality', 'Not-determined', 'Low-quality'))

dev.off()
p <- ggplot(master_table_noProtein, aes(x = factor(checkv_quality), y = sum_len, fill = checkv_quality)) + 
  geom_boxplot() + 
  #facet_wrap(~Site) + 
  xlab("Site")  +
  ylab("Genome Size") +
  ggtitle("Viral Genome Sizes") +
  guides(fill=guide_legend(title = "CheckV Quality")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_x_discrete(limits = rev(levels(master_table_noProtein$checkv_quality)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 600000), breaks = seq(0, 600000, by = 100000), expand = c(0, 0)) + #, expand = c(0, 0)
  scale_fill_viridis_d(guide = guide_legend(reverse = TRUE)) +
  #ylim(1000,565000) +
  coord_flip()
p

ggsave("output/GenomeSizeKB.png", p, width = 12, height = 6, dpi = 500)
ggsave("output/GenomeSizeKB.pdf", p, width = 12, dpi = 500)









