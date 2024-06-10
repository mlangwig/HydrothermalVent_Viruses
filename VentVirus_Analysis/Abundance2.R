########################## Compiling and parsing virus CoverM abundance ############################
library(dplyr)
library(tidyverse)
library(knitr)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(viridis)
library(pals)
library(RColorBrewer)
library(stringi)
library(ggbreak)

######################################### inputs ##############################################

#coverm abundance with read counts and genome coverage
abund <- read.delim2("input/PlumeVentVirus-vs-Reads-CoverM-Count_MinCov.tsv")

#seqkit number of reads per sample
read_length <- read.delim2(file = "../../abundance/stats.tsv")
read_length$file <- gsub("-",".",read_length$file) #replace symbol

#metadata to map site names to reads
abun_names <- read.delim2(file = "input/AssembliesToReads_Mapping.txt", header = FALSE)
abun_names <- abun_names %>%
  rename("Site" = "V1") %>%
  rename("Reads" = "V2")
abun_names$Reads <- gsub("-",".",abun_names$Reads) #replace symbol

######################################### parse ##############################################

#remove unmapped
abund <- abund[-1,]

#convert 0s to NAs so they are not included in distribution assessments
abund <- abund %>% dplyr::mutate_all(funs(ifelse(. == 0, NA, .)))

#convert from wide to long
abund_long <- melt(abund, id = c("Genome"), na.rm = TRUE)
#make a new column separating out the number type
abund_long <- separate(abund_long, 
                       variable, into = c("variable", "Type"), sep = "(?<=\\.fastq\\.gz\\.)")
abund_long$variable <- gsub(".fastq.gz.",".fastq.gz",abund_long$variable)

#map the site names to the reads
abund_long <- abun_names %>%
  dplyr::select("Site", "Reads") %>%
  right_join(abund_long, by = c("Reads" = "variable"))

#now group by genome and variable and keep rows where covered fraction is >=.7
abund_long_cov <- abund_long %>%
  group_by(Genome, Reads) %>%
  filter(any(Type == "Covered.Fraction" & value >= 0.7) | !any(Type == "Covered.Fraction")) %>%
  ungroup() %>%
  group_by(Reads, Type) %>%
  mutate(Percentile_Rank=rank(as.numeric(value))/length(as.numeric(value))) %>% #add col with percentile rank per category
  ungroup()

#map number of reads
abund_long_cov <- read_length %>%
  dplyr::select(c("file", "num_seqs")) %>%
  right_join(abund_long_cov, by = c("file" = "Reads")) %>%
  rename("Reads" = "file")
abund_long_cov$value <- as.numeric(abund_long_cov$value)

#divide read count by number of reads *100
abund_long_norm <- abund_long_cov
abund_long_norm <- abund_long_norm %>%
  mutate(abun_norm = (value/num_seqs*100)) %>%
  filter(Type == "Read.Count")
  #filter(Type == "Relative.Abundance....")

# ### with Patricia
# abund_long_norm_Lau_TM <- abund_long_norm %>%
#   filter(Site == "Lau_Basin_Tahi_Moana_min1000_2")

##################################### host metadata ##############################################

abund_long_norm_iphop <- iphop %>%
  dplyr::select(Virus, "Host.genus") %>%
  right_join(abund_long_norm, by = c("Virus" = "Genome"))   
#only keep phylum and class of the tax string
abund_long_norm_iphop <- abund_long_norm_iphop %>% separate("Host.genus", c("d", "p", "c", "o", "f", "g"), 
                                                            sep= ";")
#select only phylum and class taxonomy
abund_long_norm_iphop <- abund_long_norm_iphop %>% select(c("Virus", "Site", "Reads", "p", "c", "o", "f", "g","abun_norm"))
abund_long_norm_iphop <- abund_long_norm_iphop %>% 
  filter(str_detect(c, "c__Gammaproteobacteria") | str_detect(p, "p__Campylobacterota"))

#abund_long_norm_iphop <- abund_long_norm_iphop %>% 
#  filter(str_detect(c, "c__Gammaproteobacteria"))

abund_long_norm_iphop$Site <- gsub("min1000","Plume", abund_long_norm_iphop$Site) 
abund_long_norm_iphop$Site <- gsub("_scaffolds_Plume","", abund_long_norm_iphop$Site) 
#abun_long_iphop_p$Site <- gsub("Seawater_scaffolds_Plume","Seawater", abun_long_iphop_p$Site) 

#add column for plume vs vent
abund_long_norm_iphop$Locat <- abund_long_norm_iphop$Site
abund_long_norm_iphop$Locat <- gsub(".*Plume.*","Plume", abund_long_norm_iphop$Locat)
abund_long_norm_iphop$Locat <- gsub(".*Seawater.*","Plume", abund_long_norm_iphop$Locat)
abund_long_norm_iphop$Locat <- gsub(".*Brothers.*","Vent", abund_long_norm_iphop$Locat)
abund_long_norm_iphop$Locat <- gsub(".*ELSC.*","Vent", abund_long_norm_iphop$Locat)
abund_long_norm_iphop$Locat <- gsub(".*EPR.*","Vent", abund_long_norm_iphop$Locat)
abund_long_norm_iphop$Locat <- gsub(".*Guaymas.*","Vent", abund_long_norm_iphop$Locat)
abund_long_norm_iphop$Locat <- gsub(".*MAR.*","Vent", abund_long_norm_iphop$Locat)

#filter for viruses with viral hallmark
abund_long_norm_iphop_vh <- master_table_iphop_filt %>%
  dplyr::select(vMAG, total_hallmarks, d:g) %>%
  right_join(abund_long_norm_iphop, by = c("vMAG" = "Virus"))

#abun_long_iphop <- abun_long_iphop %>% 
#  filter(str_detect(Locat, "Plume")) #%>%
# filter(str_detect(o, "o__PS1"))

# #for Proteobacteria keep class, everything else, keep phylum
abund_long_proteo <- abund_long_norm_iphop %>% filter(grepl("p__Pseudomonadota", p))
abund_long_proteo <- abund_long_proteo %>% select(-c("p"))
abund_long_proteo <- abund_long_proteo %>% rename("Taxa" = "c")
abund_long_proteo <- abund_long_proteo %>% select(c("Virus","Taxa","Site","abun_norm", "Locat"))

abund_long_norm_iphop <- abund_long_norm_iphop %>% filter(!grepl("p__Pseudomonadota", p))
abund_long_norm_iphop <- abund_long_norm_iphop %>% select(-c("c"))
abund_long_norm_iphop <- abund_long_norm_iphop %>% rename("Taxa" = "p")
abund_long_norm_iphop <- abund_long_norm_iphop %>% select(c("Virus","Taxa","Site","abun_norm", "Locat"))
#put the data frames back together
abund_long_norm_iphop <- rbind(abund_long_norm_iphop, abund_long_proteo)

############################ Sum abundance by predicted host and site #############################

abund_long_norm_iphop_p <- abund_long_norm_iphop %>%
  group_by(Site, Taxa, Locat) %>% 
  summarise(value=sum(as.numeric(abun_norm))) %>% #summing the group
  #filter(grepl("c__Gammaproteobacteria|p__Campylobacterota", Taxa)) %>% #only grab Gamma and Campylo
  ungroup()

############################ Set order of sites for plotting #############################

# test_Camp <- abund_long_norm_iphop_p %>%
#   filter(grepl("p__Campylobacterota", Taxa)) %>%
#   group_by(Taxa) %>%
#   arrange((value))
# #cat(dQuote(test_Camp$Site, FALSE), '\n', sep = ",") #to print strings 
# 
# test_Gam <- abund_long_norm_iphop_p %>%
#   filter(!grepl("p__Campylobacterota", Taxa)) %>%
#   group_by(Taxa) %>%
#   arrange(desc(value))
#   
# #put the data frames back together
# abun_long_iphop_p <- rbind(abun_iphop_p_Camp, abun_iphop_p_Gam)
#     
# # lock in factor level order
# #abun_long_iphop_p$value <- factor(abun_long_iphop_p$value, levels = abun_long_iphop_p$value)
# abun_long_iphop_p$Site <- factor(abun_long_iphop_p$Site, levels = unique(abun_long_iphop_p$Site))

############ Change names #############################

#sulfur_VentVirus_AMGs$KO <- gsub("K23144","UAP", sulfur_VentVirus_AMGs$KO)

abund_long_norm_iphop_p$Taxa <- gsub("c__|p__","",abund_long_norm_iphop_p$Taxa)
abund_long_norm_iphop_p$Site <- gsub("_"," ",abund_long_norm_iphop_p$Site)
abund_long_norm_iphop_p$Locat <- gsub("Vent","Deposit",abund_long_norm_iphop_p$Locat)


#set order of x axis 
abund_long_norm_iphop_p$Site <- factor(abund_long_norm_iphop_p$Site, 
                                       levels=c("Axial Plume", "MAR Rainbow 354-166",
                                                "Cayman Shallow Plume 1", "Lau Basin Kilo Moana Plume 3",
                                                "Cayman Shallow Plume 3", "Cayman Shallow Plume 2",
                                                "Lau Basin Tahi Moana Plume 1", "Lau Basin Abe Plume 1",
                                                "Brothers NWCB S139",
                                                "Cayman Deep Plume 2", "Lau Basin Kilo Moana Plume 2",
                                                "Brothers NWCB S012",
                                                "Cayman Deep Plume 3", "Lau Basin Mariner Plume 1",
                                                "Cayman Deep Plume 1", "Brothers NWCB S140",
                                                "Lau Basin Kilo Moana Plume 1",
                                                "Brothers NWCB S146", 
                                                "ELSC Tui Malila 134-614", "Brothers NWCB S141",
                                                "Lau Basin Abe Plume 3", "EPR 4281-140",
                                                "ELSC Abe A3", "MAR Lucky 356-308",
                                                "Lau Basin Mariner Plume 2", "Axial Seawater",
                                                "Guaymas Basin Plume", "Lau Basin Tahi Moana Plume 2",
                                                "Lau Basin Kilo Moana Plume 4", 
                                                "ELSC Abe 128-326",
                                                "EPR PIR-30", "ELSC Vai Lili V2",
                                                "Lau Basin Tui Malila Plume", "Lau Basin Abe Plume 2", 
                                                #end Gamma
                                                "Brothers Diffuse S015", "Brothers Diffuse S009",
                                                "MAR Lucky 356-284", "Brothers NWCA S143",
                                                "ELSC Abe A1", 
                                                "Brothers NWCA S144", "ELSC Tui Malila T2",
                                                "ELSC Tui Malila T11", "Brothers LC S014",
                                                "MAR Rainbow 355-202",
                                                "Guaymas 4559-240", "ELSC Tui Malila T10",
                                                "Brothers LC S016", "Guaymas 4571-419",
                                                "ELSC Mariner 131-447", "ELSC Mariner M17",
                                                "ELSC Bowl M2", "Brothers NWCA S013",
                                                "Brothers UC S147", "ELSC Tui Malila 132-544",
                                                "ELSC Mariner M10", "Brothers NWCA S017",
                                                "ELSC Bowl M1",
                                                "Guaymas 4561-380","Brothers UC S010",
                                                "Brothers NWCA S145",
                                                "Brothers UC S011","Brothers NWCA S142","Guaymas 4561-384")) # start Camp lowest to highest 


####################### Stacked bar plot of Gamma and Campylo-infecting viral abundance ################################
## Shows functional redundancy of viruses infecting these hosts

# #distinct colors
# library(RColorBrewer)
# n <- 39
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# #pie(rep(1,n), col=sample(col_vector, n))

#factor
level_order <- c('Gammaproteobacteria', 'Campylobacterota') 

####The following produces Figure X, which was modified in Biorender
dev.off()
plot <- abund_long_norm_iphop_p %>%
  ggplot(aes(x = as.numeric(value), y = Site, fill = factor(Taxa, levels = c(level_order)))) + #y = reorder(Site, value, sum) | factor(checkv_quality, levels = level_order)
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(begin = .5,
                       end = 0) +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Virus Abundance", y = "Site",
       fill = "Predicted microbial host") +
  guides(fill=guide_legend(override.aes = list(size=8))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        strip.text.x = element_text(size = 12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_x_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  #geom_hline(yintercept = 21.5) +
  # annotate(geom="text", x=6, y=30, label="Plume",
  #          color="black") +
  #ggtitle("Viral Abundance") + #Change for top X grabbed
  facet_wrap(.~Locat) +
  scale_y_discrete(limits=rev)
#coord_flip()
plot

ggsave("output/coverm_CampGamma_normAbun.png", plot,
       height = 12, width = 15,
       bg = "transparent")

############################## Bubble plot of virus abundance summed by taxa by site ################################
## Shows distribution of different viral taxa at different sites

#Add genome quality info to remove viruses based on length or quality/both?

#prepare input data from the normalized abundance info
abund_long_norm_tax <- abund_long_norm %>%
  filter(Genome != "Lau_Basin_Tahi_Moana_vRhyme_bin_115")
abund_long_norm_tax <- genomad_tax %>%
  right_join(abund_long_norm_tax, by = c("genome" = "Genome")) %>%
  drop_na() %>%
  select(c("genome", "c", "abun_norm"))
abund_long_norm_tax$Site <- abund_long_norm_tax$genome
abund_long_norm_tax <- abund_long_norm_tax %>%
  separate(Site, c("Site", NA), sep = "_vRhyme|_NODE|_k95|_scaffold")
abund_long_norm_tax$Site_Gen <- abund_long_norm_tax$Site
abund_long_norm_tax$Site_Gen <- stri_replace_all_regex(abund_long_norm_tax$Site_Gen,
                                                pattern=c("_A[0-9]",
                                                          "_T[0-9][0-9]", "_T[0-9]", "_S0[0-9][0-9]",
                                                          "_S1[0-9][0-9]", "_[0-9][0-9][0-9]-[0-9][0-9][0-9]",
                                                          "-38[0-9]"),
                                                replacement='',
                                                vectorize=FALSE)

abund_long_norm_tax$Site_Gen <- gsub("*_M1[0-9]","",abund_long_norm_tax$Site_Gen)
abund_long_norm_tax$Site_Gen <- gsub("*_M[0-9]","",abund_long_norm_tax$Site_Gen)

abund_long_norm_tax <- abund_long_norm_tax %>%
  group_by(Site_Gen, c) %>%
  summarise(abun_norm = sum(abun_norm)) %>%
  ungroup() %>%
  mutate(log_n = log(abun_norm)) %>% #create log transformed abundance
  mutate(log_n_plus1 = log(abun_norm+1)) %>% #add 1 to log transformed #log_n+1
  mutate(c = gsub("c__", "", c)) %>% #remove c__
  mutate(c = sub("^$", "Unknown", c)) %>% #replace blanks with unknown
  mutate(Site_Gen = gsub("_", " ", Site_Gen)) %>% #remove underscore from site names
  mutate(Site_Gen = gsub("Guaymas Basin", "Guaymas Basin Plume", Site_Gen)) %>%
  mutate(Site_Gen = gsub(" V2", "", Site_Gen))
  
#set order of Sites on x axis
abund_long_norm_tax$Site_Gen <- factor(abund_long_norm_tax$Site_Gen, 
                                       levels=c("Axial Plume", "Axial Seawater",
                                                "Cayman Deep", "Cayman Shallow",
                                                "Guaymas Basin Plume", 
                                                "Lau Basin Abe", "Lau Basin Kilo Moana",
                                                "Lau Basin Mariner", "Lau Basin Tahi Moana", 
                                                "Lau Basin Tui Malila", 
                                                #end plume
                                                "Brothers Diffuse", "Brothers LC",
                                                "Brothers NWCA", "Brothers NWCB",
                                                "Brothers UC",
                                                "ELSC Abe", "ELSC Bowl",
                                                "ELSC Mariner", "ELSC Tui Malila",
                                                "ELSC Vai Lili",
                                                "EPR 4281-140", "EPR PIR-30",
                                                "Guaymas 4559-240", "Guaymas 4561",
                                                "Guaymas 4571-419",
                                                "MAR Lucky", "MAR Rainbow")) #end deposit

#add columns with metadata for faceting
abund_long_norm_tax$Site_Type <- abund_long_norm_tax$Site_Gen
abund_long_norm_tax$Site_Type <- gsub(".*Axial.*","Plume", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*Cayman.*","Plume", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub("Guaymas Basin Plume","Plume", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*Lau.*","Plume", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*Brothers.*","Deposit", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*ELSC.*","Deposit", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*EPR.*","Deposit", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*Guaymas.*","Deposit", abund_long_norm_tax$Site_Type)
abund_long_norm_tax$Site_Type <- gsub(".*MAR.*","Deposit", abund_long_norm_tax$Site_Type)

abund_long_norm_tax$Virus_Realm <- abund_long_norm_tax$c
abund_long_norm_tax$Virus_Realm <- gsub(".*Caudo.*","Duplodnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub(".*Malgrand.*","Monodnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Megaviricetes","Varidnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Faserviricetes","Monodnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Laserviricetes","Varidnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Tokiviricetes","Adnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Tectiliviricetes","Varidnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Arfiviricetes","Monodnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Herviviricetes","Duplodnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Duplopiviricetes","Riboviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Huolimaviricetes","Monodnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Pokkesviricetes","Varidnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Resentoviricetes","Riboviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Vidaverviricetes","Riboviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Polintoviricetes","Varidnaviria", abund_long_norm_tax$Virus_Realm)
abund_long_norm_tax$Virus_Realm <- gsub("Maveriviricetes","Varidnaviria", abund_long_norm_tax$Virus_Realm)

# write.table(abund_long_norm_tax, file = "output/abund_long_norm_tax.tsv", quote = FALSE,
#             row.names = FALSE, sep = "\t")

dev.off()
p <- ggplot(abund_long_norm_tax, aes(y=c, x=Site_Gen))+
  geom_point(aes(size=log_n_plus1, color = log_n_plus1))+ 
  #scale_size_continuous(breaks = c(-8, -4, 0, 4))+
  #scale_color_continuous(guide="legend", 
  #                       type = "viridis")+ #,breaks = c(-8, -4, 0, 4)
  scale_color_viridis_c(option = "G",
                        guide = "legend")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=-90),
        #panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-.025),
        text = element_text(color="black"),
        legend.position="right",
        strip.placement = "outside",
        panel.spacing.x = unit(0, "lines"),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+ #to prevent facet from being very spaced out
  scale_x_discrete(limits=rev)+ #position="top"
  scale_y_discrete(position = "right")+ #confusing it is right because of coord_flip
  xlab("")+ #sites
  ylab("")+ #viral class
  labs(size = "Viral log relative abundance",
       color = "Viral log relative abundance") +
  facet_grid(Site_Type ~ Virus_Realm, scales = "free", switch = "y") + #wow this took me awhile - remember you can't factor the y axis and then get it to facet properly
  #scale_y_discrete(limits=rev)+ # has to be this instead of factor
  coord_flip()
p  


ggsave("output/VirusAbundance_byTaxa2.png", p, dpi = 500, width = 12, height = 7) #, width = 12, height = 6,
ggsave("output/VirusAbundance_byTaxa2.pdf", p, dpi = 500, width = 12, height = 7)


########################## Compiling and parsing microbial MAG CoverM abundance ############################

######################################### inputs ##############################################

#coverm abundance with read counts and genome coverage
abund_MAGs <- read.delim2("input/PlumeVentMicrobialMAGs-vs-Reads-CoverM-Count_MinCov.tsv")

#metadata to map microbial taxonomy
gtdb <- read.delim2(file = "input/gtdbtkv2_3_2_PlumeVentMAGs.txt", header = TRUE)
gtdb <- gtdb %>%
  select(c("user_genome" , "classification")) %>%
  separate(classification, c('d', 'p', 'c', 'o', 'f', 'g', 's'), sep= ";")

######################################### parse ##############################################

#remove unmapped
abund_MAGs <- abund_MAGs[-1,]

#convert 0s to NAs so they are not included in distribution assessments
abund_MAGs <- abund_MAGs %>% dplyr::mutate_all(funs(ifelse(. == 0, NA, .)))

#convert from wide to long
abund_MAGs_long <- melt(abund_MAGs, id = c("Genome"), na.rm = TRUE)
#make a new column separating out the number type
abund_MAGs_long <- separate(abund_MAGs_long, 
                       variable, into = c("variable", "Type"), sep = "(?<=\\.fastq\\.gz\\.)")
abund_MAGs_long$variable <- gsub(".fastq.gz.",".fastq.gz",abund_MAGs_long$variable)

#map the site names to the reads
abund_MAGs_long <- abun_names %>%
  dplyr::select("Site", "Reads") %>%
  right_join(abund_MAGs_long, by = c("Reads" = "variable"))

#now group by genome and variable and keep rows where covered fraction is >=.7
abund_MAGs_long_cov <- abund_MAGs_long %>%
  group_by(Genome, Reads) %>%
  filter(any(Type == "Covered.Fraction" & value >= 0.7) | !any(Type == "Covered.Fraction")) %>%
  ungroup() %>%
  group_by(Reads, Type) %>%
  mutate(Percentile_Rank=rank(as.numeric(value))/length(as.numeric(value))) %>% #add col with percentile rank per category
  ungroup()

#map number of reads
abund_MAGs_long_cov <- read_length %>%
  dplyr::select(c("file", "num_seqs")) %>%
  right_join(abund_MAGs_long_cov, by = c("file" = "Reads")) %>%
  rename("Reads" = "file")
abund_MAGs_long_cov$value <- as.numeric(abund_MAGs_long_cov$value)

#divide read count by number of reads *100
abund_MAGs_long_norm <- abund_MAGs_long_cov
abund_MAGs_long_norm <- abund_MAGs_long_norm %>%
  mutate(abun_norm = (value/num_seqs*100)) %>%
  filter(Type == "Read.Count")

##################################### MAG taxonomy metadata ##############################################

abund_MAGs_long_norm_gtdb <- gtdb %>%
  dplyr::select(user_genome:s) %>%
  right_join(abund_MAGs_long_norm, by = c("user_genome" = "Genome"))   

#select only phylum and class taxonomy
#abund_MAGs_long_norm_gtdb <- abund_MAGs_long_norm_gtdb %>% select(c("user_genome", "Site", "Reads", d:s,"abun_norm"))
abund_MAGs_long_norm_gtdb <- abund_MAGs_long_norm_gtdb %>% 
  filter(str_detect(c, "c__Gammaproteobacteria") | str_detect(p, "p__Campylobacterota"))

abund_MAGs_long_norm_gtdb$Site <- gsub("_metaspades_scaffolds.min1000.fasta","", abund_MAGs_long_norm_gtdb$Site) 
abund_MAGs_long_norm_gtdb$Site <- gsub("min1000","Plume", abund_MAGs_long_norm_gtdb$Site) 
abund_MAGs_long_norm_gtdb$Site <- gsub("_scaffolds_Plume","", abund_MAGs_long_norm_gtdb$Site) 

#add column for plume vs vent
abund_MAGs_long_norm_gtdb$Locat <- abund_MAGs_long_norm_gtdb$Site
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*Plume.*","Plume", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*Seawater.*","Plume", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*Cayman.*","Plume", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*Brothers.*","Vent", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*ELSC.*","Vent", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*EPR.*","Vent", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*Guaymas.*","Vent", abund_MAGs_long_norm_gtdb$Locat)
abund_MAGs_long_norm_gtdb$Locat <- gsub(".*MAR.*","Vent", abund_MAGs_long_norm_gtdb$Locat)

############################ Sum abundance by predicted host and site #############################

#for Proteobacteria/Pseudomonadota keep class, everything else, keep phylum
abund_MAGs_long_proteo <- abund_MAGs_long_norm_gtdb %>% filter(grepl("p__Pseudomonadota", p))
abund_MAGs_long_proteo <- abund_MAGs_long_proteo %>% select(-c("p"))
abund_MAGs_long_proteo <- abund_MAGs_long_proteo %>% rename("Taxa" = "c")
abund_MAGs_long_proteo <- abund_MAGs_long_proteo %>% select(c("user_genome","Taxa","Site","abun_norm", "Locat"))

abund_MAGs_long_norm_gtdb <- abund_MAGs_long_norm_gtdb %>% filter(!grepl("p__Pseudomonadota", p))
abund_MAGs_long_norm_gtdb <- abund_MAGs_long_norm_gtdb %>% select(-c("c"))
abund_MAGs_long_norm_gtdb <- abund_MAGs_long_norm_gtdb %>% rename("Taxa" = "p")
abund_MAGs_long_norm_gtdb <- abund_MAGs_long_norm_gtdb %>% select(c("user_genome","Taxa","Site","abun_norm", "Locat"))
#put the data frames back together
abund_MAGs_long_norm_gtdb <- rbind(abund_MAGs_long_norm_gtdb, abund_MAGs_long_proteo)

abund_MAGs_long_norm_gtdb_p <- abund_MAGs_long_norm_gtdb %>%
  group_by(Site, Taxa, Locat) %>% 
  summarise(value=sum(as.numeric(abun_norm))) %>% #summing the group
  #filter(grepl("c__Gammaproteobacteria|p__Campylobacterota", Taxa)) %>% #only grab Gamma and Campylo
  ungroup()

# tst <- abund_MAGs_long_norm_gtdb %>%
#   filter(Site == "Brothers_UC_S011",
#          Taxa == "p__Campylobacterota")

abund_MAGs_long_norm_gtdb_p$Taxa <- gsub("c__|p__","",abund_MAGs_long_norm_gtdb_p$Taxa)
abund_MAGs_long_norm_gtdb_p$Site <- gsub("_"," ",abund_MAGs_long_norm_gtdb_p$Site)
abund_MAGs_long_norm_gtdb_p$Locat <- gsub("Vent","Deposit",abund_MAGs_long_norm_gtdb_p$Locat)
abund_MAGs_long_norm_gtdb_p$Site <- gsub(".fasta","",abund_MAGs_long_norm_gtdb_p$Site)

#set order of x axis 
abund_MAGs_long_norm_gtdb_p$Site <- factor(abund_MAGs_long_norm_gtdb_p$Site, 
                                       levels=c("Axial Plume", "MAR Rainbow 354-166",
                                                "Cayman Shallow Plume 1", "Lau Basin Kilo Moana Plume 3",
                                                "Cayman Shallow Plume 3", "Cayman Shallow Plume 2",
                                                "Lau Basin Tahi Moana Plume 1", "Lau Basin Abe Plume 1",
                                                "Brothers NWCB S139",
                                                "Cayman Deep Plume 2", "Lau Basin Kilo Moana Plume 2",
                                                "Brothers NWCB S012",
                                                "Cayman Deep Plume 3", "Lau Basin Mariner Plume 1",
                                                "Cayman Deep Plume 1", "Brothers NWCB S140",
                                                "Lau Basin Kilo Moana Plume 1",
                                                "Brothers NWCB S146", 
                                                "ELSC Tui Malila 134-614", "Brothers NWCB S141",
                                                "Lau Basin Abe Plume 3", "EPR 4281-140",
                                                "ELSC Abe A3", "MAR Lucky 356-308",
                                                "Lau Basin Mariner Plume 2", "Axial Seawater",
                                                "Guaymas Basin Plume", "Lau Basin Tahi Moana Plume 2",
                                                "Lau Basin Kilo Moana Plume 4", 
                                                "ELSC Abe 128-326",
                                                "EPR PIR-30", "ELSC Vai Lili V2",
                                                "Lau Basin Tui Malila Plume", "Lau Basin Abe Plume 2", 
                                                #end Gamma
                                                "Brothers Diffuse S015", "Brothers Diffuse S009",
                                                "MAR Lucky 356-284", "Brothers NWCA S143",
                                                "ELSC Abe A1", 
                                                "Brothers NWCA S144", "ELSC Tui Malila T2",
                                                "ELSC Tui Malila T11", "Brothers LC S014",
                                                "MAR Rainbow 355-202",
                                                "Guaymas 4559-240", "ELSC Tui Malila T10",
                                                "Brothers LC S016", "Guaymas 4571-419",
                                                "ELSC Mariner 131-447", "ELSC Mariner M17",
                                                "ELSC Bowl M2", "Brothers NWCA S013",
                                                "Brothers UC S147", "ELSC Tui Malila 132-544",
                                                "ELSC Mariner M10", "Brothers NWCA S017",
                                                "ELSC Bowl M1",
                                                "Guaymas 4561-380","Brothers UC S010",
                                                "Brothers NWCA S145",
                                                "Brothers UC S011","Brothers NWCA S142","Guaymas 4561-384")) # start Camp lowest to highest 



####################### Stacked bar plot of Gamma and Campylo MAG abundance ################################
## Shows functional redundancy of microbial MAGs

#factor
level_order <- c('Gammaproteobacteria', 'Campylobacterota') 

####The following produces Figure X, which was modified in Biorender
dev.off()
plot <- abund_MAGs_long_norm_gtdb_p %>%
  ggplot(aes(x = as.numeric(value), y = Site, fill = factor(Taxa, levels = c(level_order)))) + #y = reorder(Site, value, sum) | factor(checkv_quality, levels = level_order)
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(begin = .5,
                       end = 0) +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Microbial MAG Abundance", y = "Site",
       fill = "Microbial taxonomy") +
  guides(fill=guide_legend(override.aes = list(size=8))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        strip.text.x = element_text(size = 12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_x_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  #geom_hline(yintercept = 21.5) +
  # annotate(geom="text", x=6, y=30, label="Plume",
  #          color="black") +
  #ggtitle("Viral Abundance") + #Change for top X grabbed
  facet_wrap(.~Locat) +
  scale_y_discrete(limits=rev)
#coord_flip()
plot

ggsave("output/coverm_MAG_CampGamma_normAbun.png", plot,
       height = 12, width = 15,
       bg = "transparent")



