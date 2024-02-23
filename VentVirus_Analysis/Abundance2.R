########################## Compiling and parsing virus CoverM abundance ############################
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(viridis)
library(pals)
library(RColorBrewer)

######################################### inputs ##############################################

#coverm abundance with read counts and genome coverage
abund <- read.delim2("input/PlumeVentVirus-vs-Reads-CoverM-Count_MinCov.tsv")

#seqkit number of reads per sample
read_length <- read.delim2(file = "../../abundance/stats.tsv")
read_length$file <- gsub("-",".",read_length$file) #replace symbol

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

#now group by genome and variable and keep rows where covered fraction is >=.7
abund_long_cov <- abund_long %>%
  group_by(Genome, variable) %>%
  filter(any(Type == "Covered.Fraction" & value >= 0.7) | !any(Type == "Covered.Fraction")) %>%
  ungroup() %>%
  group_by(variable, Type) %>%
  mutate(Percentile_Rank=rank(as.numeric(value))/length(as.numeric(value))) %>% #add col with percentile rank per category
  ungroup()

#map number of reads
abund_long_cov <- read_length %>%
  dplyr::select(c("file", "num_seqs")) %>%
  right_join(abund_long_cov, by = c("file" = "variable")) %>%
  rename("Reads" = "file")
abund_long_cov$value <- as.numeric(abund_long_cov$value)

#divide read count by number of reads *100
abund_long_norm <- abund_long_cov
abund_long_norm <- abund_long_norm %>%
  mutate(abun_norm = (value/num_seqs*100))

#WHEN YOU COME BACK DROP EVERYTHING BUT READ COUNT ON ABUND LONG NORM

##################################### host metadata ##############################################

abund_long_norm_iphop <- iphop %>%
  dplyr::select(Virus, "Host genus") %>%
  right_join(abund_long_norm, by = c("Virus" = "Genome"))   
#only keep phylum and class of the tax string
abund_long_norm_iphop <- abund_long_norm_iphop %>% separate("Host genus", c("d", "p", "c", "o", "f", "g"), 
                                                            sep= ";")
#select only phylum and class taxonomy
abund_long_norm_iphop <- abund_long_norm_iphop %>% select(c("Virus","p", "c", "o", "f", "g","V1","abun_norm"))
abund_long_norm_iphop <- abund_long_norm_iphop %>% 
  filter(str_detect(c, "c__Gammaproteobacteria") | str_detect(p, "p__Campylobacterota"))

#abund_long_norm_iphop <- abund_long_norm_iphop %>% 
#  filter(str_detect(c, "c__Gammaproteobacteria"))

abund_long_norm_iphop$Site <- gsub("min1000","Plume", abund_long_norm_iphop$V1) 
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

#abun_long_iphop <- abun_long_iphop %>% 
#  filter(str_detect(Locat, "Plume")) #%>%
# filter(str_detect(o, "o__PS1"))

# #for Proteobacteria keep class, everything else, keep phylum
abund_long_proteo <- abund_long_norm_iphop %>% filter(grepl("p__Proteobacteria", p))
abund_long_proteo <- abund_long_proteo %>% select(-c("p"))
abund_long_proteo <- abund_long_proteo %>% rename("Taxa" = "c")
abund_long_proteo <- abund_long_proteo %>% select(c("Virus","Taxa","Site","abun_norm", "Locat"))

abund_long_norm_iphop <- abund_long_norm_iphop %>% filter(!grepl("p__Proteobacteria", p))
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
                                       levels=c("Lau Basin Kilo Moana Plume 3","Lau Basin Abe Plume 1",
                                                "Lau Basin Tahi Moana Plume 1","Lau Basin Kilo Moana Plume 2",
                                                "Brothers NWCB S139","Brothers NWCB S012",
                                                "ELSC Mariner M17", "Axial Plume",
                                                "MAR Rainbow 355-202","Brothers NWCB S140",
                                                "Cayman Shallow Plume 1","MAR Rainbow 354-166",
                                                "Cayman Shallow Plume 2","Lau Basin Abe Plume 2",
                                                "Brothers NWCB S141","Cayman Deep Plume 2","Brothers NWCA S144",
                                                "Lau Basin Abe Plume 3","ELSC Tui Malila 134-614",
                                                "Cayman Shallow Plume 3","EPR 4281-140","ELSC Tui Malila T2",
                                                "Brothers NWCB S146","Cayman Deep Plume 3",
                                                "Brothers NWCA S143","MAR Lucky 356-308",
                                                "Cayman Deep Plume 1","MAR Lucky 356-284","ELSC Abe 128-326",
                                                "Lau Basin Mariner Plume 2","ELSC Abe A3",
                                                "Lau Basin Kilo Moana Plume 1","Lau Basin Mariner Plume 1",
                                                "ELSC Vai Lili V2","Guaymas Basin Plume","Axial Seawater",
                                                "EPR PIR-30","Brothers Diffuse S009","Lau Basin Kilo Moana Plume 4",
                                                "Brothers Diffuse S015","Lau Basin Tui Malila Plume",
                                                "Lau Basin Tahi Moana Plume 2", #end Gamma
                                                "ELSC Abe A1","ELSC Tui Malila T11",
                                                "ELSC Tui Malila T10","Brothers LC S014","Brothers LC S016",
                                                "Guaymas 4559-240","ELSC Bowl M2",
                                                "ELSC Tui Malila 132-544","Brothers NWCA S013",
                                                "ELSC Mariner 131-447","Guaymas 4571-419","ELSC Mariner M10",
                                                "Guaymas 4561-380","Brothers UC S010","Brothers NWCA S017",
                                                "Brothers UC S147","Brothers NWCA S145","ELSC Bowl M1",
                                                "Brothers UC S011","Guaymas 4561-384","Brothers NWCA S142")) # start Camp lowest to highest 


############################ Plot  #############################

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
       height = 10, width = 13,
       bg = "transparent")
