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
########################## inputs ############################

#read input
abun <- read.delim(file = "input/PlumeVentVirus-vs-Reads-CoverM.tsv", header = TRUE)
#read in metadata to change names
abun_names <- read.delim2(file = "input/AssembliesToReads_Mapping.txt", header = FALSE)

########################## name cleaning ############################
#clean names
colnames(abun) = gsub("gz.Relative.Abundance....", "gz", colnames(abun))
#drop SWIR
abun <- abun %>% select(-c("X58P_trim_1.fastq.gz", "SWIR_B_trim_1.fastq.gz"))

#remove .fasta
#from whole df
#abun_names <- as.data.frame(apply(abun_names,2, str_remove_all, "_metaspades_scaffolds.min1000.fasta|.fasta"))
abun_names$V1 <- gsub("_metaspades_scaffolds.min1000.fasta|.fasta","", abun_names$V1) 
#change dashes to dots
abun_names$V2 <- gsub("-","\\.", abun_names$V2) 

########################## add Site name and melt long ############################

#transform from wide to long
abun_long <- melt(abun, id = c("Genome")) 
#rename variable column to Site
abun_long<-rename(abun_long, "Site" = "variable")
#drop unmapped
abun_long<-abun_long[!grepl("unmapped", abun_long$Genome),]

#change the read names to site names
abun_long <- abun_names %>%
  dplyr::select("V1", "V2") %>%
  right_join(abun_long, by = c("V2" = "Site"))
abun_long<-rename(abun_long,"Site" = "V1") %>%
  select(-"V2")

#drop 0s
abun_long<-filter(abun_long, value > 0)

########################## add host metadata ############################

abun_long_iphop <- iphop %>%
  dplyr::select(Virus, Host.genus) %>%
  right_join(abun_long, by = c("Virus" = "Genome"))   
#only keep phylum and class of the tax string
abun_long_iphop <- abun_long_iphop %>% separate(Host.genus, c("d", "p", "c", "o", "f", "g"), 
                                        sep= ";")
#select only phylum and class taxonomy
abun_long_iphop <- abun_long_iphop %>% select(c("Virus","p","c","Site","value"))

#for Proteobacteria keep class, everything else, keep phylum
abun_long_proteo <- abun_long_iphop %>% filter(grepl("p__Proteobacteria", p))
abun_long_proteo <- abun_long_proteo %>% select(-c("p"))
abun_long_proteo <- abun_long_proteo %>% rename("Taxa" = "c")
abun_long_proteo <- abun_long_proteo %>% select(c("Virus","Taxa","Site","value"))


abun_long_iphop <- abun_long_iphop %>% filter(!grepl("p__Proteobacteria", p))
abun_long_iphop <- abun_long_iphop %>% select(-c("c"))
abun_long_iphop <- abun_long_iphop %>% rename("Taxa" = "p")
#put the data frames back together
abun_long_iphop <- rbind(abun_long_iphop, abun_long_proteo)

############################ Sum abundance by predicted host and site #############################

abun_long_iphop_p <- abun_long_iphop %>%
  group_by(Site, Taxa) %>% 
  summarise(value=sum(as.numeric(value))) %>% #summing the group
  filter(grepl("c__Gammaproteobacteria|p__Campylobacterota", Taxa)) %>% #only grab Gamma and Campylo
  ungroup() 

############################ Set order of sites for plotting #############################

# abun_iphop_p_Camp <- abun_long_iphop_p %>%
#   filter(grepl("p__Campylobacterota", Taxa)) %>%
#   group_by(Taxa) %>% 
#   arrange(desc(value))
# 
# abun_iphop_p_Gam <- abun_long_iphop_p %>%
#   filter(!grepl("p__Campylobacterota", Taxa)) %>%
#   group_by(Taxa) %>% 
#   arrange(value)
#   
# #put the data frames back together
# abun_long_iphop_p <- rbind(abun_iphop_p_Camp, abun_iphop_p_Gam)
#     
# # lock in factor level order
# #abun_long_iphop_p$value <- factor(abun_long_iphop_p$value, levels = abun_long_iphop_p$value)
# abun_long_iphop_p$Site <- factor(abun_long_iphop_p$Site, levels = unique(abun_long_iphop_p$Site))

#set order of x axis 
abun_long_iphop_p$Site <- factor(abun_long_iphop_p$Site, 
                                 levels=c("Lau_Basin_Kilo_Moana_min1000_4", "Lau_Basin_Tahi_Moana_min1000_2", 
                                          "Brothers_NWCB_S139", "Brothers_NWCB_S012",
                                          "Lau_Basin_Abe_min1000_2", "Lau_Basin_Abe_min1000_3", 
                                          "Cayman_Shallow_min1000_2", "MAR_Rainbow_354-166",
                                          "Brothers_NWCB_S140", "MAR_Rainbow_355-202", "Cayman_Deep_min1000_3", 
                                          "Brothers_NWCB_S141", "EPR_4281-140", "Lau_Basin_Kilo_Moana_min1000_2", 
                                          "Cayman_Shallow_min1000_1", "Lau_Basin_Tahi_Moana_min1000_1",
                                          "Brothers_NWCB_S146", "Cayman_Deep_min1000_2", "Lau_Basin_Kilo_Moana_min1000_3",
                                          "ELSC_Tui_Malila_134-614", "ELSC_Tui_Malila_T2", "Lau_Basin_Mariner_min1000_2",
                                          "Lau_Basin_Mariner_min1000_1", "Lau_Basin_Kilo_Moana_min1000_1",
                                          "Lau_Basin_Abe_min1000_1", "Cayman_Shallow_min1000_3", "MAR_Lucky_356-308", "ELSC_Abe_128-326",
                                          "ELSC_Abe_A3", "ELSC_Vai_Lili_V2", "MAR_Lucky_356-284", "Lau_Basin_Tui_Malila_min1000",
                                          "Cayman_Deep_min1000_1", "Axial_Plume_scaffolds_min1000", "Guaymas_Basin_min1000",
                                          "Brothers_Diffuse_S009", "Axial_Seawater_scaffolds_min1000", "EPR_PIR-30", #end Gamma highest to lowest
                                          "Brothers_Diffuse_S015", "ELSC_Abe_A1", "Brothers_NWCA_S143", "Brothers_LC_S014",
                                          "ELSC_Tui_Malila_T11", "ELSC_Tui_Malila_T10", "Brothers_NWCA_S144",
                                          "Guaymas_4571-419", "Guaymas_4559-240", "Brothers_LC_S016",
                                          "Guaymas_4561-380", "ELSC_Mariner_131-447", "Brothers_NWCA_S013",
                                          "ELSC_Mariner_M17", "ELSC_Bowl_M2", "ELSC_Tui_Malila_132-544",
                                          "Brothers_UC_S147", "ELSC_Bowl_M1", "Brothers_UC_S010", "Brothers_NWCA_S017",
                                          "Brothers_UC_S011", "ELSC_Mariner_M10", "Guaymas_4561-384",
                                          "Brothers_NWCA_S145", "Brothers_NWCA_S142")) # start Camp lowest to highest 
############################ Plot  #############################


####The following produces Figure X, which was modified in Biorender
dev.off()
plot <- abun_long_iphop_p %>%
  ggplot(aes(x = as.numeric(value), y = Site, fill = Taxa)) + #y = reorder(Site, value, sum)
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(begin = .5,
                       end = 0) +
  labs(x = "% Relative Abundance", y = "Site") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
        #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_x_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  ggtitle("Viral Abundance") #+ #Change for top X grabbed
  #facet_wrap(.~Taxa) 
#scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
plot

# ggsave("Output/coverm_abun.png", plot, width = 10, height = 5, dpi = 500,
#        bg = "transparent")
# ggsave("Output/coverm_abun.pdf", plot, width = 10, height = 5, dpi = 500,
#        bg = "transparent")







