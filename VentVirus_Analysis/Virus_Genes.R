##################### Visualizing viral gene content/AMGs #####################

##################### Read in/Create major inputs #####################

sulfur_viruses <- read.csv("../iPHoP/output/iphop_VentPlume_sulfur_50comp.csv", header = TRUE)

sulfur_vMAG_master <- master_table_vMAGs[master_table_vMAGs$vMAG %in% sulfur_viruses$Virus,] #subset table using list of names
  
sulfur_unbinned_master <- master_table_unbinned[master_table_unbinned$scaffold %in% sulfur_viruses$Virus,]

non_AMGs <- read.delim("input/non_AMG_list.txt", header = FALSE)
##################### Filter for AMGs #####################
sulfur_vMAG_AMGs <- sulfur_vMAG_master %>% filter(grepl("AMG", AMG))
sulfur_unbinned_AMGs <- sulfur_unbinned_master %>% filter(grepl("AMG", AMG))

#remove non-AMGs that I know based on mmseqs protein clustering of all viruses - things that are present at all geographically distinct Vent sites
`%notin%` <- Negate(`%in%`)

sulfur_vMAG_AMGs <- sulfur_vMAG_AMGs[sulfur_vMAG_AMGs$KO %notin% non_AMGs$V1,] #subset table using list of names
sulfur_vMAG_AMGs <- sulfur_vMAG_AMGs[sulfur_vMAG_AMGs$Pfam %notin% non_AMGs$V1,] #subset table using list of names
sulfur_vMAG_AMGs <- sulfur_vMAG_AMGs[sulfur_vMAG_AMGs$VOG %notin% non_AMGs$V1,] #subset table using list of names

sulfur_unbinned_AMGs <- sulfur_unbinned_AMGs[sulfur_unbinned_AMGs$KO %notin% non_AMGs$V1,] #subset table using list of names
sulfur_unbinned_AMGs <- sulfur_unbinned_AMGs[sulfur_unbinned_AMGs$Pfam %notin% non_AMGs$V1,] #subset table using list of names
sulfur_unbinned_AMGs <- sulfur_unbinned_AMGs[sulfur_unbinned_AMGs$VOG %notin% non_AMGs$V1,]
##################### Filter for those bordered by VOGs #####################

##################### Make input for plotting #####################
sulfur_vMAG_AMGs <- sulfur_vMAG_AMGs %>% 
  select(c("vMAG", "KO.name")) %>%
  separate(KO.name, c("KO", NA), sep= ";") %>%
  rename(vMAG,"Virus" = "vMAG")

sulfur_unbinned_AMGs <- sulfur_unbinned_AMGs %>% 
  select(c("scaffold", "KO.name")) %>%
  separate(KO.name, c("KO", NA), sep= ";") %>%
  rename(scaffold,"Virus" = "scaffold")

sulfur_VentVirus_AMGs <- rbind(sulfur_vMAG_AMGs, sulfur_unbinned_AMGs)

#transform for plotting
sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>%
  group_by(Virus) %>%
  count(KO, name = "KO_count")
#add metadata
sulfur_VentVirus_AMGs <- sulfur_viruses %>%
  dplyr::select("Virus", "c") %>%
  right_join(sulfur_VentVirus_AMGs, by = c("Virus" = "Virus")) 

sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>%
  mutate(c = str_replace(c, "c__", "")) %>%
  na_if('')
sulfur_VentVirus_AMGs$c <- sulfur_VentVirus_AMGs$c %>% replace_na("Bacteria") #replace blank cell with Bac

#fix gene names
sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>% 
  separate(KO, c("KO", NA), sep= ",")
sulfur_VentVirus_AMGs$KO <- gsub("K23144","UAP", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("K13522", "nadM", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("E2.1.1.104", "CCoAOMT", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("rhnA-cobC", "rhnA", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("E2.7.7.3A", "coaD", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("K16150", "GS", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("E2.7.7.24", "rfbA", sulfur_VentVirus_AMGs$KO)
sulfur_VentVirus_AMGs$KO <- gsub("E4.2.1.46", "rfbB", sulfur_VentVirus_AMGs$KO)

#subset unique list of gene names so can map metabolic pathway onto it
sulfur_virus_KO_list <- sulfur_VentVirus_AMGs %>%
  distinct("KO") %>%
  pull()


  rename("unique(sulfur_VentVirus_AMGs$KO)","KO" = "unique(sulfur_VentVirus_AMGs$KO)")

############# plot ##################
# library(randomcoloR)
# n <- 13
# palette <- (distinctColorPalette(n))
# 
# my_colors <- (pals::kelly(n=13)[2:14])
col_vector<-c("#7FC97F", "#d9d9d9", "#FDC086",
              "#FFFF99", "#386CB0", "#F0027F",
              "#a6cee3", "#ff7f00", "#1B9E77",
              "#A6761D", "#7570B3", "#e9a3c9",
              "#000000") #"#E6AB02", "#A6761D", "#666666")

dev.off()
p <- ggplot(sulfur_VentVirus_AMGs, aes(y=Virus, x=KO))+
  #geom_point(aes(size=KO_count, colour = c))+ 
  geom_point(color='black', shape = 21, stroke = .15, aes(fill=factor(c), size=KO_count)) + # THE SHAPE = 21 IS CRITICAL TO GET THE CIRCLE BORDER
  scale_size_continuous(breaks = c(1, 2, 3), range = c(1, 4))+
  scale_fill_discrete(guide="legend", type = col_vector)+ #type = "viridis for color
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=360, size = 7),
        #panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-.10, size = 5),
        axis.text.y = element_text(size = 5),
        text = element_text(color="black"),
        legend.position="right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.title.y = element_text(size = 8))+
  scale_x_discrete(position="top")+
  xlab("")+
  ylab("Virus")+
  labs(size = "Gene count",
       fill = "Host class") +
  facet_grid(c ~ ., scales = "free_y", space = "free") + #wow this took me awhile - remember you can't factor the y axis and then get it to facet properly
  scale_y_discrete(limits=rev) # has to be this instead of factor
p  

#ggsave("output/iphop_sulfur_50comp_genes.png", width = 7, height = 6)

