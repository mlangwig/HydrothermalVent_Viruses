##################### Visualizing viral gene content/AMGs #####################

##################### Read in/Create major inputs #####################

sulfur_viruses <- read.csv("../iPHoP/output/iphop_VentPlume_sulfur_50comp.csv", header = TRUE)

sulfur_vMAG_master <- master_table_vMAGs[master_table_vMAGs$vMAG %in% sulfur_viruses$Virus,] #subset table using list of names
  
sulfur_unbinned_master <- master_table_unbinned[master_table_unbinned$scaffold %in% sulfur_viruses$Virus,]

non_AMGs <- read.delim("input/non_AMG_list.txt", header = FALSE)

##################### Filter for AMGs #####################

##################### Filter for AMGs bordered by VOGs #####################

sulfur_unbinned_master_temp <- sulfur_unbinned_master %>%
  select(c("protein", "KO", "AMG", "VOG"))

sulfur_vMAG_master_temp <- sulfur_vMAG_master %>%
  select(c("protein", "KO", "AMG", "VOG"))

sulfur_VentVirus_AMGs_filtered <- rbind(sulfur_unbinned_master_temp, sulfur_vMAG_master_temp)

library(DescTools)

sulfur_VentVirus_AMGs_filtered <- sulfur_VentVirus_AMGs_filtered[OrderMixed(sulfur_VentVirus_AMGs_filtered$protein), ]
sulfur_VentVirus_AMGs_filtered$protein2 <- sub("_[^_]+$", "", sulfur_VentVirus_AMGs_filtered$protein)

test <- sulfur_VentVirus_AMGs_filtered %>%
  mutate_all(na_if,"")

## Loop for detecting AMG and Vogs
protein.vec <- unique(test$protein2)
amg.vog.filter <- data.frame(protein = NA,
                             KO = NA,
                             AMG = NA,
                             VOG = NA,
                             protein2 = NA)

for(i in 1:length(protein.vec)){
  
  ## Grab the protein
  protein.tmp <- protein.vec[i]
  
  ## Subset by the protein
  test.tmp <- test[test$protein2 == protein.tmp, ]
  
  ## Look for any AMG
  amg.detect <- any(!is.na(test.tmp$AMG))
  
  ## If this is true then proceed
  if(amg.detect == TRUE){
  
    ## Where are the AMGs
    amg.pos <- which(!is.na(test.tmp$AMG))
    
    ## Where are the VOGS
    vog.pos <- which(!is.na(test.tmp$VOG))
    
    ## Where are they the same
    amg.vog.match <- which(!is.na(test.tmp$AMG) & !is.na(test.tmp$VOG))
    
    ## See if we have at least 1 match
    ## If so we subset the DF for these row since we definitely want them
    if(length(amg.vog.match) >= 1){
      match.df <- test.tmp[amg.vog.match, ]
      amg.vog.filter <- rbind(amg.vog.filter, match.df)
    }
    
    ## See where we have an amg that doesn't match the amg.vog.match
    amg.pos <- amg.pos[!amg.pos %in% amg.vog.match]
    
    ## If we have some non-matches
    if(length(amg.pos) >= 1){
    for(j in 1:length(amg.pos)){
      amg.pos.tmp <- amg.pos[j]
      ## Test if the vog position has an amg before or after
      if(amg.pos.tmp < max(vog.pos) & amg.pos.tmp > min(vog.pos)){
        ## if so we keep these rows
        row.keep <- c(amg.vog.match, amg.pos.tmp)
        amg.sandwich <- test.tmp[row.keep, ]
        amg.vog.filter <- rbind(amg.vog.filter, amg.sandwich)
      }
    }
  
    }
  } else {
    print(paste0("No AMGs detected for protein iteration:", i))
  }
  }

amg.vog.filter <- amg.vog.filter[-1, ]

####### remove non-AMGs that I know based on mmseqs protein #######
############ clustering of all viruses - things that are present 
############ at all geographically distinct Vent sites
`%notin%` <- Negate(`%in%`)

amg.vog.filter <- amg.vog.filter[amg.vog.filter$KO %notin% non_AMGs$V1,] 
#sulfur_vMAG_AMGs <- sulfur_vMAG_AMGs[sulfur_vMAG_AMGs$Pfam %notin% non_AMGs$V1,] 
amg.vog.filter <- amg.vog.filter[amg.vog.filter$VOG %notin% non_AMGs$V1,] 

#### map the filtered list back to meaningful virus genome names

sulfur_vMAG_AMGs <- sulfur_vMAG_master[sulfur_vMAG_master$KO %in% amg.vog.filter$KO,] 
sulfur_unbinned_AMGs <- sulfur_unbinned_master[sulfur_unbinned_master$KO %in% amg.vog.filter$KO,]

# sulfur_vMAG_AMGs <- sulfur_vMAG_master %>% filter(grepl("AMG", AMG))
# sulfur_unbinned_AMGs <- sulfur_unbinned_master %>% filter(grepl("AMG", AMG))

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
# sulfur_VentVirus_AMGs <- sulfur_viruses %>%
#   dplyr::select("Virus", "c") %>%
#   right_join(sulfur_VentVirus_AMGs, by = c("Virus" = "Virus")) %>%
#   unique()

sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>%
  mutate(c = str_replace(c, "c__", ""))
sulfur_VentVirus_AMGs$c <- sulfur_VentVirus_AMGs$c %>%
  na_if("")
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
#REDO THE TOP COMMANDS TO GET THE RIGHT FILE FORMAT BACK - up to make input for plotting
sulfur_virus_KO_list1 <- sulfur_unbinned_AMGs %>%
  distinct(KO)
sulfur_virus_KO_list2 <- sulfur_vMAG_AMGs %>%
  distinct(KO)

sulfur_virus_KO_list <- rbind(sulfur_virus_KO_list1, sulfur_virus_KO_list2) %>%
  distinct(KO)

# write.table(sulfur_virus_KO_list, file = "output/sulfur_virus_KO_list.txt", 
#             col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


############ map the KO pathway to the genes ##############
KO_pathway <- read.delim(file = "input/sulfur_virus_KO_list_mod.txt", header = TRUE, sep = "\t")

KO_pathway <- KO_pathway %>% 
  separate(Description, c("Description", NA), sep= ",")
KO_pathway$Description <- gsub("K23144","UAP", KO_pathway$Description)
KO_pathway$Description <- gsub("K13522", "nadM", KO_pathway$Description)
KO_pathway$Description <- gsub("E2.1.1.104", "CCoAOMT", KO_pathway$Description)
KO_pathway$Description <- gsub("rhnA-cobC", "rhnA", KO_pathway$Description)
KO_pathway$Description <- gsub("E2.7.7.3A", "coaD", KO_pathway$Description)
KO_pathway$Description <- gsub("K16150", "GS", KO_pathway$Description)
KO_pathway$Description <- gsub("E2.7.7.24", "rfbA", KO_pathway$Description)
KO_pathway$Description <- gsub("E4.2.1.46", "rfbB", KO_pathway$Description)

#map to plotting file
sulfur_VentVirus_AMGs <- KO_pathway %>%
  dplyr::select("Description", "Pathway") %>%
  right_join(sulfur_VentVirus_AMGs, by = c("Description" = "KO")) %>%
  unique()

# #removing genes I know are not good context by manual examination
# sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>%
#   filter(!grepl("CS|fadJ|egtC|ATPF0A", Description))
#   
# sulfur_VentVirus_AMGs <- unique(sulfur_VentVirus_AMGs)

### add column for site of virus to try and make plot smaller
sulfur_VentVirus_AMGs$VirusSite <- sulfur_VentVirus_AMGs$Virus # copy column
sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>% 
  separate(VirusSite, c("VirusSite", NA), sep= "(?=_NODE|_k95|_scaffold|_vRhyme)") %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*ELSC.*","Lau_Basin")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*Brothers.*","Brothers_Volcano")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*Guaymas.*","Guaymas_Basin")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*MAR.*","Mid_Atlantic_Ridge")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*EPR.*","East_Pacific_Rise")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*Cayman.*","Mid_Cayman_Rise")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*Lau.*","Lau_Basin")) %>%
  mutate_at(vars(VirusSite), ~ str_replace(., ".*Axial.*","Axial_Seamount")) 

#to count by Site rather than by Virus
sulfur_VentVirus_AMGs <- sulfur_VentVirus_AMGs %>%
  group_by(VirusSite, c, Description, Pathway) %>%
  summarise("KO_count_perSite" = sum(KO_count))

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
p <- ggplot(sulfur_VentVirus_AMGs, aes(y=Description, x=VirusSite))+
  #geom_point(aes(size=KO_count, colour = c))+ 
  geom_point(color='black', shape = 21, stroke = .15, aes(fill=factor(c), size=KO_count_perSite)) + # THE SHAPE = 21 IS CRITICAL TO GET THE CIRCLE BORDER
  scale_size_continuous(breaks = c(1, 6, 12), range = c(4, 14))+
  scale_fill_discrete(guide="legend", type = col_vector)+ #type = "viridis for color
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=360, size = 22),
        strip.text.x = element_text(angle = 90, size = 22),
        #panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, size = 19, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 20),
        text = element_text(color="black"),
        legend.position="right",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        panel.spacing = unit(.15, "lines"))+
  scale_x_discrete(position="bottom")+
  guides(fill = guide_legend(override.aes = list(size = 10)))+
  xlab("Viruses by Site")+
  ylab("Gene")+
  labs(size = "Gene count",
       fill = "Host class") +
  facet_grid(Pathway ~ c, scales = "free", space = "free") + #wow this took me awhile - remember you can't factor the y axis and then get it to facet properly
  scale_y_discrete(limits=rev) #+ # has to be this instead of factor 
  # coord_flip()
p  

#ggsave("output/VentVirus_AVGs_wide_QC.png", p, width = 30, height = 14)

ggsave("output/VentVirus_AVGs_long_QC.png", p, width = 22, height = 25)
ggsave("output/VentVirus_AVGs_long_QC.svg", p, width = 20, height = 25)

################### Spencer helping with