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


