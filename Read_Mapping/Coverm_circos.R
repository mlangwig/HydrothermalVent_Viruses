#install.packages("circlize")
library(circlize)
library(stringi)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(viridis)
library(pals)
library(RColorBrewer)

########################### major inputs ################################################

#coverm count data
coverm_rc <- read.delim2(file = "../../abundance/PlumeVentVirus-vs-Reads-CoverM-Count_MinCov.tsv")

#metadata to change names
abun_names <- read.delim2(file = "../VentVirus_Analysis/input/AssembliesToReads_Mapping.txt", header = FALSE)

########################### convert wide to long coverm ################################################

#drop SWIR
coverm_rc <- coverm_rc %>% select(-contains("X58P_trim_1.fastq.gz"))
coverm_rc <- coverm_rc %>% select(-contains("SWIR_B_trim_1.fastq.gz"))
#remove unmapped
coverm_rc <- coverm_rc[-1,]

#remove contam viruses
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

#making file for Spencer
coverm_rc <- coverm_rc %>% select(-contains("Relative.Abundance...."))
colnames(coverm_rc) = gsub(".Read.Count", "", colnames(coverm_rc))
coverm_rc <- coverm_rc %>%
  filter(!str_detect(Genome, remove.list))
coverm_rc$Genome_Site <- coverm_rc$Genome
coverm_rc <- coverm_rc %>% 
  separate(Genome_Site, c("Genome_Site", NA), sep = "_NODE|_scaffold|_vRhyme|_k95")

#faster replace all the naming patterns
#test <- coverm_rc_long
coverm_rc$Genome_Site <- stri_replace_all_regex(coverm_rc$Genome_Site,
                                                     pattern=c("_A[0-9]",
                                                               "_T[0-9][0-9]", "_T[0-9]", "_S0[0-9][0-9]",
                                                               "_S1[0-9][0-9]", "_[0-9][0-9][0-9]-[0-9][0-9][0-9]",
                                                               "-38[0-9]"),
                                                     replacement='',
                                                     vectorize=FALSE)

coverm_rc$Genome_Site <- gsub("*_M1[0-9]","",coverm_rc$Genome_Site)
coverm_rc$Genome_Site <- gsub("*_M[0-9]","",coverm_rc$Genome_Site)

#change order to bring Genome_Site to front
coverm_rc <- coverm_rc %>%
  select(Genome_Site, everything())

#convert to numeric
# Identify columns ending with "Covered.Fraction"
fraction_columns <- grep("Covered\\.Fraction$", names(coverm_rc), value = TRUE)
#convert those columns
coverm_rc <- coverm_rc %>%
  mutate_at(vars(all_of(fraction_columns)), as.numeric)

#write_delim(coverm_rc, file = "Output/coverm_rc_output_CovFrac.tsv", delim = "\t")


######################## change read names to site names for intra vent comparison ################################################

#clean up names file for mapping
abun_names$V1 <- gsub("_metaspades_scaffolds.min1000.fasta|.fasta|_min1000_[0-9]|_min1000|_scaffolds","", abun_names$V1) 
abun_names$V1 <- gsub("_metaspades_scaffolds.min1000.fasta|.fasta|_min1000_[0-9]|_min1000|_scaffolds","", abun_names$V1) 
#change dashes to dots - CRUCIAL OR YOU DROP BROTHERS NAMES
abun_names$V2 <- gsub("-","\\.", abun_names$V2)

#make shorter
abun_names$V1 <- gsub("_A[0-9]","",abun_names$V1)
abun_names$V1 <- gsub("_M10","",abun_names$V1)
abun_names$V1 <- gsub("_T[0-9][0-9]","",abun_names$V1)
abun_names$V1 <- gsub("_T[0-9]","",abun_names$V1)
abun_names$V1 <- gsub("_M10","",abun_names$V1)
abun_names$V1 <- gsub("*_S0[0-9][0-9]","",abun_names$V1)
abun_names$V1 <- gsub("*_S1[0-9][0-9]","",abun_names$V1)
abun_names$V1 <- gsub("*_M1[0-9]","",abun_names$V1)
abun_names$V1 <- gsub("*_M[0-9]","",abun_names$V1)
abun_names$V1 <- gsub("*_[0-9][0-9][0-9]-[0-9][0-9][0-9]","",abun_names$V1)
abun_names$V1 <- gsub("*-380","",abun_names$V1)
abun_names$V1 <- gsub("*-384","",abun_names$V1)

#write_delim(abun_names, file = "Output/read_name_mapping.tsv", delim = "\t")

######################## Make symmetrical matrix input ################################################
########### the following code was written by Spencer R. Keyser (skeyser@wisc.edu)

coverm <- coverm_rc
read.map <- abun_names
colnames(read.map) <- c("Site", "UID")

## filter any column with the string "Covered.Fraction" >= 0.7
## filter any column with the read count greater than 5
for(i in 1:length(colnames(coverm))){
  name.tmp <- colnames(coverm)[i]
  if(str_detect(name.tmp, "Genome|Genome_Site|Covered.Fraction")){ next } else {
    colnames(coverm)[i] <- paste0(name.tmp, ".ReadCount")
  }
}
colnames(coverm)

## Make the dataframe long format and filter for >= 5 reads, >=70% covered fraction
coverm.long <- coverm %>%
  pivot_longer(
    cols = -1:-2
  ) %>%
  mutate(Id.var = case_when(str_detect(name, ".Covered.Fraction") ~ "CF",
                            str_detect(name, ".ReadCount") ~ "RC")) %>%
  mutate(RC_Ind = if_else(Id.var == "RC" & value >= 5, T, F)) %>%
  mutate(CF_Ind = if_else(Id.var == "CF" & value >= 0.7, T, F))

#add genome size info to coverm long to filter out <3kb


## Pull the read count part out
coverm.reads <- coverm.long %>%
  filter(Id.var == "RC") %>%
  select(Genome_Site, Genome, name, Reads = value) %>%
  mutate(name = gsub(".ReadCount", "", name))

## Pull the cover fraction part out
coverm.cf <- coverm.long %>%
  filter(Id.var == "CF") %>%
  select(Genome_Site, Genome, name, CoverFrac = value) %>%
  mutate(name = gsub(".Covered.Fraction", "", name))

## Join these two back
cover.full <- left_join(coverm.cf, coverm.reads)

## filter by the criteria
coverm.filter <- cover.full %>%
  filter(CoverFrac >= 0.7 & Reads >= 5)

##Read Mapping AND REMOVES WHEN GEN SITE=READ SITE
coverm.map <- coverm.filter %>%
  left_join(read.map, by = c("name" = "UID")) %>%
  rename(Read_Site = Site) %>%
  filter(Genome_Site != Read_Site) %>%
  select(c("Genome_Site", "Read_Site", "Genome", "name", "CoverFrac", "Reads"))

## Add the genome sizes and remove >3kb
coverm.map <- gensize %>%
  dplyr::select(file, sum_len) %>%
  right_join(coverm.map, by = c("file" = "Genome"))
coverm.map <- rename(coverm.map, "Genome" = "file")

#filter out <3kb
coverm.map <- coverm.map %>%
  filter(sum_len >= 3000)

## Reduce to just the genome_site and read_site
# REDO THIS PART HERE AND SKIP GD PART BELOW TO HAVE ALL IN ONE PLOT 
coverm.simple <- coverm.map %>%
  select(Genome_Site, Read_Site) %>%
  group_by(Genome_Site, Read_Site) %>%
  count()


## Make the Genome_Site x Read_Site Matrix with # of connections
######################### the following code was written by Spencer R. Keyser (skeyser@wisc.edu)
## make a zero-filled matrix
t.mat <- matrix(nrow = length(unique(coverm.simple$Genome_Site)), ncol = length(unique(coverm.simple$Read_Site)), data = 0)
colnames(t.mat) <- unique(coverm.simple$Genome_Site)
rownames(t.mat) <- unique(coverm.simple$Genome_Site)

## Holder dataframe
for(i in 1:nrow(t.mat)){
  g.site.tmp <- rownames(t.mat)[i]
  for(j in 1:ncol(t.mat)){
    r.site.tmp <- colnames(t.mat)[j]
    val.tmp <- coverm.simple %>%
      filter(Genome_Site == g.site.tmp & Read_Site == r.site.tmp) %>%
      pull(n)
    if(is.integer(val.tmp) && length(val.tmp) == 0L){
      print("No connections")
    } else {
      t.mat[i,j] <- val.tmp
    }
  }
}

## Make the matrix symmetric by adding the transverse of the asymmetric mat
t.sym <- t.mat + t(t.mat)

######### create matrix for color mapping border ##########
#The following loop was written by ChatGPT:
result_matrix <- matrix("black", nrow = nrow(t.sym), ncol = ncol(t.sym), dimnames = dimnames(t.sym))

# Iterate through the matrix and update the new matrix
for (i in seq_len(nrow(t.sym))) {
  for (j in seq_len(ncol(t.sym))) {
    # Extract the first part of the row and column names before "_"
    row_prefix <- sub("_.*", "", rownames(t.sym)[i])
    col_prefix <- sub("_.*", "", colnames(t.sym)[j])
    
    # Check if the prefixes are not identical and the value is greater than 0
    if (row_prefix != col_prefix && t.sym[i, j] > 0) {
      result_matrix[i, j] <- "red"
    }
  }
}

# Print the result
print(result_matrix)

######################## plot it in circlize ################################################

nm = unique(unlist(dimnames(t.sym)))
#group = structure(gsub("\\d", "", nm), names = nm)
#df for messing with names
group <- as.data.frame(nm, row.names = NULL) %>%
  rename("vent" = "nm")
group$group = group$vent
#adjust group names
group$group <- gsub(".*Lau_Basin.*","Lau Basin Plume", group$group)
group$group <- gsub(".*ELSC.*","Lau Basin Deposit", group$group)
group$group <- gsub(".*Brothers.*","Brothers Volcano", group$group)
group$group <- gsub("Guaymas_[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9]","Guaymas Basin Deposit", group$group)
group$group <- gsub("Guaymas_Basin","Guaymas Basin Plume", group$group)
group$group <- gsub("Guaymas_[0-9][0-9][0-9][0-9]","Guaymas Basin Deposit", group$group)
group$group <- gsub(".*Cayman.*","Mid-Cayman Rise", group$group)
group$group <- gsub(".*Axial.*","Axial Seamount", group$group)
group$group <- gsub(".*EPR.*","East Pacific Rise", group$group)
group$group <- gsub(".*MAR.*","Mid-Atlantic Ridge", group$group)
#adjust vent names
#group <- tibble::rownames_to_column(group, "vent")
group$vent <- gsub("Lau_Basin_","P_", group$vent)
group$vent <- gsub("ELSC_","D_", group$vent)
group$vent <- gsub("Brothers_","", group$vent)
group$vent <- gsub("Guaymas_Basin","P_Guaymas", group$vent)
group$vent <- gsub("Guaymas_","", group$vent)
group$vent <- gsub("Cayman_","", group$vent)
group$vent <- gsub("Axial_","", group$vent)
group$vent <- gsub("EPR_","", group$vent)
group$vent <- gsub("MAR_","", group$vent)
group$vent <- gsub("_"," ", group$vent)

m.connect_test <- t.sym
rownames(m.connect_test) <- paste(group$vent)
colnames(m.connect_test) <- paste(group$vent)

#also change col names on mapping matrix
rownames(result_matrix) <- paste(group$vent)
colnames(result_matrix) <- paste(group$vent)


# library(data.table)
# m <- data.table(as.vector(m.connect_test), as.vector(col(m.connect_test)), as.vector(row(m.connect_test)))
# m <- m[ order(m[,1]), ]

#vector again
#group <- dplyr::pull(group)
group2 <- setNames(as.character(group$group), group$vent)

#set colors
grid.col = c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478", "East Pacific Rise" = "#CE9A28",
             "Guaymas Basin Deposit" = "#28827A", "Guaymas Basin Plume" = "#499e97",
             "Lau Basin Plume" = "#3F78C1", "Lau Basin Deposit" = "#72a0db",
             "Mid-Atlantic Ridge" = "#8c510a", "Mid-Cayman Rise" = "#000000",
             #distinct vent locations above
             "4281-140" = "#CE9A28", "4559-240" = "#28827A",
             "4561" = "#28827A", "4571-419" = "#28827A", 
             #Guaymas deposit
             "D Mariner" = "#72a0db", "D Abe" = "#72a0db", "D Vai Lili V2" = "#72a0db", 
             "D Tui Malila" = "#72a0db", "D Bowl" = "#72a0db",
             #Lau deposit
             "Deep" = "#000000", "Shallow" = "#000000", 
             #MCR
             "Diffuse" = "#B56478", "NWCA" = "#B56478", "LC" = "#B56478", 
             "NWCB" = "#B56478", "UC" = "#B56478", 
             #Brothers Volcano
             "Lucky" = "#8c510a", "Rainbow" = "#8c510a",
             #MAR
             "P Abe" = "#3F78C1", "P Kilo Moana" = "#3F78C1", "P Mariner" = "#3F78C1",
             "P Tahi Moana" = "#3F78C1", "P Tui Malila" = "#3F78C1", 
             #Lau plume
             "Seawater" = "#4F508C", "Plume" = "#4F508C", 
             #Axial Seamount
             "PIR-30" = "#CE9A28", 
             #EPR
             "P Guaymas" = "#499e97" 
             #Guaymas plume
) #specify colors of sites/outer ring

col_fun = "grey" #specify color of links between them

#reset before running
dev.off()
circos.clear()
svg("Output/Circos_coverm_gd_iv_3kb_70mincov.svg")
#set font size
par(cex = 1.8, mar = c(0, 0, 0, 0))
#set gaps between blocks
#circos.par(gap.after = c(rep(.5, length(unique(colnames(t.mat_test)))))) #this is crucial to plot

#chordDiagram
chordDiagram(m.connect_test, 
             group = group2, 
             grid.col = grid.col,
             annotationTrack = c("grid"),
             annotationTrackHeight = c(.06, .06),
             #link.border = "black",
             link.border = result_matrix,
             link.sort = FALSE,
             transparency = 0.3,
             symmetric = TRUE,
             big.gap = 8,
             col = col_fun,
             preAllocateTracks = list(
               track.height = mm_h(6),
               track.margin = c(mm_h(1), 0)
             ))
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  # circos.text(CELL_META$xcenter,
  #             ylim[1] + cm_h(3), 
  #             labels = gsub("\\..*","", sector.name),
  #             col = "black", cex = 0.5, facing = "bending", 
  #             adj = c(0.5, 0.5), niceFacing = TRUE, font = 0.5)
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.5, niceFacing = TRUE,
              col = "white")
}, bg.border = NA)

highlight.sector(sector.index = c("P Abe", "P Kilo Moana",
                                  "P Mariner", "P Tahi Moana",
                                  "P Tui Malila"), track.index = 1, col = "#3F78C1", 
                 text = "Lau Basin Plume", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("D Mariner", "D Tui Malila",
                                  "D Bowl", "D Abe",
                                  "D Vai Lili V2"), track.index = 1, col = "#72a0db", 
                 text = "Lau Basin Deposit", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("NWCA", "NWCB", "LC",
                                  "UC", "Diffuse"), track.index = 1, col = "#B56478", 
                 text = "Brothers Volcano", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4561", "4571-419",
                                  "4559-240"), track.index = 1, col = "#28827A", 
                 text = "Guaymas Basin Deposit", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("P Guaymas"), track.index = 1, col = "#499e97", 
                 text = "Guaymas Basin Plume", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Deep", "Shallow"), track.index = 1, col = "#000000", 
                 text = "Mid-Cayman Rise", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Plume", "Seawater"), track.index = 1, col = "#4F508C", 
                 text = "Axial Seamount", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4281-140", "PIR-30"), track.index = 1, col = "#CE9A28", 
                 text = "East Pacific Rise", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Lucky", "Rainbow"), track.index = 1, col = "#8c510a", 
                 text = "Mid-Atlantic Ridge", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")

dev.off()




######################## make circos plot of geo distinct ##################


#create a version of just geo distinct to plot
coverm.simple_gd <- coverm.simple %>%
  mutate(gen_gd = Genome_Site) %>%
  mutate(read_gd = Read_Site) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*Lau_Basin.*", "Lau_Basin")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*Brothers.*", "Brothers")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*ELSC.*", "ELSC")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*Axial.*", "Axial")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*Cayman.*", "Cayman")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*EPR.*", "EPR")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*Guaymas_[0-9].*", "Guaymas_D")) %>%
  mutate_at(vars(gen_gd, read_gd), ~ str_replace(., ".*MAR.*", "MAR"))

#filter it
coverm.simple_gd <- coverm.simple_gd %>%
  filter(gen_gd != read_gd)

#drop the other columns to plot (broad names will be incorporated as group)
coverm.simple_gd <- coverm.simple_gd %>%
  select(c('Genome_Site', 'Read_Site', 'n')) %>%
  ungroup()
#add the missing site to Genome_Site column so name distribution is even
#this is a dummy line and won't contribute to count
coverm.simple_gd <- coverm.simple_gd %>%
  add_row(Genome_Site = 'Guaymas_4559-240',
          Read_Site = 'Lau_Basin_Abe',
          n = as.integer(0L))
  
# #only plot ones that are geo distinct like Fig 2A in manuscript
# #DON'T DO THIS IF YOU WANT TO PLOT ALL IN SAME FIGURE
# coverm.simple <- coverm.simple_gd\

#using coverm.simple_gd as input
coverm.simple_gd_sum <- coverm.simple_gd %>%
  ungroup() %>%
  select(c('gen_gd', 'read_gd', 'n')) %>%
  group_by(gen_gd, read_gd) %>%
  summarise(n = sum(n)) %>%
  ungroup()

## Make the Genome_Site x Read_Site Matrix with # of connections
######################### the following code was written by Spencer R. Keyser (skeyser@wisc.edu)
## make a zero-filled matrix
t.mat <- matrix(nrow = length(unique(coverm.simple_gd$Genome_Site)), ncol = length(unique(coverm.simple_gd$Read_Site)), data = 0)
colnames(t.mat) <- unique(coverm.simple_gd$Genome_Site)
rownames(t.mat) <- unique(coverm.simple_gd$Genome_Site)

## Holder dataframe
for(i in 1:nrow(t.mat)){
  g.site.tmp <- rownames(t.mat)[i]
  for(j in 1:ncol(t.mat)){
    r.site.tmp <- colnames(t.mat)[j]
    val.tmp <- coverm.simple_gd %>%
      filter(Genome_Site == g.site.tmp & Read_Site == r.site.tmp) %>%
      pull(n)
    if(is.integer(val.tmp) && length(val.tmp) == 0L){ #| is.numeric(val.tmp)
      print("No connections")
    } else {
      t.mat[i,j] <- val.tmp
    }
  }
}

## Make the matrix symmetric by adding the transverse of the asymmetric mat
t.sym <- t.mat + t(t.mat)


######################## plot it in circlize ################################################

nm = unique(unlist(dimnames(t.sym)))
group = structure(gsub("\\d", "", nm), names = nm)
#df for messing with names
group <- as.data.frame(nm, row.names = NULL) %>%
  rename("vent" = "nm")
group$group = group$vent
#adjust group names
group$group <- gsub(".*Lau_Basin.*","Lau Basin Plume", group$group)
group$group <- gsub(".*ELSC.*","Lau Basin Deposit", group$group)
group$group <- gsub(".*Brothers.*","Brothers Volcano", group$group)
group$group <- gsub("Guaymas_[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9]","Guaymas Basin Deposit", group$group)
group$group <- gsub("Guaymas_Basin","Guaymas Basin Plume", group$group)
group$group <- gsub("Guaymas_[0-9][0-9][0-9][0-9]","Guaymas Basin Deposit", group$group)
group$group <- gsub(".*Cayman.*","Mid-Cayman Rise", group$group)
group$group <- gsub(".*Axial.*","Axial Seamount", group$group)
group$group <- gsub(".*EPR.*","East Pacific Rise", group$group)
group$group <- gsub(".*MAR.*","Mid-Atlantic Ridge", group$group)
#adjust vent names
#group <- tibble::rownames_to_column(group, "vent")
group$vent <- gsub("Lau_Basin_","P_", group$vent)
group$vent <- gsub("ELSC_","D_", group$vent)
group$vent <- gsub("Brothers_","", group$vent)
group$vent <- gsub("Guaymas_Basin","P_Guaymas", group$vent)
group$vent <- gsub("Guaymas_","", group$vent)
group$vent <- gsub("Cayman_","", group$vent)
group$vent <- gsub("Axial_","", group$vent)
group$vent <- gsub("EPR_","", group$vent)
group$vent <- gsub("MAR_","", group$vent)
group$vent <- gsub("_"," ", group$vent)

m.connect_test <- t.sym
rownames(m.connect_test) <- paste(group$vent)
colnames(m.connect_test) <- paste(group$vent)

#vector again
group <- dplyr::pull(group)
group2 <- setNames(as.character(group$group), group$vent)


#set colors
grid.col = c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478", "East Pacific Rise" = "#CE9A28",
             "Guaymas Basin Deposit" = "#28827A", "Guaymas Basin Plume" = "#499e97",
             "Lau Basin Plume" = "#3F78C1", "Lau Basin Deposit" = "#72a0db",
             "Mid-Atlantic Ridge" = "#8c510a", "Mid-Cayman Rise" = "#000000",
             #distinct vent locations above
             "4281-140" = "#CE9A28", "4559-240" = "#28827A",
             "4561" = "#28827A", "4571-419" = "#28827A", 
             #Guaymas deposit
             "D Mariner" = "#72a0db", "D Abe" = "#72a0db", "D Vai Lili V2" = "#72a0db", 
             "D Tui Malila" = "#72a0db", "D Bowl" = "#72a0db",
             #Lau deposit
             "Deep" = "#000000", "Shallow" = "#000000", 
             #MCR
             "Diffuse" = "#B56478", "NWCA" = "#B56478", "LC" = "#B56478", 
             "NWCB" = "#B56478", "UC" = "#B56478", 
             #Brothers Volcano
             "Lucky" = "#8c510a", "Rainbow" = "#8c510a",
             #MAR
             "P Abe" = "#3F78C1", "P Kilo Moana" = "#3F78C1", "P Mariner" = "#3F78C1",
             "P Tahi Moana" = "#3F78C1", "P Tui Malila" = "#3F78C1", 
             #Lau plume
             "Seawater" = "#4F508C", "Plume" = "#4F508C", 
             #Axial Seamount
             "PIR-30" = "#CE9A28", 
             #EPR
             "P Guaymas" = "#499e97" 
             #Guaymas plume
) #specify colors of sites/outer ring

col_fun = "grey" #specify color of links between them

#reset before running
dev.off()
circos.clear()
svg("../Read_Mapping/Output/Circos_coverm_gd_3kb_70mincov_largeFont.svg")
#set font size
par(cex = 1.8, mar = c(0, 0, 0, 0))
#set gaps between blocks
#circos.par(gap.after = c(rep(.5, length(unique(colnames(t.mat_test)))))) #this is crucial to plot

#chordDiagram
chordDiagram(m.connect_test, 
             group = group2, 
             grid.col = grid.col,
             annotationTrack = c("grid"),
             annotationTrackHeight = c(.06, .06), #change height of inner group track
             link.border = "black",
             transparency = 0.3,
             symmetric = TRUE,
             big.gap = 8,
             col = col_fun, #turn this off to get colors of groups for bands
             preAllocateTracks = list(
               track.height = mm_h(6),
               track.margin = c(mm_h(1), 0)
             ))
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  # circos.text(CELL_META$xcenter,
  #             ylim[1] + cm_h(3), 
  #             labels = gsub("\\..*","", sector.name),
  #             col = "black", cex = 0.5, facing = "bending", 
  #             adj = c(0.5, 0.5), niceFacing = TRUE, font = 0.5)
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.5, niceFacing = TRUE,
              col = "white")
}, bg.border = NA)

highlight.sector(sector.index = c("P Abe", "P Kilo Moana",
                                  "P Mariner", "P Tahi Moana",
                                  "P Tui Malila"), track.index = 1, col = "#3F78C1", 
                 text = "Lau Basin Plume", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("D Mariner", "D Tui Malila",
                                  "D Bowl", "D Abe",
                                  "D Vai Lili V2"), track.index = 1, col = "#72a0db", 
                 text = "Lau Basin Deposit", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("NWCA", "NWCB", "LC",
                                  "UC", "Diffuse"), track.index = 1, col = "#B56478", 
                 text = "Brothers Volcano", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4561", "4571-419",
                                  "4559-240"), track.index = 1, col = "#28827A", 
                 text = "Guaymas Basin Deposit", cex = 0.6, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("P Guaymas"), track.index = 1, col = "#499e97", 
                 text = "Guaymas Basin Plume", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Deep", "Shallow"), track.index = 1, col = "#000000", 
                 text = "Mid-Cayman Rise", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Plume", "Seawater"), track.index = 1, col = "#4F508C", 
                 text = "Axial Seamount", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4281-140", "PIR-30"), track.index = 1, col = "#CE9A28", 
                 text = "East Pacific Rise", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Lucky", "Rainbow"), track.index = 1, col = "#8c510a", 
                 text = "Mid-Atlantic Ridge", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")

dev.off()


################### make circos plot of intra vents ##########################################

#remove coverm.simple_gd from coverm.simple and make new iv

coverm.simple_iv <- coverm.simple %>%
  anti_join(coverm.simple_gd)

## Make the Genome_Site x Read_Site Matrix with # of connections
######################### the following code was written by Spencer R. Keyser (skeyser@wisc.edu)
## make a zero-filled matrix
t.mat <- matrix(nrow = length(unique(coverm.simple_iv$Genome_Site)), ncol = length(unique(coverm.simple_iv$Read_Site)), data = 0)
colnames(t.mat) <- unique(coverm.simple_iv$Genome_Site)
rownames(t.mat) <- unique(coverm.simple_iv$Genome_Site)

## Holder dataframe
for(i in 1:nrow(t.mat)){
  g.site.tmp <- rownames(t.mat)[i]
  for(j in 1:ncol(t.mat)){
    r.site.tmp <- colnames(t.mat)[j]
    val.tmp <- coverm.simple_iv %>% #change here
      filter(Genome_Site == g.site.tmp & Read_Site == r.site.tmp) %>%
      pull(n)
    if(is.integer(val.tmp) && length(val.tmp) == 0L){ #| is.numeric(val.tmp)
      print("No connections")
    } else {
      t.mat[i,j] <- val.tmp
    }
  }
}

## Make the matrix symmetric by adding the transverse of the asymmetric mat
t.sym <- t.mat + t(t.mat)


######################## plot IV in circlize ################################################

nm = unique(unlist(dimnames(t.sym)))
group = structure(gsub("\\d", "", nm), names = nm)
#df for messing with names
group <- as.data.frame(nm, row.names = NULL) %>%
  rename("vent" = "nm")
group$group = group$vent
#adjust group names
group$group <- gsub(".*Lau_Basin.*","Lau Basin Plume", group$group)
group$group <- gsub(".*ELSC.*","Lau Basin Deposit", group$group)
group$group <- gsub(".*Brothers.*","Brothers Volcano", group$group)
group$group <- gsub("Guaymas_[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9]","Guaymas Basin Deposit", group$group)
group$group <- gsub("Guaymas_Basin","Guaymas Basin Plume", group$group)
group$group <- gsub("Guaymas_[0-9][0-9][0-9][0-9]","Guaymas Basin Deposit", group$group)
group$group <- gsub(".*Cayman.*","Mid-Cayman Rise", group$group)
group$group <- gsub(".*Axial.*","Axial Seamount", group$group)
group$group <- gsub(".*EPR.*","East Pacific Rise", group$group)
group$group <- gsub(".*MAR.*","Mid-Atlantic Ridge", group$group)
#adjust vent names
#group <- tibble::rownames_to_column(group, "vent")
group$vent <- gsub("Lau_Basin_","P_", group$vent)
group$vent <- gsub("ELSC_","D_", group$vent)
group$vent <- gsub("Brothers_","", group$vent)
group$vent <- gsub("Guaymas_Basin","P_Guaymas", group$vent)
group$vent <- gsub("Guaymas_","", group$vent)
group$vent <- gsub("Cayman_","", group$vent)
group$vent <- gsub("Axial_","", group$vent)
group$vent <- gsub("EPR_","", group$vent)
group$vent <- gsub("MAR_","", group$vent)
group$vent <- gsub("_"," ", group$vent)

m.connect_test <- t.sym
rownames(m.connect_test) <- paste(group$vent)
colnames(m.connect_test) <- paste(group$vent)

#vector again
group <- dplyr::pull(group)
group2 <- setNames(as.character(group$group), group$vent)


#set colors
grid.col = c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478", "East Pacific Rise" = "#CE9A28",
             "Guaymas Basin Deposit" = "#28827A", "Guaymas Basin Plume" = "#499e97",
             "Lau Basin Plume" = "#3F78C1", "Lau Basin Deposit" = "#72a0db",
             "Mid-Atlantic Ridge" = "#8c510a", "Mid-Cayman Rise" = "#000000",
             #distinct vent locations above
             "4281-140" = "#CE9A28", "4559-240" = "#28827A",
             "4561" = "#28827A", "4571-419" = "#28827A", 
             #Guaymas deposit
             "D Mariner" = "#72a0db", "D Abe" = "#72a0db", "D Vai Lili V2" = "#72a0db", 
             "D Tui Malila" = "#72a0db", "D Bowl" = "#72a0db",
             #Lau deposit
             "Deep" = "#000000", "Shallow" = "#000000", 
             #MCR
             "Diffuse" = "#B56478", "NWCA" = "#B56478", "LC" = "#B56478", 
             "NWCB" = "#B56478", "UC" = "#B56478", 
             #Brothers Volcano
             "Lucky" = "#8c510a", "Rainbow" = "#8c510a",
             #MAR
             "P Abe" = "#3F78C1", "P Kilo Moana" = "#3F78C1", "P Mariner" = "#3F78C1",
             "P Tahi Moana" = "#3F78C1", "P Tui Malila" = "#3F78C1", 
             #Lau plume
             "Seawater" = "#4F508C", "Plume" = "#4F508C", 
             #Axial Seamount
             "PIR-30" = "#CE9A28", 
             #EPR
             "P Guaymas" = "#499e97" 
             #Guaymas plume
) #specify colors of sites/outer ring

col_fun = "grey" #specify color of links between them

#reset before running
dev.off()
circos.clear()
svg("../Read_Mapping/Output/Circos_coverm_iv_3kb_70mincov_largeFont.svg")
#set font size
par(cex = 1.8, mar = c(0, 0, 0, 0))
#set gaps between blocks
#circos.par(gap.after = c(rep(.5, length(unique(colnames(t.mat_test)))))) #this is crucial to plot

#chordDiagram
chordDiagram(m.connect_test, 
             group = group2, 
             grid.col = grid.col,
             annotationTrack = c("grid"),
             annotationTrackHeight = c(.06, .06), #change height of inner group track
             link.border = "black",
             transparency = 0.3,
             symmetric = TRUE,
             big.gap = 8,
             col = col_fun, #turn this off to get colors of groups for bands
             preAllocateTracks = list(
               track.height = mm_h(6),
               track.margin = c(mm_h(1), 0)
             ))
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  # circos.text(CELL_META$xcenter,
  #             ylim[1] + cm_h(3), 
  #             labels = gsub("\\..*","", sector.name),
  #             col = "black", cex = 0.5, facing = "bending", 
  #             adj = c(0.5, 0.5), niceFacing = TRUE, font = 0.5)
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.5, niceFacing = TRUE,
              col = "white")
}, bg.border = NA)

highlight.sector(sector.index = c("P Abe", "P Kilo Moana",
                                  "P Mariner", "P Tahi Moana",
                                  "P Tui Malila"), track.index = 1, col = "#3F78C1", 
                 text = "Lau Basin Plume", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("D Mariner", "D Tui Malila",
                                  "D Bowl", "D Abe",
                                  "D Vai Lili V2"), track.index = 1, col = "#72a0db", 
                 text = "Lau Basin Deposit", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("NWCA", "NWCB", "LC",
                                  "UC", "Diffuse"), track.index = 1, col = "#B56478", 
                 text = "Brothers Volcano", cex = 0.8, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4561", "4571-419",
                                  "4559-240"), track.index = 1, col = "#28827A", 
                 text = "Guaymas Basin Deposit", cex = 0.6, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("P Guaymas"), track.index = 1, col = "#499e97", 
                 text = "Guaymas Basin Plume", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Deep", "Shallow"), track.index = 1, col = "#000000", 
                 text = "Mid-Cayman Rise", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Plume", "Seawater"), track.index = 1, col = "#4F508C", 
                 text = "Axial Seamount", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4281-140", "PIR-30"), track.index = 1, col = "#CE9A28", 
                 text = "East Pacific Rise", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("Lucky", "Rainbow"), track.index = 1, col = "#8c510a", 
                 text = "Mid-Atlantic Ridge", cex = 0.5, text.col = "white", niceFacing = TRUE,
                 facing = "bending")

dev.off()



######################## understanding the coverm filtered output ################################################

#make same file to manipulate
coverm.anno <- coverm.map
#creat PV mapping
pv_map <- as.data.frame(unique(sort(coverm.anno$Read_Site))) %>%
  rename("sites" = "unique(sort(coverm.anno$Read_Site))")
pv_map$PV <- c("Plume", "Plume", "Vent", "Vent",
               "Vent", "Vent", "Vent", "Plume",
               "Plume", "Vent", "Vent", "Vent",
               "Vent", "Vent", "Vent", "Vent",
               "Vent", "Vent", "Vent", "Plume",
               "Plume", "Plume", "Plume", "Plume",
               "Plume", "Vent", "Vent") 
#add this info to the large file of filtered read mapping
coverm.anno <- pv_map %>%
  dplyr::select('sites', 'PV') %>%
  right_join(coverm.anno, by = c("sites" = "Genome_Site"))
coverm.anno <- coverm.anno %>%
  rename("Genome_Site" = "sites")

coverm.anno <- pv_map %>%
  dplyr::select('sites', 'PV') %>%
  right_join(coverm.anno, by = c("sites" = "Read_Site"))
coverm.anno <- coverm.anno %>%
  rename("Read_Site" = "sites")

#add genome size
#read in gen size calculated with seqkit for accurate values
gensize <- read.table(file = "../../Seqkit_GenSize/VentPlume_seqkit_stats_renamed_space.txt",
                      header = TRUE)
coverm.anno <- gensize %>%
  dplyr::select('file', 'sum_len') %>%
  right_join(coverm.anno, by = c("file" = "Genome"))

coverm.anno.pv <- coverm.anno %>%
  filter(PV.x != PV.y) %>%
  filter(sum_len >= 3000)

coverm.anno.pv <- coverm.anno.pv %>%
  mutate(Read_Genome = paste(coverm.anno.pv$Read_Site, "_", coverm.anno.pv$Genome_Site))


############# plot all read mapping in circlize, red outline for GD ################################################





############################# unused #############################

# #convert from wide to long
# coverm_rc_long <- melt(coverm_rc, id = c("Genome"), na.rm = TRUE)
# #remove 0s
# coverm_rc_long <- coverm_rc_long %>%  
#   filter(value > 0) %>%
#   rename("Reads_Site" = "variable") %>%
#   rename("Reads_count" = "value")
# 
# coverm_rc_long <- coverm_rc_long %>%
#   filter(!str_detect(Genome, remove.list))

# #rename site names on coverm file
# coverm_rc_long$Reads_Site <- gsub(".Read.Count","", coverm_rc_long$Reads_Site)
# coverm_rc_long <- abun_names %>%
#   dplyr::select("V1", "V2") %>%
#   right_join(coverm_rc_long, by = c("V2" = "Reads_Site")) %>%
#   rename("Reads_Site" = "V1") %>%
#   rename("Reads" = "V2")
# 
# ######################## make genome site column ################################################
# 
# coverm_rc_long$Genome_Site <- coverm_rc_long$Genome 
# coverm_rc_long <- coverm_rc_long %>% 
#   separate(Genome_Site, c("Genome_Site", NA), sep = "_NODE|_scaffold|_vRhyme|_k95") %>%
#   select(c("Reads_Site", "Reads", "Genome_Site", "Genome", "Reads_count"))
# 
# #faster replace all the naming patterns
# #test <- coverm_rc_long
# coverm_rc_long$Genome_Site <- stri_replace_all_regex(coverm_rc_long$Genome_Site,
#                                                      pattern=c("_A[0-9]",
#                                                                "_T[0-9][0-9]", "_T[0-9]", "_S0[0-9][0-9]",
#                                                                "_S1[0-9][0-9]", "_[0-9][0-9][0-9]-[0-9][0-9][0-9]",
#                                                                "-38[0-9]"),
#                                                      replacement='',
#                                                      vectorize=FALSE)
# 
# coverm_rc_long$Genome_Site <- gsub("*_M1[0-9]","",coverm_rc_long$Genome_Site)
# coverm_rc_long$Genome_Site <- gsub("*_M[0-9]","",coverm_rc_long$Genome_Site)
# 
# #remove self to self matches
# coverm_rc_long <- coverm_rc_long %>%
#   filter(Reads_Site != Genome_Site)
# 
# ######################## generate the matrix ################################################
# 
# #Spencer
# m.connect <- matrix(ncol = length(unique(coverm_rc_long$Reads_Site)),
#                     nrow = length(unique(coverm_rc_long$Reads_Site)),
#                     data = 0)
# rownames(m.connect) <- unique(coverm_rc_long$Reads_Site)
# colnames(m.connect) <- unique(coverm_rc_long$Reads_Site)
# 
# for(i in 1:length(unique(coverm_rc_long$Reads_Site))){
#   site.tmp <- unique(coverm_rc_long$Reads_Site)[i]
#   tmp.df <- coverm_rc_long %>%
#     filter(Reads_Site == site.tmp)
#   for(j in 1:length(unique(tmp.df$Genome_Site))){
#     site.tmp2 <- unique(tmp.df$Genome_Site)[j]
#     tmp.df2 <- tmp.df %>%
#       filter(Genome_Site == site.tmp2)
#     con.ct <- nrow(tmp.df2)
#     m.connect[rownames(m.connect) == site.tmp, colnames(m.connect) == site.tmp2] <- con.ct
#   } #j
# } #i
# 
# test.spencer <- coverm_rc_long %>% 
#   filter(Reads_Site %in% c("Brothers_Diffuse", "Brothers_LC")) %>%
#   filter(Genome_Site %in% c("Brothers_Diffuse", "Brothers_LC"))


# ######################### old chaotic way ######################### 


# ######################### repeating Circos plot with coverm read mapping data ######################### 
# 
# coverm_rc <- read.delim2(file = "../VentVirus_Analysis/output/coverm_readCount_top50.tsv")
# 
# #make site names the same between reads and genome
# coverm_rc <- coverm_rc %>% separate(Reads_Site, c("Reads_Site", NA),
#                                     sep = "_min")
# #remove weird spaces
# coverm_rc$Reads_Site <- gsub(" ","", coverm_rc$Reads_Site) #remove spaces
# coverm_rc$Genome_Site <- gsub(" ","", coverm_rc$Genome_Site)
# 
# #drop rows where read and genome site are the same
# coverm_rc <- coverm_rc %>% 
#   mutate(duplicates = if_else(Reads_Site == Genome_Site,
#                               TRUE,
#                               FALSE)) %>% 
#   filter(duplicates == FALSE)
# 
# coverm_rc <- coverm_rc %>% select(Reads_Site, Genome_Site, Genome, CovermReadCount)
# 
# #test to see number of combos
# coverm_rc.temp <- coverm_rc %>%
#   filter(grepl("Axial", Reads_Site)) %>%
#   filter(grepl("Guaymas", Genome_Site))
# 
# coverm_rc_long <- coverm_rc %>% 
#   pivot_longer(cols = c('Reads_Site', 'Genome_Site')) %>%
#   rename(value, "Site"= "value")
# 
# ################################ create intra vent file ##########################################
# 
# # #replace strings with Vent or Plume
# # coverm_rc_long$Site2 <- gsub(".*Lau_Basin.*","Plume",coverm_rc_long$Site) #the placement of the periods is crucial for replacing whole string
# # coverm_rc_long$Site2 <- gsub(".*Cayman.*","Plume",coverm_rc_long$Site2)
# # coverm_rc_long$Site2 <- gsub(".*Guaymas_Basin.*","Plume",coverm_rc_long$Site2)
# # coverm_rc_long$Site2 <- gsub(".*Axial.*","Plume",coverm_rc_long$Site2)
# # #separate the plume and vent dataframes to make life easier
# # temp_plume <- coverm_rc_long %>% filter(Site2 == "Plume")
# # #remove Plume from mcl clusters so can make all names left Vent
# # coverm_rc_long <- coverm_rc_long %>% filter(Site2 != "Plume")
# # coverm_rc_long$Site2 <- "Vent"
# # #Now put them back together
# # coverm_rc_long <- rbind(coverm_rc_long, temp_plume)
# 
# #number of viruses from Plume and Vent:
# table(coverm_rc_long$Site)
# #406 plume and 1042 vent with skani parameters, 3kb, and 50AF
# 
# #group by cluster, count occurrences of Site
# temp_count <- coverm_rc_long %>% group_by(Genome,CovermReadCount) %>% count(Site2)
# 
# #see if any cluster now occurs twice FOR SPENCER
# temp_count.temp <- temp_count %>%
#   group_by(Genome, CovermReadCount) %>%
#   tally() %>%
#   filter(n == 2) %>%
#   mutate(unique_id = paste0(Genome, "_", CovermReadCount))
# 
# #get clusters that are vent and plume
# coverm_rc_long_gd <- coverm_rc_long %>%
#   mutate(unique_id = paste0(Genome, "_", CovermReadCount)) %>%
#   filter(unique_id %in% temp_count.temp$unique_id)
# 
# 
# ######################### creating input matrix coverm ################################################
# ######################### intra vent coverm #####################################
# 
# #sort the coverm_long_gd in alph order to avoid headaches
# coverm_rc_long_gd <- coverm_rc_long_gd %>%
#   group_by(CovermReadCount,Genome) %>%
#   arrange(Site2, .by_group = T) %>%
#   ungroup()
# 
# test2 <- coverm_rc_long_gd %>%
#   select(Site,unique_id) %>%
#   group_by(Site,unique_id) %>%
#   distinct()
# 
# #mod names to be less specific for gd
# #adjust group names
# test2$Site <- gsub(".*Lau_Basin.*","Lau_Basin_Plume", test2$Site)
# test2$Site <- gsub(".*ELSC.*","Lau_Basin_Deposit", test2$Site)
# test2$Site <- gsub(".*Brothers.*","Brothers_Volcano", test2$Site)
# test2$Site <- gsub(".*Guaymas.*","Guaymas_Basin", test2$Site)
# test2$Site <- gsub(".*Cayman.*","Mid-Cayman_Rise", test2$Site)
# test2$Site <- gsub(".*Axial.*","Axial_Seamount", test2$Site)
# test2$Site <- gsub(".*EPR.*","East_Pacific_Rise", test2$Site)
# test2$Site <- gsub(".*MAR.*","Mid-Atlantic_Ridge", test2$Site)
# 
# #add site to the unique_id to make it truly unique
# 
# #length(unique(test$Genome))
# 
# ######################### the following code was written by Spencer R. Keyser (skeyser@wisc.edu)
# ## make a zero-filled matrix
# t.mat2 <- matrix(nrow = length(unique(test2$Site)), ncol = length(unique(test2$Site)), data = NA)
# colnames(t.mat2) <- unique(test2$Site)
# rownames(t.mat2) <- unique(test2$Site)
# 
# ## Holder dataframe
# df.hold2 <- data.frame(from = NA, to = NA, unique_id = NA)
# cluster.list2 <- unique(test2$unique_id)
# for(i in 1:length(cluster.list2)){
#   clust.tmp2 <- cluster.list2[i]
#   id.tmp2 <- test2[test2$unique_id == clust.tmp2, ]
#   pw2 <- apply(combn(id.tmp2$Site,2),2,paste,collapse=' ')
#   pw2 <- data.frame(SitePair = pw2)
#   pw2 <- tidyr::separate(pw2, col = SitePair, into = c("from", "to"), sep = " ")
#   pw2$unique_id <- rep(unique(id.tmp2$unique_id), nrow(pw2))
#   df.hold2 <- rbind(df.hold2, pw2)
# }
# 
# df.hold2 <- df.hold2[complete.cases(df.hold2),]
# 
# ## Take the pairs and calculate the number of times they are the same
# df.hold2$Unique <- paste0(df.hold2$from, "__", df.hold2$to)
# 
# test2 <- df.hold2 %>%
#   group_by(Unique) %>%
#   count() %>%
#   separate(col = Unique, into = c("from", "to"), sep = "__")
# 
# for(i in 1:nrow(test2)){
#   tmp2 <- test2[i,]
#   to.tmp2 <- test2$to[i]
#   from.tmp2 <- test2$from[i]
#   val.tmp2 <- test2$n[i]
#   t.mat2[colnames(t.mat2) == from.tmp2, rownames(t.mat2) == to.tmp2] <- val.tmp2
#   t.mat2[colnames(t.mat2) == to.tmp2, rownames(t.mat2) == from.tmp2] <- val.tmp2
#   
# }
# 
# t.mat2[is.na(t.mat2)] <- 0
# 























