#install.packages("circlize")
library(circlize)
library(tidyverse)
library(dplyr)

########################### major inputs ################################################

#sort the ani_long_gd in alph order to avoid headaches
ani_long_meta_gd <- ani_long_meta_gd %>%
  group_by(id) %>%
  arrange(Site, .by_group = T)

######################### creating input matrix ################################################
######################### geo distinct #####################################

#transform ani_long_meta_gd from ANI_clust.R to input for circlize
circ_mat <- ani_long_meta_gd %>%
  select(c("id", "Site")) %>%
  unique() %>%
  group_by(id) %>%
  #count(Site)
  mutate(Site = toString(Site)) %>%
  unique() %>%
  ungroup()

#get the counts of unique strings
circ_mat2 <- as.data.frame(table(circ_mat$Site))

#separate strings for summing same site match in different order
circ_mat2 <- circ_mat2 %>% separate(Var1, c("Site1", "Site2"), 
                                  sep= ",")
#remove weird spaces
circ_mat2$Site1 <- gsub(" ","", circ_mat2$Site1) #remove spaces
circ_mat2$Site2 <- gsub(" ","", circ_mat2$Site2) #remove spaces

#tidy complete Site1 first
df1 <- circ_mat2 %>%
  complete(Site1, Site2, fill = list(Freq = 0))
#get the ones with values and flip Site1 and Site2 order
df_temp <- df1 %>% 
  filter(Freq > 0) %>%
  select(c("Site2", "Site1", "Freq")) %>%
  rename("Site" = "Site2") %>%
  rename("Site2" = "Site1") %>%
  rename("Site1" = "Site")
df1 <- bind_rows(df1, df_temp)

#tidy complete Site2 first, then switch the name order and drop repeats
df2 <- circ_mat2 %>%
  complete(Site2, Site1, fill = list(Freq = 0)) %>%
  rename("Site" = "Site2") %>%
  rename("Site2" = "Site1") %>%
  rename("Site1" = "Site") %>%
  filter(Freq < 1)

#combine the two complete test dfs, sort uniq, and complete again for full one
df<-bind_rows(df1, df2) %>%
  unique() %>%
  complete(Site1, Site2, fill = list(Freq = 0))  
df$Site1 <- gsub("_"," ", df$Site1) #add spaces
df$Site2 <- gsub("_"," ", df$Site2) #add spaces

#create the matrix
mat <- as.data.frame.matrix(xtabs(Freq ~ Site1 + Site2, df))
mat <- as.matrix(mat)

######################### visualize ################################################


dev.off()
png("Output/Circos_gd.png", res = 300, width = 350, height = 350, units="mm", 
    pointsize = 13)
#pdf(file = "Output/Circos_gd.pdf")
#tiff("Output/Circos_gd.tiff", width = 3, height = 3, units = 'in', res = 300)
circos.clear()
#set font size
par(cex = 1, mar = c(0, 0, 0, 0))
#set gaps between blocks
circos.par(gap.after = c(rep(15, length(unique(colnames(mat)))))) #this is crucial to plot

#set colors
grid.col = c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478",
             "Guaymas Basin" = "#28827A", "Lau Basin" = "#3F78C1",
            "Mid Atlantic Ridge" = "#8c510a", "Mid Cayman Rise" = "#000000") #specify colors of sites/outer ring

col_fun = "grey" #specify color of links between them

chordDiagram(as.matrix(mat),
             grid.col = grid.col,
             col = col_fun,
             annotationTrack = "grid", 
             #symmetric = TRUE,
             preAllocateTracks = 1,
             link.target.prop = FALSE,
             symmetric = TRUE,
             link.border = "black",
             transparency = 0.3)
             # order = c("Axial_Seamount", "Brothers_Volcano", "Guaymas_Basin",
             #           "Lau_Basin", "Mid_Atlantic_Ridge", "Mid_Cayman_Rise"),
             #directional = 1,
             #direction.type = c("diffHeight", "arrows"),
             #link.arr.type = "big.arrow", diffHeight = -mm_h(2))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), font = 2) # 
  # circos.axis(h = "bottom",
  #             labels.cex = .6,
  #             sector.index = sector.name
  # )
}, bg.border = NA)

dev.off()

#dev.copy(png,"Output/Circos_gd.png")
## all of your plotting code
#dev.off()

######################### creating input matrix ################################################
######################### intra vent #####################################

#sort the ani_long_gd in alph order to avoid headaches
ani_long_meta_iv <- ani_long_meta_iv %>%
  group_by(id) %>%
  arrange(Site, .by_group = T)

#transform ani_long_meta_gd from ANI_clust.R to input for circlize
circ_mat <- ani_long_meta_iv %>%
  select(c("id", "Site")) %>%
  unique() %>%
  group_by(id) %>%
  #count(Site)
  mutate(Site = toString(Site)) %>%
  unique() %>%
  ungroup()

#get the counts of unique strings
circ_mat2 <- as.data.frame(table(circ_mat$Site))

#separate strings for summing same site match in different order
circ_mat2 <- circ_mat2 %>% separate(Var1, c("Site1", "Site2"), 
                                    sep= ",")
#remove weird spaces
circ_mat2$Site1 <- gsub(" ","", circ_mat2$Site1) #remove spaces
circ_mat2$Site2 <- gsub(" ","", circ_mat2$Site2) #remove spaces

#tidy complete Site1 first
df1 <- circ_mat2 %>%
  complete(Site1, Site2, fill = list(Freq = 0))
#get the ones with values and flip Site1 and Site2 order
df_temp <- df1 %>% 
  filter(Freq > 0) %>%
  select(c("Site2", "Site1", "Freq")) %>%
  rename("Site" = "Site2") %>%
  rename("Site2" = "Site1") %>%
  rename("Site1" = "Site")
df1 <- bind_rows(df1, df_temp)

#tidy complete Site2 first, then switch the name order and drop repeats
df2 <- circ_mat2 %>%
  complete(Site2, Site1, fill = list(Freq = 0)) %>%
  rename("Site" = "Site2") %>%
  rename("Site2" = "Site1") %>%
  rename("Site1" = "Site") %>%
  filter(Freq < 1)

#combine the two complete test dfs, sort uniq, and complete again for full one
df<-bind_rows(df1, df2) %>%
  unique() %>%
  complete(Site1, Site2, fill = list(Freq = 0))

#create the matrix
mat <- as.data.frame.matrix(xtabs(Freq ~ Site1 + Site2, df))
mat <- as.matrix(mat)


####################### unused

# temp <- df %>%
#   pivot_wider(names_from = "Site1", values_from = "Freq") #%>%
#   remove_rownames() %>%
#   column_to_rownames(var="Site2") %>%
#   replace(is.na(.), 0)

#mat <- read.delim2(file = "Input/circ_input.tsv", sep = "\t")
# mat <- mat %>%
#   remove_rownames() %>%
#   column_to_rownames(var="X")
