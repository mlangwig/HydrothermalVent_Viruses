#install.packages("circlize")
library(circlize)
library(tidyverse)
library(dplyr)

########################### major inputs ################################################

#sort the ani_long_gd in alph order to avoid headaches
ani_long_meta_gd <- ani_long_meta_gd %>%
  group_by(id) %>%
  arrange(Site, .by_group = T)

######################### creating inpit matrix ################################################

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

#create the matrix
mat <- as.data.frame.matrix(xtabs(Freq ~ Site1 + Site2, df))
mat <- as.matrix(mat)

######################### visualize ################################################

dev.off()
circos.clear()
#set font size
par(cex = .5, mar = c(0, 0, 0, 0))
#set gaps between blocks
circos.par(gap.after = c(rep(10, length(unique(colnames(mat)))))) #this is crucial to plot

#set colors
grid.col = c("Axial_Seamount" = "#4F508C", "Brothers_Volcano" = "#B56478",
             "Guaymas_Basin" = "#28827A", "Lau_Basin" = "#3F78C1",
            "Mid_Atlantic_Ridge" = "#8c510a", "Mid_Cayman_Rise" = "#000000") #specify colors of sites/outer ring

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
  circos.text(CELL_META$xcenter, 
              ylim[1] + cm_h(2), 
              sector.name, 
              facing = "clockwise",
              niceFacing = TRUE, 
              adj = c(0.5, 0.5),
              cex = 1,
              font = 2)
  # circos.axis(h = "bottom",
  #             labels.cex = .6,
  #             sector.index = sector.name
  # )
}, bg.border = NA)


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
