#install.packages("circlize")
library(circlize)
library(tidyverse)
library(dplyr)

#sort the ani_long_gd in alph order to avoid headaches
ani_long_meta_gd <- ani_long_meta_gd %>%
  group_by(id) %>%
  arrange(Site, .by_group = T)

#transform ani_long_meta_gd from ANI_clust.R to input for circlize
circ_mat <- ani_long_meta_gd %>%
  select(c("id", "Site")) %>%
  unique() %>%
  group_by(id) %>%
  #count(Site)
  mutate(Site = toString(Site)) %>%
  unique()

#get the counts of unique strings
circ_mat2 <- as.data.frame(table(circ_mat$Site))

#separate strings for summing same site match in different order
circ_mat2 <- circ_mat2 %>% separate(Var1, c("Site1", "Site2"), 
                                  sep= ",")

mat <- circ_mat2 %>%
  pivot_wider(names_from = "Site1", values_from = "Freq") #%>%
  remove_rownames() %>%
  column_to_rownames(var="Site2") %>%
  replace(is.na(.), 0)

mat <- read.delim2(file = "Input/circ_input.tsv", sep = "\t")
mat <- mat %>%
  remove_rownames() %>%
  column_to_rownames(var="X")
mat <- as.matrix(mat)


#visualize
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
