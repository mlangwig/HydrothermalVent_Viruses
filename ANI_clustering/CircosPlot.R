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
# png("Output/Circos_gd.png", res = 300, width = 350, height = 350, units="mm", 
#     pointsize = 13)
svg("Output/Circos_gd.svg")
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

#sort the ani_long_iv in alphabetical order to avoid headaches
ani_long_meta_iv <- ani_long_meta_iv %>%
  group_by(id) %>%
  arrange(Site, .by_group = T)

test <- ani_long_meta_iv %>%
  select(Site, id) %>%
  group_by(Site, id) %>%
  distinct()

length(unique(test$id))

######################### the following code was written by Spencer R. Keyser (skeyser@wisc.edu)
## make a zero-filled matrix
t.mat <- matrix(nrow = length(unique(test$Site)), ncol = length(unique(test$Site)), data = NA)
colnames(t.mat) <- unique(test$Site)
rownames(t.mat) <- unique(test$Site)

## Holder dataframe
df.hold <- data.frame(from = NA, to = NA, id = NA)
cluster.list <- unique(test$id)
for(i in 1:length(cluster.list)){
  clust.tmp <- cluster.list[i]
  id.tmp <- test[test$id == clust.tmp, ]
    pw <- apply(combn(id.tmp$Site,2),2,paste,collapse=' ')
    pw <- data.frame(SitePair = pw)
    pw <- tidyr::separate(pw, col = SitePair, into = c("from", "to"), sep = " ")
    pw$id <- rep(unique(id.tmp$id), nrow(pw))
    df.hold <- rbind(df.hold, pw)
  }

df.hold <- df.hold[complete.cases(df.hold),]

## Take the pairs and calculate the number of times they are the same
df.hold$Unique <- paste0(df.hold$from, "__", df.hold$to)

test <- df.hold %>%
  group_by(Unique) %>%
  count() %>%
  separate(col = Unique, into = c("from", "to"), sep = "__")

for(i in 1:nrow(test)){
  tmp <- test[i,]
  to.tmp <- test$to[i]
  from.tmp <- test$from[i]
  val.tmp <- test$n[i]
  t.mat[colnames(t.mat) == from.tmp, rownames(t.mat) == to.tmp] <- val.tmp
  t.mat[colnames(t.mat) == to.tmp, rownames(t.mat) == from.tmp] <- val.tmp
  
}

t.mat[is.na(t.mat)] <- 0

#############

######################### plotting ################################################
######################### intra vent #####################################
nm = unique(unlist(dimnames(t.mat)))
group = structure(gsub("\\d", "", nm), names = nm)
#df for messing with names
group <- as.data.frame(group, row.names = NULL)
#adjust group names
group$group <- gsub(".*Lau_Basin.*","Lau Basin Plume", group$group)
group$group <- gsub(".*ELSC.*","Lau Basin Deposit", group$group)
group$group <- gsub(".*Brothers.*","Brothers Volcano", group$group)
group$group <- gsub(".*Guaymas.*","Guaymas Basin", group$group)
group$group <- gsub(".*Cayman.*","Mid-Cayman Rise", group$group)
group$group <- gsub(".*Axial.*","Axial Seamount", group$group)
group$group <- gsub(".*EPR.*","East Pacific Rise", group$group)
#adjust vent names
group <- tibble::rownames_to_column(group, "vent")
group$vent <- gsub("Lau_Basin_","P_", group$vent)
group$vent <- gsub("ELSC_","D_", group$vent)
group$vent <- gsub("Brothers_","", group$vent)
group$vent <- gsub("Guaymas_","", group$vent)
group$vent <- gsub("Cayman_","", group$vent)
group$vent <- gsub("Axial_","", group$vent)
group$vent <- gsub("EPR_","", group$vent)
group$vent <- gsub("_"," ", group$vent)

t.mat_test <- t.mat
rownames(t.mat_test) <- paste(group$vent)
colnames(t.mat_test) <- paste(group$vent)

#vector again
#group <- dplyr::pull(group)
group2 <- setNames(as.character(group$group), group$vent)


#set colors
grid.col = c("Axial Seamount" = "#4F508C", "Brothers Volcano" = "#B56478",
             "Guaymas Basin" = "#28827A", "Lau Basin" = "#3F78C1",
             "Mid Atlantic Ridge" = "#8c510a", "Mid Cayman Rise" = "#000000",
             "Seawater" = "#4F508C", "Plume" = "#4F508C", #axial sites
             "Diffuse" = "#B56478", "NWCA" = "#B56478", 
             "NWCB" = "#B56478", "UC" = "#B56478", #brothers sites
             "4561" = "#28827A", "4571-419" = "#28827A", "4559-240" = "#28827A", #Guaymas
             "P Abe" = "#3F78C1", "P Kilo Moana" = "#3F78C1", "P Mariner" = "#3F78C1",
             "P Tahi Moana" = "#3F78C1", "P Tui Malila" = "#3F78C1", #Lau plume
             "D Mariner" = "#72a0db", "D Abe" = "#72a0db", "D Vai Lili V2" = "#72a0db", #Lau deposit
             "D Tui Malila" = "#72a0db", "D Bowl" = "#72a0db",
             "Deep" = "#000000", "Shallow" = "#000000", #MCR
             "4281-140" = "#CE9A28", "PIR-30" = "#CE9A28" #EPR
             ) #specify colors of sites/outer ring

col_fun = "grey" #specify color of links between them

#reset before running
dev.off()
circos.clear()
svg("Output/Circos_iv.svg")
#set font size
par(cex = 1, mar = c(0, 0, 0, 0))
#set gaps between blocks
#circos.par(gap.after = c(rep(.5, length(unique(colnames(t.mat_test)))))) #this is crucial to plot

#chordDiagram
chordDiagram(t.mat_test, 
             group = group2, 
             grid.col = grid.col,
             annotationTrack = c("grid"),
             link.border = "black",
             transparency = 0.3,
             symmetric = TRUE,
             big.gap = 5,
             col = col_fun,
             preAllocateTracks = list(
               track.height = mm_h(4),
               track.margin = c(mm_h(2), 0)
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
highlight.sector(sector.index = c("NWCA", "NWCB",
                                  "UC", "Diffuse"), track.index = 1, col = "#B56478", 
                 text = "Brothers Volcano", cex = 0.6, text.col = "white", niceFacing = TRUE,
                 facing = "bending")
highlight.sector(sector.index = c("4561", "4571-419",
                                  "4559-240"), track.index = 1, col = "#28827A", 
                 text = "Guaymas Basin", cex = 0.8, text.col = "white", niceFacing = TRUE,
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

dev.off()

####################### unused

# circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#   xlim = get.cell.meta.data("xlim")
#   ylim = get.cell.meta.data("ylim")
#   sector.name = get.cell.meta.data("sector.index")
#   
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#               facing = "inside", niceFacing = TRUE, adj = c(0, 0.5), font = 1)
# }, bg.border = NA)

# temp <- df %>%
#   pivot_wider(names_from = "Site1", values_from = "Freq") #%>%
#   remove_rownames() %>%
#   column_to_rownames(var="Site2") %>%
#   replace(is.na(.), 0)

#mat <- read.delim2(file = "Input/circ_input.tsv", sep = "\t")
# mat <- mat %>%
#   remove_rownames() %>%
#   column_to_rownames(var="X")

# #transform ani_long_meta_gd from ANI_clust.R to input for circlize
# circ_mat <- ani_long_meta_iv %>%
#   select(c("id", "Site")) %>%
#   unique() %>%
#   group_by(id) %>%
#   #count(Site)
#   mutate(Site = toString(Site)) %>%
#   unique() %>%
#   ungroup()
# 
# # #get the counts of unique strings
# # circ_mat2 <- as.data.frame(table(circ_mat$Site))
# # 
# # #separate strings for summing same site match in different order
# # circ_mat2 <- circ_mat2 %>% separate(Var1, c("Site1", "Site2", "Site3", "Site4", "Site5"), 
# #                                     sep= ",")
# 
# #remove weird spaces
# circ_mat2$Site1 <- gsub(" ","", circ_mat2$Site1) #remove spaces
# circ_mat2$Site2 <- gsub(" ","", circ_mat2$Site2) #remove spaces
# circ_mat2$Site3 <- gsub(" ","", circ_mat2$Site3) #remove spaces
# circ_mat2$Site4 <- gsub(" ","", circ_mat2$Site4) #remove spaces
# circ_mat2$Site5 <- gsub(" ","", circ_mat2$Site5) #remove spaces
# 
# 
# #tidy complete Site1 first
# df1 <- circ_mat2 %>%
#   complete(Site1, Site2, fill = list(Freq = 0))
# #get the ones with values and flip Site1 and Site2 order
# df_temp <- df1 %>% 
#   filter(Freq > 0) %>%
#   select(c("Site2", "Site1", "Freq")) %>%
#   rename("Site" = "Site2") %>%
#   rename("Site2" = "Site1") %>%
#   rename("Site1" = "Site")
# df1 <- bind_rows(df1, df_temp)
# 
# #tidy complete Site2 first, then switch the name order and drop repeats
# df2 <- circ_mat2 %>%
#   complete(Site2, Site1, fill = list(Freq = 0)) %>%
#   rename("Site" = "Site2") %>%
#   rename("Site2" = "Site1") %>%
#   rename("Site1" = "Site") %>%
#   filter(Freq < 1)
# 
# #combine the two complete test dfs, sort uniq, and complete again for full one
# df<-bind_rows(df1, df2) %>%
#   unique() %>%
#   complete(Site1, Site2, fill = list(Freq = 0))
# 
# #create the matrix
# mat <- as.data.frame.matrix(xtabs(Freq ~ Site1 + Site2, df))
# mat <- as.matrix(mat)


