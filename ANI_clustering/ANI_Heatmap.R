library(pheatmap)
library(reshape2)

#rm envmt
#rm(list = ls())

#read input
mat<-read.delim2("Input/skani_ani_sulfur_FullMatrix_renamed_AFvalues.txt", header = TRUE)
mat_df<-as.data.frame(mat)

mat_df_2 <- mat_df[,-1]
rownames(mat_df_2) <- mat_df[,1]

#convert data to numeric
mat_df_2 <- mat_df_2 %>% mutate_at(1:28, as.numeric)


###### Filter out rows where the only ANI is 100 to itself #####################

r <- rowSums(mat_df_2) #get sum of rows
rind <- which(r > 100) #get those that have sum >100
csum <- colSums(mat_df_2) 
cind <- which(csum > 100) #same for cols as above
ani_sulfur_filter <- mat_df_2[rind, cind] #filter for rows and cols that were >100 and make new df

###stats
min(ani_sulfur_filter[ani_sulfur_filter > 0])
# Lowest ANI among viruses without hit to only itself is 82.54% not caring about % aligned fraction
# This number is 49.75% when only results for ≥50% aligned fraction are used + correcting ANI value for
# aligned fraction

# 80 sulfur viruses that have ≥82.54% ANI to each other
# X number are between geographically distinct vent sites

# 26 viruses that have ≥95.09% ANI to each other and ≥50% AF
# X number between geographically distinct vents

#write the table
# write.table(ani_sulfur_filter, file = "Output/skani_matrix_SulfurViruses_50comp_filtered.txt", quote = FALSE,
#             col.names = NA, sep = "\t") #colnames = NA to keep first blank cell

#remove .fasta from row and col names
colnames(ani_sulfur_filter) = gsub(".fasta", "", colnames(ani_sulfur_filter))
rownames(ani_sulfur_filter) = gsub(".fasta", "", rownames(ani_sulfur_filter))

##Add Metadata Labels
metadata<-read.csv("Input/metadata_filtered.csv")
#metadata<-metadata %>% select("Virus", "Host")
rownames(metadata) <- metadata[,1]
metadata_tax<-as.data.frame(metadata[,c(2,3,4)])
#rownames(metadata_tax) <- metadata[,1]
#colnames(metadata_tax) <- c("Class")

#rename column
metadata_tax <- rename(metadata_tax, "Class" = "Host")

metadata_tax <- metadata_tax %>%
  mutate(Class = str_replace(Class, "c__", "")) %>%
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
   na_if('')
metadata_tax$Genus <- metadata_tax$Genus %>% replace_na("Unknown")

#this switches them between no blanks for rows or cols?
metadata_tax2 <- metadata_tax
rownames(metadata_tax) = colnames(ani_sulfur_filter)
rownames(metadata_tax2) = rownames(ani_sulfur_filter)

# library(randomcoloR)
# n <- 8
# palette <- rev(distinctColorPalette(n))
# 
# my_colour = list(Class = palette)

my_colour = list(
  Class = c(Alphaproteobacteria = '#8c510a', Aquificae = "#762a83", 
            Campylobacteria = "#234673", Dissulfuribacteria = '#43a2ca', Gammaproteobacteria = "#DE6B7E"), # , Bacteroidia = "#e7d4e8", Thermococci = "#e0f3db", Thermoproteia = '#99d8c9' 
  Genus = c("Caminibacter" = "#800000", "Dissulfuribacter" = "#ffd8b1", "GCA-2747325" = "#f58231", 
            "JAADCZ01" = "#fffac8", "Marinosulfonomonas" = "#808000", "Nautilia" = "#42d4f4", 
            "Persephonella" = "#aaffc3", "Sulfurovum" = "#FA8072", "SZUA-116" = "#bfef45", "SZUA-152" = "#4363d8", 
            "Thermopetrobacter" = "#a9a9a9",  "UBA2013" = "#469990", "Unknown" = "#dcbeff"), 
  Completeness = c("50.28" = "#ffffff", "75.42" = "#99d8c9", "100.0" = "#2ca25f")) 

dev.off()
#heatmap
yep_cor<-pheatmap::pheatmap(ani_sulfur_filter, 
                            fontsize = 6, 
                            fontsize_row = 5, 
                            fontsize_col = 5,
                            show_colnames = TRUE,
                            show_rownames = TRUE,
                            annotation_row = metadata_tax2,
                            annotation_col = metadata_tax,
                            clustering_method = "ward.D",
                            annotation_colors = my_colour,
                            cellwidth = 10,
                            cellheight = 10,
                            border_color = "grey",
                            display_numbers = round(ani_sulfur_filter,1),
                            number_color = "black")


#export
save_pheatmap_png <- function(x, filename, width=2700, height=2300, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(yep_cor, "Output/sulfur_VirusANI_pheatmap.png")

##pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(yep_cor, "Output/sulfur_VirusANI_pheatmap.pdf")


dev.off()


#save as svg
ggsave(filename = "Output/sulfur_VirusANI_pheatmap.svg", plot = yep_cor, height = 14)


