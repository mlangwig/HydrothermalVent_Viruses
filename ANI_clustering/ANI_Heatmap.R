library(pheatmap)
library(reshape2)

#rm envmt
#rm(list = ls())

#read input
mat<-read.delim2("Input/skani_ani_matrix_full_SulfurViruses_50comp_renamed.txt", header = TRUE)
mat_df<-as.data.frame(mat)

mat_df_2 <- mat_df[,-1]
rownames(mat_df_2) <- mat_df[,1]

#convert data to numeric
mat_df_2 <- mat_df_2 %>% mutate_at(2:322, as.numeric)

##Add Metadata Labels
metadata<-read.csv("metadata2.csv")
rownames(metadata) <- metadata[,1]
metadata_phylum<-as.data.frame(metadata[,c(2)])
rownames(metadata_phylum) <- metadata[,1]
colnames(metadata_phylum) <- c("Phylum")

my_colour = list(
  Phylum = c(Korarchaeota = "#234673", Banfieldarchaeota = "#B9707E", 
             Crenarchaeota = '#A4D5B2'))

dev.off()
#heatmap
yep_cor<-pheatmap::pheatmap(mat_df_2, 
                            fontsize = 6, 
                            #fontsize_row = 5, 
                            #fontsize_col = 5,
                            show_colnames = FALSE,
                            show_rownames = FALSE,
                            #annotation_row = metadata_phylum,
                            #annotation_col = metadata_phylum,
                            clustering_method = "ward.D",
                            #annotation_colors = my_colour,
                            cellwidth = 2,
                            cellheight = 2,
                            border_color = "grey")

#export
save_pheatmap_png <- function(x, filename, width=1400, height=1200, res = 200) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(yep_cor, "NovelArchaea_pheatmap.png")

##pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(yep_cor, "NovelArchaea_pheatmap.pdf")


dev.off()
