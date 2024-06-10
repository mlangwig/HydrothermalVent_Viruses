#################################################### Beta Diversity test ######################################################

library(tidyverse)
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(phyloseq)
library(vegan)
library(csv)
library(ape)
#install.packages('devtools')
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
set.seed(81)

#################################################### create input ######################################################

######## abundance matrix aka OTU matrix
abund_norm_wide <- abund_long_norm
abund_norm_wide$Site <- gsub("_metaspades_scaffolds.min1000.fasta","", abund_norm_wide$Site) 
abund_norm_wide$Site <- gsub("min1000","Plume", abund_norm_wide$Site) 
abund_norm_wide$Site <- gsub("_scaffolds_Plume","", abund_norm_wide$Site)
abund_norm_wide$Site <- gsub(".fasta","", abund_norm_wide$Site)  
  
abund_norm_wide <- abund_norm_wide %>%  
  select(Genome, Site, abun_norm) %>%
  pivot_wider(names_from = Site, values_from = abun_norm, values_fill = list(abun_norm = 0)) %>%
  filter(!(Genome == "Lau_Basin_Tahi_Moana_vRhyme_bin_115")) #remove bin 115 because we suspect problems with read mapping/read duplication

#convert abund_norm_wide to data_table
abund_norm_wide <- as.data.table(abund_norm_wide)
#create the "OTU matrix" aka viral genome relative abundances matrix  
vmagmat <- as.matrix(abund_norm_wide, rownames = "Genome")

#make sure the column order is alphabetical
vmagmat <- vmagmat[sort(rownames(vmagmat)), sort(colnames(vmagmat))]

# write.table(vmagmat, file = "output/vmag_mat.tsv", quote = FALSE,
#             row.names = FALSE, sep = "\t")

######## taxonomy metadata matrix

missing_genomes <- abund_norm_wide$Genome[!(abund_norm_wide$Genome %in% genomad_tax$genome)]
new_rows <- data.frame(genome = missing_genomes,
                       r = "Unknown",
                       k = "Unknown",
                       p = "Unknown",
                       c = "Unknown",
                       o = "Unknown",
                       f = "Unknown",
                       g = "Unknown",
                       s = "Unknown")
genomad_tax2 <- rbind(genomad_tax, new_rows)

# find the difference in Genome names between genomad_tax2 and abund_norm_wide besides Lau_Basin_Tahi_Moana_vRhyme_bin_115 
unique_in_df1 <- setdiff(genomad_tax2$genome, abund_norm_wide$Genome)
print(unique_in_df1)
# these two viruses are absent in abund_norm_wide because they do not have coverage >70 so didn't pass the quality requirements

genomad_tax2 <- genomad_tax2 %>%
  filter(!(genome == "Brothers_NWCB_S146_NODE_255717_length_1218_cov_0.416132"|
           genome == "ELSC_Tui_Malila_T2_NODE_265519_length_1602_cov_0.296271"))

#convert to data_table and then matrix
genomad_tax2 <- as.data.table(genomad_tax2)
taxmat <- as.matrix(genomad_tax2, rownames = "genome")

#make sure the column order is alphabetical
taxmat <- taxmat[sort(rownames(taxmat)),]

######## sample metadata matrix

sampledata <- read.delim(file = "input/sample_metadata_div.txt", header = TRUE, row.names = "Site")
#make sure the column order is alphabetical
sampledata <- sampledata[sort(rownames(sampledata)),]

# #convert to data_table and then matrix
# sampledata <- as.data.table(sampledata, rownames = "Site")
#samplemat <- as.matrix(sampledata, rownames = "Site")

#################################################### phyloseq ######################################################

#loading phyloseq object code taken from https://joey711.github.io/phyloseq/import-data

#create OTU and TAX objects
OTU = otu_table(vmagmat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAMP = sample_data(sampledata)

#create phyloseq object
physeq = phyloseq(OTU, TAX, SAMP)

#check
otu_table(physeq)[1:5, 1:5]
ntaxa(physeq)
nsamples(physeq)
sample_variables(physeq)

# rest based on this guide: https://rpubs.com/lgschaerer/betadiv

#bray cutis ordination
#ordination
bray <- ordinate(
  physeq = physeq, #change this to your phyloseq
  method = "PCoA", 
  distance = "bray" 
)

head(sample_data(physeq))

###################################### phyloseq plot NMDS ordination ######################################################

dev.off()
plot_ordination(
  physeq = physeq,                                                          
  ordination = bray) + #,label = "Hydrothermal_Vent"
  geom_point(aes(color = General_Site, shape = Type), size = 3) + #shape size on NMDS
  scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1",
                            "#8c510a", "#000000")) +
  theme_linedraw() +
  geom_text(mapping = aes(label = Hydrothermal_Vent, hjust = 1, vjust = 1.5), size = 3) +
  theme(                             
    legend.title = element_blank(),                                          #removes legend title
    #legend.position = "bottom",
    legend.text = element_text(size = 15, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 15),
    strip.text = element_text(face = "bold", size = 15)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) #+          #fills legend points based on the fill command
  #facet_wrap(~Type)  
  
#have to separate plume and vent to see things


############################################ phyloseq PERMANOVA ######################################################

#Calculate distance and save as a matrix
bray <- phyloseq::distance(physeq, method = "bray")
#sam <- data.frame(sample_data(physeq))
#Run PERMANOVA on distances.
adonis2(bray ~ General_Site*Type*Latitude_DD*Longitude_DD*Depth_m, data = sampledata, permutations = 1000)



