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
library(ggrepel)
library(plyr)
set.seed(81)

######################################## see distribution of abundance data ######################################################
abund_long_norm_no <- abund_long_norm %>%
  filter(!(Genome == "Lau_Basin_Tahi_Moana_vRhyme_bin_115"))

dev.off()
p<-ggplot(abund_long_norm_no, aes(x=value)) + 
  geom_histogram(color="black", fill="white") #binwidth = .0001, 
p

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

#subset input to separate plume and deposit
physeq_p <- subset_samples(physeq, Type == "Plume")
physeq_d <- subset_samples(physeq, Type == "Deposit")

# rest based on this guide: https://rpubs.com/lgschaerer/betadiv

#bray cutis ordination
#ordination
bray <- ordinate(
  physeq = physeq, #change this to your phyloseq
  method = "NMDS", #PCoA
  k = 3,
  distance = "bray" 
)

head(sample_data(physeq))

###################################### phyloseq plot NMDS ordination ######################################################

dev.off()
plot_ordination(
  physeq = physeq,                                                          
  ordination = bray) + #,label = "Hydrothermal_Vent" #bray
  geom_point(aes(color = General_Site, shape = Type), size = 3) + #shape size on NMDS #, shape = Type
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
  #stat_ellipse(type = "norm", linetype = 2) +
  #stat_ellipse(type = "t")
  #facet_wrap(~Type)


############################################ phyloseq PERMANOVA ######################################################

#Calculate distance and save as a matrix
bray <- phyloseq::distance(physeq, method = "bray")
#sam <- data.frame(sample_data(physeq))
#Run PERMANOVA on distances.
adonis2(bray ~ General_Site*Type*Latitude_DD*Longitude_DD*Depth_m, data = sampledata, permutations = 1000)

#################################### phyloseq loop through all distance metrics ######################################################

#get list of the dist methods
dist_methods <- unlist(distanceMethodList)
#remove the ones that require a phylogenetic tree
dist_methods <- dist_methods[-(1:3)]
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]

#loop through each distance method and save to plot list
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(physeq_p, method=i)
  # Calculate ordination
  iMDS  <- ordinate(physeq_p, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeq_p, iMDS, color="General_Site", shape="Type")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
#forgot, need to remove method morisita, chao, and cao because non-integer data
#and maybe raup because empty species?

#combine and plot
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=General_Site, shape=Type))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on distinct distance metrics for Vent Virus dataset")
p











