library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)

setwd("~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/iPHoP/")

######################################### Read the input ##################################################

#iphop Vent
iphop_vent <- read.csv(file = "input/Host_prediction_to_genus_m90.csv", header = TRUE)
#653 unique hosts without any filtering

#iphop results Plume
iphop_plume <- read.csv(file = "~/Google Drive/My Drive/Faith/PlumeViruses/AMG_homology/AMG_homology_network/input/Host_prediction_to_genus_m90_vMAGs_vUnbinned.csv")
#######combine the iphop results for both Vent and Plume
iphop_VentPlume <- rbind(iphop_vent, iphop_plume)

#CheckV
#checkv <- read.delim(file = "input/CheckV_quality_vMAGs_vUnbinned.tsv", header = TRUE, sep = "\t")

#Genome Size and CheckV quality for Plume and Vent viruses
sites_iphop <- read.delim(file = "../VentVirus_Analysis/output/gensize_VentPlume.tsv", header = TRUE, sep = "\t")

#GTDB MAG taxonomy for Vent
MAG_tax <- read.delim(file = "input/gtdbtk_v1.5.0_VentMAGs.tsv", header = TRUE)
#GTDB MAG tax for Plume
MAG_Tax_Plume <- read.delim(file = "~/Google Drive/My Drive/Faith/Plume_MAGs/GTDBtk_v1.5.0/all_PlumeMAGs_gtdbtkv1.5.0.txt", header = TRUE)
##########combine them into 1
mag_gtdb_VentPlume <- rbind(MAG_tax, MAG_Tax_Plume)

#virus taxonomy from genomad for Vent and Plume viruses
virus_tax <- read.delim(file = "../../genomad/vUnbinned_vMAGs_1500Ns_PlumeVent_genomad_tax_parsed.txt", header = TRUE)

################################## Map quality data so you can filter ########################################

#map
# iphop <- checkv %>%
#   dplyr::select(contig_id, contig_length, checkv_quality, provirus, gene_count, viral_genes,
#                 completeness, contamination, warnings) %>%
#   right_join(iphop, by = c("contig_id" = "Virus"))

#map
iphop_VentPlume <- sites_iphop %>%
  dplyr::select(contig_id, KB, checkv_quality, provirus,
                completeness, contamination, warnings) %>%
  right_join(iphop_VentPlume, by = c("contig_id" = "Virus"))

####################### Remove iphop results that don't match MAG taxonomy ########################################
#I will only keep predictions that match the MAG data that I have

#remove ;s_ in gtdbtk classification for mapping
mag_gtdb_VentPlume <- mag_gtdb_VentPlume %>% separate(classification, c("classification", NA), sep= "(?=;s__)")
#vlookup mapping
iphop_VentPlume <- mag_gtdb_VentPlume %>%
  dplyr::select(classification, user_genome) %>%
  right_join(iphop_VentPlume, by = c("classification" = "Host.genus"))
iphop_VentPlume <- virus_tax %>%
  dplyr::select(genome, lineage) %>%
  right_join(iphop_VentPlume, by = c("genome" = "contig_id"))
#drop NAs
iphop_VentPlume <- iphop_VentPlume %>%
  drop_na(user_genome)
#rename classification column
iphop_VentPlume<-rename(iphop_VentPlume,"Host.genus" = "classification")
iphop_VentPlume<-rename(iphop_VentPlume,"Virus" = "genome")
#258 unique hosts when filter by matching MAG taxonomy

################################### Quality control iphop results ########################################
#I am removing viruses whose checkv quality was Not determined because in my manual inspections,
#this gets rid of a lot of junk

iphop_VentPlume<-iphop_VentPlume[!grepl("Not-determined", iphop_VentPlume$checkv_quality),]

#Removing viruses with the warning "no viral genes detected" because my manual inspections suggest these are not viral
#or are poor enough quality that I don't want to keep

iphop_VentPlume<-iphop_VentPlume[!grepl("no viral genes detected", iphop_VentPlume$warnings),]

#Now removing viruses ≤5kb because I am not sure I trust host predictions to viral fragments
#And I'd like the potential for more genomic context from the virus

#filter for viruses with genome >5 KB
iphop_VentPlume<-subset(iphop_VentPlume, iphop_VentPlume$KB>=5)

#Remove viruses with contamination >20% because these don't look great

iphop_VentPlume<-subset(iphop_VentPlume, iphop_VentPlume$contamination<=20)
#194 unique hosts when filtering by all these quality metrics

##########################remove user_genome so I can see the table##########################
iphop_VentPlume<-select(iphop_VentPlume, -c("user_genome"))
iphop_VentPlume<-unique(iphop_VentPlume)

write.table(file = "output/iphop_VentPlumeresults_qc.txt", iphop_VentPlume, quote = FALSE, sep = "\t", col.names = TRUE,
            row.names = FALSE)

####################### Generate the network, starting with the nodes file ########################################

##subset virus and host genus
virus_host_iphop<-iphop_VentPlume[,c("Virus", "Host.genus", "List.of.methods")]

###USE THIS SEGMENT OF CODE IF YOU WANT EACH HOST MATCH TO BE UNIQUE (don't group by taxa)
##Add unique columns with x and y coordinates
###Add number to host.genus column so each one is unique
# virus_host_iphop$count <- 1:nrow(virus_host_iphop)
# virus_host_iphop$Host.genus <- paste0(virus_host_iphop$Host.genus, ".", 
#                                      virus_host_iphop$count)
# virus_host_iphop <- subset(virus_host_iphop, select = -c(count))


###USE THIS SEGMENT OF CODE IF YOU WANT NO NUMBER ASSOCIATED WITH HOST.GENUS TO MAKE THEM UNIQUE
###add columns where x and y coordinates will go, mapping nodes to each other
virus_host_iphop<-virus_host_iphop %>% mutate(x = NA, y = NA)

#add the site name to the microbial host name because you'll want it later for faceting
virus_host_iphop$Site <- virus_host_iphop$Virus # copy column
virus_host_iphop <- virus_host_iphop %>% separate(Site, c("Site", NA),
                                                  sep= "(?=_NODE|_scaffold|_vRhyme)") #separate by NODE and k95
virus_host_iphop$Host.genus<-paste0(virus_host_iphop$Site,
                                    "_", virus_host_iphop$Host.genus)

###generate unique list of number ids and map them back to virus_host_iphop so can see
###which viruses go with which MAGs
temp <- c(unique(virus_host_iphop$Virus), unique(virus_host_iphop$Host.genus))

###use temp to generate the nodes file
nodes <- data.frame(node = sort(temp), u_id = seq(1, length(temp)))

for (i in 1:nrow(virus_host_iphop)){
  bact.tmp <- virus_host_iphop[i, ]$Virus
  virus.tmp <- virus_host_iphop[i, ]$Host.genus
  
  virus_host_iphop$x[i] <- nodes[nodes$node == bact.tmp,2]
  virus_host_iphop$y[i] <- nodes[nodes$node == virus.tmp,2]
}

#check that input_network_uniqIDs is the length you expect (unique names of all viruses and MAGs)
#the sum of these values should be the num of observations in nodes
n_distinct(unique(virus_host_iphop$Virus))
n_distinct(unique(virus_host_iphop$Host.genus))

################################# Add metadata to the nodes

#Map Order (taxonomy) as metadata
nodes$node2 <- nodes$node # copy column
nodes <- nodes %>% separate(node2, c("d", "p", "c", "o", "f", "g"), 
                            sep= ";") #split by ;

nodes<-nodes[,c("node","u_id","p")] # CHANGE THIS FOR WHICH TAXONOMY LEVEL YOU WANT TO VISUALIZE
nodes$p <- nodes$p %>% replace_na('Virus') # AND CHANGE HERE FOR TAXONOMY

#Create category for virus or bac or arc
nodes$type <- nodes$p # copy column # AND CHANGE HERE FOR TAXONOMY LEVEL
nodes$type[grepl("p__", nodes$type, ignore.case=FALSE)] <- "Host" # AND CHANGE HERE FOR TAXONOMY LEVEL #if contains o__ replace with Host


#Create column of Site names and then change taxa site name to NA
sites_iphop <- nodes %>% separate(node, c("Site", NA),
                                  sep= "(?=_NODE|_scaffold|_d|_vRhyme|_k95)")
#replace specific site to get general sites
sites_iphop$Site <- gsub(".*Lau_Basin.*","Lau_Basin",sites_iphop$Site) #the placement of the periods is crucial for replacing whole string
sites_iphop$Site <- gsub(".*Cayman.*","Mid_Cayman_Rise",sites_iphop$Site)
sites_iphop$Site <- gsub(".*Axial.*","Axial_Seamount",sites_iphop$Site)
sites_iphop$Site <- gsub(".*ELSC.*","Lau_Basin",sites_iphop$Site)
sites_iphop$Site <- gsub(".*Brothers.*","Brothers_Volcano",sites_iphop$Site)
sites_iphop$Site <- gsub(".*Guaymas.*","Guaymas_Basin",sites_iphop$Site)
sites_iphop$Site <- gsub(".*MAR.*","Mid_Atlantic_Ridge",sites_iphop$Site)
sites_iphop$Site <- gsub(".*EPR.*","East_Pacific_Rise",sites_iphop$Site)

#sites_iphop$Site[grepl("d__", sites_iphop$Site, ignore.case=FALSE)] <- "NA"

#map it back to the nodes
nodes <- sites_iphop %>%
  dplyr::select(Site, u_id) %>%
  right_join(nodes, by = c("u_id" = "u_id")) %>%
  select(node, u_id, p, type, Site) # AND CHANGE HERE FOR TAXONOMY LEVEL

# Generate the edges file from virus_host_iphop

##First clean up the method column so you can keep it as metadata
virus_host_iphop$Method <- virus_host_iphop$List.of.methods
virus_host_iphop <- virus_host_iphop %>% separate(Method, c("Method", NA),
                                                  sep= ";")
#make the edges file
edges <- as.data.frame(virus_host_iphop[, c("x","y","Method")])

########################################## plot the network ###############################################

######graph
#mods to Site in nodes for simpler plotting
#nodes_SiteShort <- nodes %>% separate(Site, c("Site", NA), sep= "(?<=Brothers|ELSC|EPR|Guaymas|MAR)")
# make a tidy graph
network <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) 
network
# setting theme_graph 
set_graph_style()

######colors

#choose the colors you want
my_colors <- rev(pals::kelly(n=16)[2:17]) #change n = for the number of colors you want for nodes
#edge_colors <- pals::kelly(n=2)[c(2,7)]

library(randomcoloR)
n <- 42
palette <- rev(distinctColorPalette(n))

library(RColorBrewer)
n <- 16
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector<-col_vector[1:16]
col_vector<-c("#7FC97F", "#BEAED4", "#FDC086",
              "#FFFF99", "#386CB0", "#F0027F",
              "#a6cee3", "#fb8072", "#1B9E77",
              "#D95F02", "#7570B3", "#e9a3c9",
              "#66A61E", "#E6AB02", "#A6761D", "#666666")

#####network

# ggraph plot of network
dev.off()
plot<- ggraph(network, layout = 'kk') +  #kk is pretty #dh better for large networks. fr for small. --> play around with this layout. See here: https://www.data-imaginist.com/2017/ggraph-introduction-layouts/
  #geom_edge_arc(strength = 0.2) + #mapping = NULL to turn off #aes(colour = Method) to color by method
  geom_node_point(aes(shape = Site, color = p), size = 2.5, alpha = .9) + #Site is the color of the nodes, shape = for a shape
  #Aesthetic: change size = # for larger text within node when fewer nodes
  #geom_node_text(aes(label = u_id), size = 3, color = "white") + #I add the number label associated with u_id --> remove it if you have too many
  scale_color_manual("Phyla", values = rev(palette)) +
  scale_edge_color_discrete() + #change the first number in FALSE = depending on what percent identity cut off you used
  scale_shape_manual(name="Hydrothermal Vent Site", labels = c("NA" = ""),
                     values = c("Axial_Seamount"=15, "Brothers_Volcano"=17,
                                "Mid_Cayman_Rise"=18, "Lau_Basin"=25,
                                "East_Pacific_Rise"=7, "Guaymas_Basin"=9, "Lau_Basin_Tahi_Moana"=12,
                                "Mid_Atlantic_Ridge"=11)) +
  theme_bw() +
  # theme(plot.background = element_blank(),
  #       panel.background = element_blank(),
  #       panel.border = element_blank(),
  #       panel.grid = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank()) +
  #facet_nodes(~ Site, scales = "free") +
  th_foreground(border = TRUE) +
  xlab(NULL) +
  ylab(NULL)
plot

#export image of network --> change the name depending on what percent identity cut off you used
ggsave("output/iphop_network_vUnbinned_vMAGs_PhylaFacet_SiteColor.png", plot, width = 14,
       height = 10, units = "in")


########################## make some plots about the iphop output ################################

########################## make the input file from virus_host_iphop ################################
#use the virus_host_iphop file to make some plots and get some metadata

virus_host_iphop_plot<-virus_host_iphop
virus_host_iphop_plot$Site<-virus_host_iphop_plot$Host.genus

virus_host_iphop_plot <- virus_host_iphop_plot %>% separate(Site, c("Site", NA), 
                                                            sep= "_d")

virus_host_iphop_plot <- virus_host_iphop_plot %>% separate(Host.genus, c(NA, "Host.genus"), 
                                                            sep= "(?=d__)")

virus_host_iphop_plot <- virus_host_iphop_plot %>% separate(Host.genus, c("d", "p", "c", "o", "f", "g"), 
                                                            sep= ";") #split by ;

virus_host_iphop_plot_bac <- virus_host_iphop_plot %>% filter(grepl('d__Bacteria', d))
virus_host_iphop_plot_arc <- virus_host_iphop_plot %>% filter(grepl('d__Archaea', d))

########################## plot hosts by class ################################

for_plotting <- virus_host_iphop_plot_bac %>%
  select(Virus, Site, c, d) %>%
  group_by(d) %>%
  count(c) %>%
  arrange(c) %>%
  mutate(c = str_replace(c, "c__", ""))

for_plotting$c <- factor(for_plotting$c)

dev.off()
p <- ggplot(for_plotting, aes(x = c, y = n, fill = c)) + 
  geom_bar(stat = "identity") + 
  #facet_wrap(~Site) + 
  xlab("Class")  +
  ylab("Number of Predicted Hosts") +
  scale_fill_viridis_d(name = "GTDBtk class") +
  ggtitle("Bacteria") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400)) + #
  scale_x_discrete(limits = rev(levels(for_plotting$c)))
#geom_text(aes(label = paste0(n), y = n),
#         vjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(legend.position = "none") +
  coord_flip()
  #facet_wrap(~ d, scales = "free")
p

for_plotting <- virus_host_iphop_plot_arc %>%
  select(Virus, Site, c, d) %>%
  group_by(d) %>%
  count(c) %>%
  arrange(c) %>%
  mutate(c = str_replace(c, "c__", ""))

for_plotting$c <- factor(for_plotting$c)

dev.off()
p2 <- ggplot(for_plotting, aes(x = c, y = n, fill = c)) + 
  geom_bar(stat = "identity") + 
  #facet_wrap(~Site) + 
  xlab("")  +
  ylab("Number of Predicted Hosts") +
  scale_fill_viridis_d(name = "GTDBtk class") +
  ggtitle("Archaea") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 40)) + #, limits = c(0, 400)
  scale_x_discrete(limits = rev(levels(for_plotting$c)))
#geom_text(aes(label = paste0(n), y = n),
#         vjust = -.5, size = 2.5, color = "black" )
p2 <- p2 + theme_bw() + 
  theme(legend.position = "none") +
  coord_flip()
#facet_wrap(~ d, scales = "free")
p2

library(patchwork)
p_final <- p + p2
p_final

ggsave("output/iphop_hosts_barplot_Class_VentPlume.png", p_final, width = 13, height = 9, units = "in")

########################## zoom in on Campylo and Gamma ################################

#drop Pseudomonas contam
virus_host_iphop_plot <- virus_host_iphop_plot %>% filter(!str_detect(g, "g__Pseudomonas_D"))
#only keep Gamma and Camp
virus_host_iphop_plot_camp_gam_alph <- virus_host_iphop_plot_bac %>% filter(str_detect(c, "c__Gammaproteobacteria|c__Campylobacteria|c__Alphaproteobacteria"))
virus_host_iphop_plot_thermos <- virus_host_iphop_plot_arc %>% filter(str_detect(c, "c__Thermococci|c__Thermoproteia|c__Thermoplasmata"))

 
# for_plotting <- virus_host_iphop_plot_camp_gam %>%
#   select(Virus, Site, c, o) %>%
#   group_by(o) %>%
#   mutate(count=n()) %>%
#   arrange(o)

for_plotting <- virus_host_iphop_plot_camp_gam_alph %>%
  select(Virus, c, g) %>%
  group_by(c) %>%
  count(g) %>%
  arrange(g) 

for_plotting$g <- factor(for_plotting$g)

# count.df <- as.data.frame(table(virus_host_iphop_plot_camp_gam$o))
# colnames(count.df) <- c("o", "Count")

dev.off()
p <- ggplot(for_plotting, aes(x = g, y = n, fill = g)) + 
  geom_bar(stat = "identity") + 
  #facet_wrap(~Site) + 
  xlab("Genus")  +
  ylab("Number of Predicted Hosts") +
  scale_fill_viridis_d(name = "GTDBtk class") +
  ggtitle("Bacteria") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 250)) + #
  scale_x_discrete(limits = rev(levels(for_plotting$g)))
#geom_text(aes(label = paste0(n), y = n),
#         vjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 12)) +
  coord_flip()
  #facet_wrap(~ c)
p

for_plotting <- virus_host_iphop_plot_thermos %>%
  select(Virus, c, g) %>%
  group_by(c) %>%
  count(g) %>%
  arrange(g) 

for_plotting$g <- factor(for_plotting$g)

dev.off()
p <- ggplot(for_plotting, aes(x = g, y = n, fill = g)) + 
  geom_bar(stat = "identity") + 
  #facet_wrap(~Site) + 
  xlab("Genus")  +
  ylab("Number of Predicted Hosts") +
  scale_fill_viridis_d(name = "GTDBtk class") +
  ggtitle("Archaea") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 40)) + #
  scale_x_discrete(limits = rev(levels(for_plotting$g)))
#geom_text(aes(label = paste0(n), y = n),
#         vjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 12)) +
  coord_flip()
#facet_wrap(~ c)
p

ggsave("output/iphop_hosts_barplot_GammaCamp_VentPlume.png", p, width = 12, height = 9, units = "in")

########################## plot sites separately ################################











