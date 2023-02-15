library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)

setwd("~/Google Drive/My Drive/PhD_Projects/VentViruses/HydrothermalVent_Viruses/iPHoP/")

######################################### Read the input ##################################################

#iphop
iphop <- read.csv(file = "input/Host_prediction_to_genus_m90.csv", header = TRUE)
#653 unique hosts without any filtering

#CheckV
checkv <- read.delim(file = "input/CheckV_quality_vMAGs_vUnbinned.tsv", header = TRUE, sep = "\t")

#Genome Size
gensize_kb <- read.delim(file = "input/Vent_vUnbinned_vMAG_GenSize_KB.tsv", header = TRUE, sep = "\t")

MAG_tax<-read.delim(file = "input/gtdbtk_v1.5.0_VentMAGs.tsv", header = TRUE)

################################## Map quality data so you can filter ########################################

#map
iphop <- checkv %>%
  dplyr::select(contig_id, contig_length, checkv_quality, provirus, gene_count, viral_genes,
                completeness, contamination, warnings) %>%
  right_join(iphop, by = c("contig_id" = "Virus"))

#map
iphop <- gensize_kb %>%
  dplyr::select(Genome, KB) %>%
  right_join(iphop, by = c("Genome" = "contig_id"))

####################### Remove iphop results that don't match MAG taxonomy ########################################
#I will only keep predictions that match the MAG data that I have

#remove ;s_ in gtdbtk classification for mapping
MAG_tax <- MAG_tax %>% separate(classification, c("classification", NA), sep= "(?=;s__)")
#vlookup mapping
iphop <- MAG_tax %>%
  dplyr::select(classification, user_genome) %>%
  right_join(iphop, by = c("classification" = "Host.genus"))
#drop NAs
iphop <- iphop %>%
  drop_na(user_genome)
#rename classification column
iphop<-rename(iphop,"Host.genus" = "classification")
iphop<-rename(iphop,"Virus" = "Genome")
#258 unique hosts when filter by matching MAG taxonomy

################################### Quality control iphop results ########################################
#I am removing viruses whose checkv quality was Not determined because in my manual inspections,
#this gets rid of a lot of junk

iphop<-iphop[!grepl("Not-determined", iphop$checkv_quality),]

#Removing viruses with the warning "no viral genes detected" because my manual inspections suggest these are not viral
#or are poor enough quality that I don't want to keep

iphop<-iphop[!grepl("no viral genes detected", iphop$warnings),]

#Now removing viruses ≤5kb because I am not sure I trust host predictions to viral fragments
#And I'd like the potential for more genomic context from the virus

#filter for viruses with genome >5 KB
iphop<-subset(iphop, iphop$KB>=5000)

#Remove viruses with contamination >20% because these don't look great

iphop<-subset(iphop, iphop$contamination<=20)
#194 unique hosts when filtering by all these quality metrics

##########################remove user_genome so I can see the table##########################
iphop<-select(iphop, -c("user_genome"))
iphop<-unique(iphop)

write.table(file = "output/iphop_results_qc.txt", iphop, quote = FALSE, sep = "\t", col.names = TRUE,
            row.names = FALSE)

####################### Generate the network, starting with the nodes file ########################################

##subset virus and host genus
virus_host_iphop<-iphop[,c("Virus", "Host.genus", "List.of.methods")]

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
                                  sep= "(?=_NODE|_scaffold|_d|_vRhyme)")
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
n <- 36
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
  geom_edge_link(mapping = NULL) + #mapping = NULL to turn off #aes(colour = Method) to color by method
  geom_node_point(aes(shape = type, color = p), size = 3, alpha = .9) + #Site is the color of the nodes
  #Aesthetic: change size = # for larger text within node when fewer nodes
  #geom_node_text(aes(label = u_id), size = 3, color = "white") + #I add the number label associated with u_id --> remove it if you have too many
  scale_color_manual("GTDBtk Phyla", values = rev(palette)) +
  #scale_edge_color_discrete() + #change the first number in FALSE = depending on what percent identity cut off you used
  # scale_shape_manual(name="Hydrothermal Vent Site", labels = c("NA" = ""),
  #                    values = c("Cayman_Deep"=15, "Cayman_Shallow"=17, 
  #                               "Guaymas_Basin"=18, "Lau_Basin_Abe"=25, 
  #                               "Lau_Basin_Kilo_Moana"=7, "Lau_Basin_Mariner"=9, "Lau_Basin_Tahi_Moana"=12, 
  #                               "Lau_Basin_Tui_Malila"=11, "SWIR"=8, "NA"=19)) + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_nodes(~ Site, scales = "free") +
  th_foreground(border = TRUE) +
  xlab(NULL) +
  ylab(NULL)
plot

#export image of network --> change the name depending on what percent identity cut off you used
ggsave("output/iphop_network_vUnbinned_vMAGs_phyla_facet.png", plot, width = 14,
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

########################## plot hosts by class ################################

for_plotting <- virus_host_iphop_plot %>%
  select(Virus, Site, c) %>%
  group_by(Site) %>%
  count(c) %>%
  arrange(c) 

for_plotting$c <- factor(for_plotting$c)

dev.off()
p <- ggplot(for_plotting, aes(x = c, y = n, fill = c)) + 
  geom_bar(stat = "identity") + 
  #facet_wrap(~Site) + 
  xlab("Taxonomy")  +
  ylab("Count") +
  scale_fill_viridis_d(name = "GTDBtk class") +
  ggtitle("iPHoP-Predicted Microbial Hosts") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400)) +
  scale_x_discrete(limits = rev(levels(for_plotting$c)))
#geom_text(aes(label = paste0(n), y = n),
#         vjust = -.5, size = 2.5, color = "black" )
p <- p + theme_bw() + 
  theme(legend.position = "none") +
  coord_flip()
p

ggsave("output/iphop_hosts_barplot_Class_nolegend.png", p, width = 10, height = 8, units = "in")














