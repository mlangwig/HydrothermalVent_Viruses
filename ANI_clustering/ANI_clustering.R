####### Determine the proportion of Plume:Vent viruses in clusters #######
library(tidyr)
library(dplyr)
library(stringr)
library(pals)
library(ggraph)
library(tidygraph)
library(reshape2)

#### Read input
dRep_clusters <- read.csv2(file = "Input/Cdb.csv", sep = ",")

mcl_clusters <- read.csv2(file = "Input/vUnbinned_vMAGs_PlumeVents.clusters.csv", sep = ",", header = FALSE)

############################################ dRep ####################################################

#### Transforming names so can get proportions
#separate names by NODE, scaffold, and vRhyme and remove that part of the name
dRep_clusters <- dRep_clusters %>% separate(genome, c("Site", NA), sep= "(?=_NODE|_scaffold|__|_k95)")
unique(dRep_clusters$Site)
#gsub to replace Plume names
dRep_clusters$Site <- gsub(".*Lau_Basin.*","Plume",dRep_clusters$Site) #the placement of the periods is crucial for replacing whole string
dRep_clusters$Site <- gsub(".*Cayman.*","Plume",dRep_clusters$Site)
dRep_clusters$Site <- gsub(".*Guaymas_Basin.*","Plume",dRep_clusters$Site)
dRep_clusters$Site <- gsub(".*Axial.*","Plume",dRep_clusters$Site)
#separate the plume and vent dataframes to make life easier
dRep_clusters_Plume <- dRep_clusters %>% filter(Site == "Plume")
#remove Plume from dRep so can make all names left Vent
dRep_clusters <- dRep_clusters %>% filter(Site != "Plume")
dRep_clusters$Site <- "Vent"
#Now put them back together
dRep_clusters <- rbind(dRep_clusters, dRep_clusters_Plume)

#Add counts to filter for clusters that occur more than once (not a singleton)
dRep_clusters_noSingle <- dRep_clusters %>%
  add_count(secondary_cluster) %>% 
  filter(n!=1) %>%
  select(-n)
#number of clusters with more than 1 rep according to dRep
length(unique(dRep_clusters_noSingle$secondary_cluster))
#641

#sort by unique
dRep_clusters_noSingle <- unique(dRep_clusters_noSingle)
#group by secondary cluster, count occurrences of Site
dRep_count <- dRep_clusters_noSingle %>% group_by(secondary_cluster) %>% count(Site)
#see if any cluster now occurs twice
dRep_count <-  dRep_count %>% group_by(secondary_cluster) %>% filter(n()>1)


############################################ mcl ####################################################

mcl_clusters <- read.csv2(file = "Input/vUnbinned_vMAGs_PlumeVents.clusters.csv", sep = ",", header = FALSE)

#add id number to rows
mcl_clusters <- mcl_clusters %>% mutate(id = row_number())
#melt by id
mcl_clusters <- melt(mcl_clusters, id.vars = "id")
#drop variable column
mcl_clusters <- select(mcl_clusters, -variable)
#drop blanks in value column
mcl_clusters <- mcl_clusters %>%
  na_if("") %>%
  na.omit()
#rename column
mcl_clusters <- rename(mcl_clusters, "Site" = "value")
#replace strings with Vent or Plume
mcl_clusters$Site <- gsub(".*Lau_Basin.*","Plume",mcl_clusters$Site) #the placement of the periods is crucial for replacing whole string
mcl_clusters$Site <- gsub(".*Cayman.*","Plume",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Guaymas_Basin.*","Plume",mcl_clusters$Site)
mcl_clusters$Site <- gsub(".*Axial.*","Plume",mcl_clusters$Site)
#separate the plume and vent dataframes to make life easier
mcl_clusters_Plume <- mcl_clusters %>% filter(Site == "Plume")
#remove Plume from mcl clusters so can make all names left Vent
mcl_clusters <- mcl_clusters %>% filter(Site != "Plume")
mcl_clusters$Site <- "Vent"
#Now put them back together
mcl_clusters <- rbind(mcl_clusters, mcl_clusters_Plume)

#Add counts to filter for clusters that occur more than once (not a singleton)
mcl_clusters_noSingle <- mcl_clusters %>%
  add_count(id) %>% 
  filter(n!=1) %>%
  select(-n)
#number of clusters with more than 1 rep according to dRep
length(unique(mcl_clusters_noSingle$id))
#2797

#group by secondary cluster, count occurrences of Site
mcl_count <- mcl_clusters_noSingle %>% group_by(id) %>% count(Site)
#see if any cluster now occurs twice
mcl_count <-  mcl_count %>% group_by(id) %>% filter(n()>1)


