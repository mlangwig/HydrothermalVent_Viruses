################################### Process abundance data for correlation ####################################

library(dplyr)

#read input
vir_abun <- read_delim(file = "input/PlumeVentVirus-vs-Reads-CoverM-Count_MinCov.tsv")
mic_abun <- read_delim(file = "input/PlumeVentMicrobialMAGs-vs-Reads-CoverM-Count_MinCov.tsv")

#remove unmapped and unnecessary columns
vir_abun <- vir_abun[-1,]
mic_abun <- mic_abun[-1,]

################################## Grabbing relative abundance ##################################

#grab relative abundance value
vir_abun <- vir_abun %>%
  select(matches("Genome|Relative"))

mic_abun <- mic_abun %>%
  select(matches("Genome|Relative"))

#add microbe taxonomy to sum abundance by it
mic_abun <- gtdb %>%
  select(c("user_genome", "c")) %>%
  right_join(mic_abun, by = c("user_genome" = "Genome"))
#filter for just Gamma and Camp
mic_abun <- mic_abun %>%
  filter(c %in% c("c__Gammaproteobacteria", "c__Campylobacteria"))

#group by and sum abundance
mic_abun_sum <- mic_abun %>%
  group_by(c) %>%
  summarize(across(contains("Relative"), sum, na.rm = TRUE))

#add virus host prediction to sum abundance by it
iphop_gc <- iphop %>%
  separate("Host.genus", c("d", "p", "c", "o", "f", "g"), sep= ";")
#filter for just Gamma and Camp infecting
iphop_gc <- iphop_gc %>%
  filter(c %in% c("c__Gammaproteobacteria", "c__Campylobacteria"))
#remove redundancy
iphop_gc <- iphop_gc %>%
  select(c("Virus", "c")) %>%
  unique()
#print the viruses that have multiple predictions
length(unique(iphop_gc$Virus))
#2567 tells me there are viruses that occur more than once
iphop_gc %>%
  count(Virus) %>%
  filter(n > 1) %>%
  pull(Virus) %>%
  print()
#this shows 2 viruses that have multiple predictions, Brothers_NWCB_S012_vRhyme_bin_152 and Brothers_NWCB_S141_vRhyme_bin_38
# bin 152 has CRISPR to Campylo so I'm keeping that one
# bin 38 has CRISPR to Gamma so I'm keeping that one
iphop_gc <- iphop_gc %>%
  filter(!(Virus == "Brothers_NWCB_S012_vRhyme_bin_152" & c == "c__Gammaproteobacteria")) %>%
  filter(!(Virus == "Brothers_NWCB_S141_vRhyme_bin_38" & c == "c__Campylobacteria"))
#check that viruses are now unique
length(unique(iphop_gc$Virus))
length(iphop_gc$Virus)

#map microbe the virus infects and only keep Gammas and Camps
vir_abun <- iphop_gc %>%
  inner_join(vir_abun, by = c("Virus" = "Genome"))

#erase unique virus ID
vir_abun2 <- vir_abun
vir_abun2$Virus <- "virus"
#join names
vir_abun2$virus_host <- paste(vir_abun2$Virus, vir_abun2$c, sep = "_")
# Move the merged column to be the first column
vir_abun2 <- vir_abun2 %>%
  select(virus_host, everything())

#group by and sum abundance
vir_abun_sum <- vir_abun2 %>%
  group_by(virus_host) %>%
  summarize(across(contains("Relative"), sum, na.rm = TRUE))

#remove first column name
colnames(vir_abun_sum)[1] <- ""
colnames(mic_abun_sum)[1] <- ""

#fix col names
#remove rel abun from colnames
colnames(vir_abun_sum) <- gsub(" Relative Abundance \\(\\%)", "", colnames(vir_abun_sum))
colnames(mic_abun_sum) <- gsub(" Relative Abundance \\(\\%)", "", colnames(mic_abun_sum))

#fix abun_names
read_mapping <- abun_names %>%
  select(c("Reads", "Site"))

#map
colnames(vir_abun_sum) <- ifelse(colnames(vir_abun_sum) %in% names(abun_names),
                        name_mapping[colnames(vir_abun_sum)], 
                        colnames(vir_abun_sum))

#write output
write_csv(vir_abun_sum, file = "output/class_VentVirus_abundance_GamCamp.csv")
write_csv(mic_abun_sum, file = "output/class_VentMicrobe_abundance_GamCamp.csv")

########################### Grabbing reads mapped normalized by # reads in sample ###########################


