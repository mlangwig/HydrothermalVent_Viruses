################################### Process abundance data for correlation ####################################

#read input
vir_abun <- read_delim(file = "input/PlumeVentVirus-vs-Reads-CoverM-Count_MinCov.tsv")
mic_abun <- read_delim(file = "input/PlumeVentMicrobialMAGs-vs-Reads-CoverM-Count_MinCov.tsv")

#remove unmapped and unnecessary columns