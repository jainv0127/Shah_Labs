library(Biostrings)
library(tidyverse)
ecoli <- readDNAStringSet("CDS/Pseudoalteromonas_translucida/ptranslucida_tac125_cds_filtered.fasta")

ecoli.no.pseudo <- ecoli[which(!str_detect(names(ecoli),"pseudo=true"))]
ecoli.complete <- ecoli.no.pseudo[which(width(ecoli.no.pseudo) %% 3 == 0)]
ecoli.complete.long <- ecoli.complete[which(width(ecoli.complete)/3 > 50)]

aa <- translate(ecoli.complete.long)
nstopcodon <- vcountPattern("*", aa, max.mismatch=0)

ecoli.complete.long.no.stop <- ecoli.complete.long[which(nstopcodon == 1)]
writeXStringSet(ecoli.complete.long.no.stop,"CDS/Pseudoalteromonas_translucida/ptranslucida_tac125_cds_filtered.fasta")
