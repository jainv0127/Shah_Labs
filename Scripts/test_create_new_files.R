library(tidyverse)
library(Biostrings)
bacteria <- read_tsv("bacteriatgcn_2023_03_29.tsv")
classification <- read_tsv("bacteriaClassification_2023_03_29.tsv")

classification <- classification %>% 
  dplyr::select(Organism_Name, Temperature_Range) %>% # this just gets the columns of interest
  filter(!is.na(Temperature_Range)) %>%
  mutate(Temperature_Range = ifelse(Temperature_Range == "Hyperthermophilic","Thermophilic",Temperature_Range))
  

bacteria.tgcn <- bacteria %>% 
  inner_join(classification, by="Organism_Name")

bacteria.tgcn.long <- bacteria.tgcn %>% 
  pivot_longer(c(-Organism_Name,-Temperature_Range), names_to="Codon",values_to = "tRNA") 


bacteria.tgcn.long <- bacteria.tgcn.long %>% 
  mutate(Codon = as.character(Biostrings::reverseComplement(DNAStringSet(Codon))),
          AA = as.character(translate(DNAStringSet(Codon),no.init.codon = T))
         ) %>% 
  filter(!Codon %in% c("TAG","TAA","TGA")) %>% # Remove stop codons 
  mutate(AA = ifelse(Codon %in% c("AGC","AGT"),"Z",AA),
         Organism_Name = str_replace_all(Organism_Name,"/","_"))

bacteria.tgcn.long.split <- bacteria.tgcn.long %>% 
  group_by(Organism_Name) %>% 
  group_split()

lapply(bacteria.tgcn.long.split, function(x) {
  species <- paste0(unname(unlist(x[1,"Organism_Name"])),".tsv")
  type <- unname(unlist(x[1,"Temperature_Range"]))
  if (type == "Mesophilic")
  {
    type <- "Mesophile"
  } else if (type == "Thermophilic") {
    type <- "Thermophile"
  } else if (type == "Psychrophilic") {
    type <- "Psychrophile"
  } else {
    print(species)
  }
  x <- x %>% 
    dplyr::select(Codon,tRNA,AA)
  write_tsv(x,file.path("tGCN_New",type,species))
})
