library(tidyverse)
library(Biostrings)
bacteria <- read_tsv("tGCN bacteria list.txt",skip = 1)
classification <- read_tsv("classification.txt",skip=1)

bacteria <- bacteria %>% 
  dplyr::rename(Organism_Name = `Organism Name`) # The ` is because of the space in the column name. We'll get rid of them. 

classification <- classification %>% 
  dplyr::rename(Organism_Name = `Organism Name`, Temperature_Range = `Temperature Range`) %>%
  dplyr::select(Organism_Name, Temperature_Range) %>% # this just gets the columns of interest
  filter(!is.na(Temperature_Range)) %>%
  mutate(Temperature_Range = ifelse(Temperature_Range == "Hyperthermophilic","Thermophilic",Temperature_Range))
  

bacteria.tgcn <- bacteria %>% 
  inner_join(classification, by="Organism_Name")

bacteria.tgcn.long <- bacteria.tgcn %>% 
  pivot_longer(c(-Organism_Name,-Temperature_Range), names_to="Codon",values_to = "tRNA") 


bacteria.tgcn.long <- bacteria.tgcn.long %>% 
  mutate(Codon = as.character(reverseComplement(DNAStringSet(Codon))),
          AA = as.character(translate(DNAStringSet(Codon),no.init.codon = T))
         ) %>% 
  filter(!Codon %in% c("TAG","TAA","TGA")) %>% # Remove stop codons 
  mutate(AA = ifelse(Codon %in% c("AGC","AGT"),"Z",AA))

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
  write_tsv(x,file.path("tGCN",type,species))
})
 