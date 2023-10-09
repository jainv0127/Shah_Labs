library(tidyverse)
library(ape)
trnagcn <- read_tsv("20161122_int_file_trnagcn_complete.tsv")
temp <- read_csv("bacteria_ncbi_temperatures (1).csv")
trnagcn <- mutate(trnagcn, ncbi_acc = str_remove(ncbi_acc,"\\.[0-9]"))
temp <- mutate(temp, Assembly = str_remove(Assembly,"\\.[0-9]"))

trna.temp <- inner_join(x = temp,y = trnagcn, by = c("Assembly" = "ncbi_acc"))
trna.temp <- filter(trna.temp, !is.na(Temperature))


write_tsv(trna.temp, "Phylogeny/targetBacteriaForLargerAnyl.tsv")

tree <-  read.tree("Phylogeny/phyloTree.dated.smooth.1")

tipNames <- trna.temp$phylotree_tip_label
todrop <- tree$tip.label[which(!tree$tip.label %in% tipNames)]
tree.filt <- drop.tip(tree,todrop)


trna.temp.filt <- trna.temp %>% 
  dplyr::slice(match(tree.filt$tip.label,phylotree_tip_label))
tree.filt$tip.label <- trna.temp.filt$`#Organism Name`
tree.filt$tip.label <- str_replace_all(tree.filt$tip.label,"/","_")
write.tree(tree.filt, "Phylogeny/DatedAgashiPhyloFiltered.nwk")


trna.temp <- mutate(trna.temp, Category = case_when(
  Temperature < 15 ~ "Psychrophile",
  15 <= Temperature & Temperature <= 45 ~ "Mesophile",
  Temperature > 45 ~ "Thermophile"
))

temp.2 <- select(trna.temp, `#Organism Name`, Category)
tgcn.2 <- select(trna.temp, `#Organism Name`, matches("^[AGCT]{3}$"))
temp.2 <- dplyr::rename(temp.2, Organism_Name = `#Organism Name`, Temperature_Range = `Category`)
tgcn.2 <- dplyr::rename(tgcn.2, Organism_Name = `#Organism Name`)

write_tsv(temp.2, "bacteriaClassification_2023_03_29.tsv")
write_tsv(tgcn.2, "bacteriatgcn_2023_03_29.tsv")
