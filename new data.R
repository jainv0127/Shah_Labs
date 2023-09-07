library(tidyverse)
library(ape)
trnagcn <- read_tsv("20161122_int_file_trnagcn_complete.tsv")
temp <- read_csv("bacteria_ncbi_temperatures (1).csv")
trnagcn <- mutate(trnagcn, ncbi_acc = str_remove(ncbi_acc,"\\.[0-9]"))
temp <- mutate(temp, Assembly = str_remove(Assembly,"\\.[0-9]"))

a <- inner_join(x = temp,y = trnagcn, by = c("Assembly" = "ncbi_acc"))
a <- filter(a, !is.na(Temperature))


write_tsv(a, "Phylogeny/targetBacteriaForLargerAnyl.tsv")

tree = read.tree("phyloTree.dated.smooth.1")
plot(tree)

tipNames <- a$phylotree_tip_label
todrop <- tree$tip.label[which(!tree$tip.label %in% tipNames)]
tree.filt <- drop.tip(tree,todrop)
plot(tree.filt)

a.filt <- a %>% dplyr::slice(match(tree.filt$tip.label,phylotree_tip_label))
tree.filt$tip.label <- a.filt$`#Organism Name`
tree.filt$tip.label <- str_replace_all(tree.filt$tip.label,"/","_")
write.tree(tree.filt, "Phylogeny/DatedAgashiPhyloFiltered.nwk")
#turn tree^ in terms of time


b <- mutate(a, Category = case_when(
  Temperature < 15 ~ "Psychrophile",
  15 <= Temperature & Temperature <= 45 ~ "Mesophile",
  Temperature > 45 ~ "Thermophile"
))

btemp <- select(b, `#Organism Name`, Category)
btgcn <- select(b, `#Organism Name`, matches("^[AGCT]{3}$"))
btemp <- dplyr::rename(btemp, Organism_Name = `#Organism Name`, Temperature_Range = `Category`)
btgcn <- dplyr::rename(btgcn, Organism_Name = `#Organism Name`)

write_tsv(btemp, "bacteriaClassification_2023_03_29.tsv")
write_tsv(btgcn, "bacteriatgcn_2023_03_29.tsv")
