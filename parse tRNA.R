library(tidyverse)
library(seqinr)

complementRules <- function(nucleotide)
{
  if (nucleotide == "A")
  {
    return("T")
  } else if (nucleotide == "C"){
    return("G")
  } else if(nucleotide == "G"){
    return("C")
  } else if (nucleotide == "T"){
    return("A")
  }
}

reverseComplement <- function(seq)
{
  seq.split <- unlist(strsplit(seq,split=""))
  reverse <- rev(seq.split)
  trna <- unlist(lapply(reverse,complementRules))
  trna <- paste(trna,collapse='')
  return(trna)
}
#parsefile <- function(f){
df <- read_tsv("C:/Users/vtslj/Downloads/baciSubt2-tRNAs (1)/baciSubt2-tRNAs.out")
ecoli <- read_tsv("C:/Users/vtslj/OneDrive/Documents/scatterplot/project/tGCN/ecoli_tRNA_2010.tsv")
codon.df <- ecoli %>% dplyr::select(AA,Codon)

tgcn <- df %>% dplyr::count(Anti, tRNA_1)

tgcn <- tgcn %>% rowwise() %>% mutate(Codon = reverseComplement(Anti),
                              Codon = ifelse(tRNA_1=="Ile2","ATA",Codon)) %>%
  filter(tRNA_1 != "fMet")

tgcn = codon.df %>% left_join(tgcn,by="Codon")

tgcn <- tgcn %>% mutate(n = ifelse(is.na(n),0, n)) %>% dplyr::select(AA,Codon,n) %>% rename(tRNA = n) %>% select(Codon, tRNA, AA)

write_tsv(tgcn,"baciSubt2_tRNA.tsv")
#str_extract("../this/is/a/test/species-tRNAs.out","(?<=/)[A-Z0-9a-z]+-tRNAs(?=\\.out)")

#put all in function, list of files u want, for loop -> run through file list and use function
#list.files("directory to out files",full.names = TRUE)
#naming files