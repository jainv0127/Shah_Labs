library(tidyverse)
library(ggplot2)
library(stringr)
#input <- read.csv(file = 'ecoli_tRNA_2010', header = TRUE, stringsAsFactors = TRUE)
input <- read_tsv(file = "ecoli_tRNA_2010.tsv")#, header = F, sep = ',', stringsAsFactors = F)
key <- read_tsv(file = 'shah_gilchrist_2010_ecoli.tsv')#, header = T, sep = ' ', stringsAsFactors = F)

#print(input)
#ggplot(data = key, aes(x = as.integer(key[,6]), y = as.integer(key[,8])) +
  #geom_point()
#print(input[x,])
for(val in c(1:nrow(input)
)){
  print(val)

  }
wobble <- function(codon){
  
  i <- substring(codon,3,3)
  if(i == "A"){
    substring(codon,3,3) <- "U"
  }
  if(i == "G"){
    substring(codon,3,3) <- "U"
  }
  if(i == "C"){
    substring(codon,3,3) <- "U"
  }
  if(i == "A"){
    substring(codon,3,3) <- "U"
  }
}
anticodonGenerator <- function(codon){

  string_split = strsplit(codon, split = "")
  rev_order = nchar(codon):1
  reversed_chars = string_split[[1]][rev_order]
  paste(reversed_chars, collapse = "")
  #print(reversed_chars)
  for(val in c(1:3)){
  if (reversed_chars[val] == "A"){
    reversed_chars[val] <- 'T'
  }
  else if (reversed_chars[val] == 'T'){
    reversed_chars[val] <- "A"
  }
  else if (reversed_chars[val] == "G"){
    reversed_chars[val] <- "C"
  }
  else if (reversed_chars[val] == "C"){
    reversed_chars[val] <- "G"
  }
  }
  return(reversed_chars)
}
  #codon <- input[val,1]
  #print(codon)
  #print(nrow(key))
  #for(val in c(1:nrow(key))){

    #k <- key[val,2]  
    #print(k)
    #if(k == codon){
      #print(key[val,])

    #}
  #}


#use join function(tidyverse)
  #give key, codon
#