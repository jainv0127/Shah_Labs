#git add
#git commit
#git push

library(tidyverse)
library(ggplot2)
library(stringr)
#input <- read.csv(file = 'ecoli_tRNA_2010', header = TRUE, stringsAsFactors = TRUE)
input <- read_tsv(file = "ecoli_tRNA_2010.tsv")#, header = F, sep = ',', stringsAsFactors = F)
key <- read_tsv(file = 'shah_gilchrist_2010_ecoli.tsv')#, header = T, sep = ' ', stringsAsFactors = F)
rownames(input) <- input$Codon
removeAnticodon("AAG")

#print(input)
#ggplot(data = key, aes(x = as.integer(key[,6]), y = as.integer(key[,8])) +
  #geom_point()
#print(input[x,])
for(val in c(1:nrow(input)
)){
  codon <- input[val,1]
  print(codon)
  print(wobble(codon$Codon))
  #print(wobble(anticodonGenerator(codon)))
  
}

makeLists <- function(codon){
  
  cognates <- list()
  ncognates <- list()
  pcognates <- list()
  
  
}

wobble("AAA")

wobble <- function(anticodon){  
  
  makeLists(codon)
  
  val <- 1

  print(anticodon)
    if(substring(anticodon,val,val) == "A"){
      unlist(strsplit(anticodon,""))
      anticodon.split <- unlist(strsplit(anticodon,""))
      wobble <- c("U","C","G", "A")
      a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
      cognates <- c(a)
     # print(cognates)
    }
    if(substring(anticodon,val,val) == "G"){
      
      unlist(strsplit(anticodon,""))
      anticodon.split <- unlist(strsplit(anticodon,""))
      wobble <- c("U","C")
      a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
      cognates <- c(a)
     # print(cognates)
      
      } 
    
    if(substring(anticodon,val,val) == "U"){
      unlist(strsplit(anticodon,""))
      anticodon.split <- unlist(strsplit(anticodon,""))
      wobble <- c("U","C","G", "A")
      a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
      cognates <- c(a)
     #print(cognates)
      
      }
  return(a)
}

#wobble rules function
  #parameter nucleotide
  #call wobble wobble function
  #paste specific spot


anticodonGenerator <- function(codon){

  string_split = strsplit(codon, split = "")
  rev_order = nchar(codon):1
  reversed_chars = string_split[[1]][rev_order]
  paste(reversed_chars, collapse = "")
  #print(reversed_chars)
  for(val in c(1:3)){
  if (reversed_chars[val] == "A"){
    reversed_chars[val] <- 'T'
    s<- paste(reversed_chars)
  }
  else if (reversed_chars[val] == 'T'){
    reversed_chars[val] <- "A"
    s<- paste(reversed_chars)
  }
  else if (reversed_chars[val] == "G"){
    reversed_chars[val] <- "C"
    s<- paste(reversed_chars)
  }
  else if (reversed_chars[val] == "C"){
    reversed_chars[val] <- "G"
    s<- paste(reversed_chars)
  }
  }
  return(s)
}
  

removeAnticodon <- function(codon){
  
  if(input[codon, tRNA]==0){
    print("nope")
  }
  else{
    print("valid")
  }
}

check.same.codons<- function(set.1,set.2)
{
  set.1 <- unlist(strsplit(set.1,","))
  set.2 <- unlist(strsplit(set.2,","))
  if (length(set.1) != length(set.2))
  {
    return(F)
  } else {
    disagree <- setdiff(set.1,set.2)
    if (length(disagree) >0)
    {
      return(F)
    }
    disagree <- setdiff(set.2,set.1)
    if (length(disagree) >0)
    {
      return(F)
    }
  }
  return(T)
}

pseudoCognate <- function(anticodon){

  
}
