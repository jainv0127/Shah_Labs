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

PNCognate <- function(codon){

  a <- strsplit(codon,",")
  neighbor <- list()
  for(val in c(1,3)){
    if (a[[val]] == "A"){
      a[[val]] <- "G"
      neighbor.append(a)
      a[[val]] <- "C"
      neighbor.append(a)
      a[[val]] <- "T"
      neighbor.append(a)
    }
    else if (a[[val]] == "G"){
      a[[val]] <- "A"
      neighbor.append(a)
      a[[val]] <- "C"
      neighbor.append(a)
      a[[val]] <- "T"
      neighbor.append(a)
    }
    else if (a[[val]] == "C"){
      a[[val]] <- "A"
      neighbor.append(a)
      a[[val]] <- "G"
      neighbor.append(a)
      a[[val]] <- "T"
      neighbor.append(a)
    }
    
  else if (a[[val]] == "T"){
    a[[val]] <- "A"
    neighbor.append(a)
    a[[val]] <- "G"
    neighbor.append(a)
    a[[val]] <- "C"
    neighbor.append(a)
  }
    for(val in c(0,lenght(a))){
      if (key[val, 2] != 0){
        if(a[[val,3]] == a[[codon,3]]){
          if(wobble(codon) != val){
            
            pcognate.append(val)
          }
          else{
            
            ncognate.append(val)
          }
        }
      
      
      
      }
    }
  
  
  }
}



checkWobble <- function(codon,trna){
  
  reverse <- reverse.complement(codon)
  trna.split <- unlist(strsplit(trna,""))
  codon.split <- unlist(strsplit(codon,""))
  reverse.split <- unlist(strsplit(reverse,""))
  if(reverse.complement(codon) == trna) 
  {
    return(1)
  } else{
    index <- which(reverse.split != trna.split)
    if(index == 1)
    {
      t.index <- 1
      c.index <- 3
    } else if (index == 2){
      t.index <- 2
      c.index <- 2
    } else if (index == 3){
      t.index <- 3
      c.index <- 1
    } else{
      print("Error")
    }
    if (trna.split[t.index] == "C" && codon.split[c.index] == "T")
    {
      return(canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "C"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "A"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "G"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "T"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "G"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "C"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "C" && codon.split[c.index] == "A"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "A"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "C" && codon.split[c.index] == "C"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "G"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "T"){
      return(canonical.wobble)
    }
    else{
      print("Nope")
    }
  }
}


EquationThree <- function(anticodon){
  w <- 0
  if (lenght(cognates) > 1){
    
  }
  input[anticodon,2] * (.652)
  cognates
  
}


