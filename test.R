#git add
#git commit
#git push

library(tidyverse)
library(ggplot2)
library(stringr)

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

## Here is a function you can use to get the reverse complement of a codon 
reverseComplement <- function(seq)
{
  seq.split <- unlist(strsplit(seq,split=""))
  reverse <- rev(seq.split)
  trna <- unlist(lapply(reverse,complementRules))
  trna <- paste(trna,collapse='')
  return(trna)
  
}


## I would work with T rather than U. I think most of our files use T.
wobble <- function(anticodon){  
  
  makeLists(codon)
  
  val <- 1
  
  if(substring(anticodon,val,val) == "A"){
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("T","C","G", "A")
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
    # print(cognates)
  }
  if(substring(anticodon,val,val) == "G"){
    
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("G","A") ## when you do the reverse complement, will give correct wobble rules between tRNA and codon
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
    # print(cognates)
    
  } 
  
  if(substring(anticodon,val,val) == "T"){
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("T","C","G", "A")
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
    #print(cognates)
    
  }
  
  if(substring(anticodon,val,val) == "C"){
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("C")
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

## Always make sure to pass in inputs as function arguments. Don't rely on global variables.
## Can get you in trouble. 
removeAnticodon <- function(input,codon)
{
  tgcn <- input %>% filter(Codon == codon) %>% dplyr::select(tRNA) %>% deframe()
  if(tgcn==0){
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
  canonical.wobble <- 0.6
  non.canonical.wobble <- 0.64
  reverse <- reverseComplement(codon)
  trna.split <- unlist(strsplit(trna,""))
  codon.split <- unlist(strsplit(codon,""))
  reverse.split <- unlist(strsplit(reverse,""))
  if(reverseComplement(codon) == trna) 
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

EquationOne <- function(codon){
  
  print(EquationFour(codon)/(EquationThree(codon)+EquationFour(codon)+0.0003146))
  
}

EquationTwo <- function(codon){
  #return not print
  #adijust equation 1,2
  #compare cognates and values
  print(0.0003146/(EquationThree(codon)+EquationFour(codon)+0.0003146))
  
}

EquationThree <- function(codon, cognates, psuedocogantes, input){
  rc <- 0
  
  for (val in length(cognates)){
    anticodon = reverseComplement(val)
    wj <- checkWobble(codon,val)
    rc = rc + input[anticodon,2] * (.652) * wj
    
      
  }
  
  for (val in length(psuedocognates)){
    anticodon = reverseComplement(val)
    wj <- checkWobble(codon,val)
    rc = rc + input[anticodon,2] * (.00062) * wj
  }
   rc = 10.992*rc
   return(rc)
}

EquationFour <- function(codon, nearcognates, input){
  rn <- 0
  
  for (val in length(nearcognates)){
    anticodon = reverseComplement(val)
    wj <- checkWobble(codon, val)
    rn = rn + input[anticodon,2] * (.00062) * wj
  }
  
  rn = 10.992*rn
  return(rn)
}

EquationOne("AGA")

makeLists <- function(codon){
  
  cognates <- list()
  ncognates <- list()
  pcognates <- list()
  
  
}


#input <- read.csv(file = 'ecoli_tRNA_2010', header = TRUE, stringsAsFactors = TRUE)
input <- read_tsv(file = "ecoli_tRNA_2010.tsv")#, header = F, sep = ',', stringsAsFactors = F)
key <- read_tsv(file = 'shah_gilchrist_2010_ecoli.tsv')#, header = T, sep = ' ', stringsAsFactors = F)
## the rownames() function with tibbles (dataframes read in from read_tsv) is deprectated, so let's remove that for now.

removeAnticodon(input,"AAG")

#print(input)
#ggplot(data = key, aes(x = as.integer(key[,6]), y = as.integer(key[,8])) +
#geom_point()
#print(input[x,])
for(val in 1:nrow(input))
{
  codon <- input[val,1]
  print(paste("Current codon:", codon$Codon))
  ## remember to reverse the coding sequence to get anticodon first
  anticodon <- reverseComplement(codon$Codon)
  print(paste("Current anticodon:",anticodon))
  wobbles <- unlist(lapply(wobble(anticodon),reverseComplement))
  #print(wobble(anticodonGenerator(codon)))
  print(paste0("Can wobble with codons:",paste(wobbles,collapse=",")))
}



wobble("AAA")


