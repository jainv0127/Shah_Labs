library(tidyverse)
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


getNeighbors <- function(codon)
{
  neighbors <- character(length=12)
  nuc <- unlist(strsplit(codon,split='',fixed=T))
  neighbors[1:4] <- paste0(nuc[1],c("A","G","T","C"),nuc[3])
  neighbors[5:8] <- paste0(c("A","G","T","C"),nuc[2],nuc[3])
  neighbors[9:12] <- paste0(nuc[1],nuc[2],c("A","G","T","C"))
  return(neighbors)
}

#returns wobble
getWobble <- function(anticodon,aa){  
  
  val <- 1
  if(substring(anticodon,val,val) == "A"){
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("T","C","G", "A")
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
  }
  if(substring(anticodon,val,val) == "G"){
    
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("G","A") ## when you do the reverse complement, will give correct wobble rules between tRNA and codon
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
  } 
  
  if(substring(anticodon,val,val) == "T"){
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    wobble <- c("T","C","G", "A")
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
  }
  
  if(substring(anticodon,val,val) == "C"){
    unlist(strsplit(anticodon,""))
    anticodon.split <- unlist(strsplit(anticodon,""))
    if (aa == "I")
    {
      wobble <- c("T")
    } else {
      wobble <- c("C")
    }
    
    a <- paste0(wobble, anticodon.split[2], anticodon.split[3])
    cognates <- c(a)
  }
  
  return(a)
}

#used to remove anticodon
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


check.same.codons<- function(set.1,set.2){
  
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


EquationThree <- function(codon, cognates, pseudo.cognates, input, aa){
  rc <- 0
  cognates <- unlist(str_split(cognates, ","))
  pseudo.cognates <- unlist(str_split(pseudo.cognates, ","))
  cognates = cognates[which(cognates != "")]
  pseudo.cognates = pseudo.cognates[which(pseudo.cognates != "")]
  if(length(cognates) != 0){
    for (val in cognates){
      tRNA <- input %>% filter(val == AntiCodon & aa == AA) %>% select(tRNA) %>% deframe()
      anticodon = reverseComplement(val)
      wj <- checkWobble(codon,val) 
      rc = rc + tRNA * (.652) * wj
    }
    
  }
  
  if(length(pseudo.cognates) != 0){
    
    for (val in pseudo.cognates){
      tRNA <- input %>% filter(val == AntiCodon & aa == AA) %>% select(tRNA) %>% deframe()
      anticodon = reverseComplement(val)
      wj <- checkWobble(codon,val)
      rc = rc + tRNA * (.00062) * wj
    }
  }
  
  rc = 10.992*rc
  return(rc)
  
}


EquationFour <- function(codon, nearcognates, input, aa){
  rn <- 0
  
  nearcognates <- unlist(str_split(nearcognates, ","))
  nearcognates = nearcognates[which(nearcognates != "")]
  
  if(length(nearcognates) != 0){
    for (val in nearcognates){
      tRNA <- input %>% filter(val == AntiCodon & aa != AA) %>% select(tRNA) %>% deframe()
      anticodon = reverseComplement(val)
      wj <- checkWobble(codon, val)
      for(t in tRNA){
        rn = rn + t * (.00062) * wj
        
      }
    }
  }
  rn = 10.992*rn
  return(rn)
}

calculateElongationMetrics <- function(input.file,target.dir="tGCN_new/",output.dir="Rates_Updated/")
{
  input <- read_tsv(file = file.path(target.dir,input.file))
  if (!dir.exists(output.dir)){
    dir.create(output.dir)
  }
  
  input <- input %>% rowwise() %>% 
    mutate(AntiCodon = reverseComplement(Codon)) %>% 
    mutate(AntiCodon=ifelse(Codon == "ATA","CAT",AntiCodon)) ## This is because of some weirdness with bacteria. 
  
  df <- data.frame(Codon = input$Codon,AA=input$AA,Cognate=NA,Pseudo.cognate=NA,Near.cognate=NA)
  rownames(df) <- df$Codon
  
  for(val in 1:nrow(input))
  {
    codon <- input[val,"Codon"] %>% deframe()
    codon.aa <- input[val,"AA"] %>% deframe()
    anticodon <- input[val,"AntiCodon"] %>% deframe()
    reverse.Comp <- reverseComplement(codon)
    one.step.neighbors <- getNeighbors(reverse.Comp)
    neighbors.df <- input %>% 
      filter((AntiCodon %in% one.step.neighbors))
    anticodon.neighbors <- neighbors.df$AntiCodon
    names(anticodon.neighbors) <- anticodon.neighbors
    
    
    wobbles <- purrr::map(anticodon.neighbors,function(x)
    {
      anticodon.aa <- neighbors.df %>% filter(AntiCodon == x) %>% 
        dplyr::select(AA) %>% deframe
      tmp <- data.frame(AA=anticodon.aa,AntiCodon=x,Wobble=rep(FALSE,length(anticodon.aa)))
      for (aa in anticodon.aa)
      {
        wobble.list <- unlist(lapply(getWobble(x,aa),reverseComplement))
        if (codon %in% wobble.list) 
        {
          tmp[which(tmp$AA == aa),"Wobble"] <- TRUE
        } 
      }
      tmp
    }) %>% bind_rows() %>% distinct()
    neighbors.df <- neighbors.df %>% left_join(wobbles,by=c("AntiCodon","AA"))
    
    neighbors.df <- neighbors.df %>% mutate(Category = case_when(
      AA == codon.aa & tRNA > 0 & Wobble ~"Cognate",
      AA == codon.aa & tRNA > 0 ~ "Pseudo-Cognate",
      AA != codon.aa & tRNA > 0 ~ "Near-Cognate"
    ) 
    )
    cognates <- neighbors.df %>% filter(Category == "Cognate") %>% dplyr::select(AntiCodon) %>% deframe() %>% unique()
    if (length(cognates) > 1)
    {
      cognates <- paste(cognates,collapse=",")
    } else if (length(cognates) == 0) {
      cognates <- "" 
    }
    
    
    pseudo.cognates <- neighbors.df %>% filter(Category == "Pseudo-Cognate") %>% dplyr::select(AntiCodon) %>% deframe()
    if (length(pseudo.cognates) > 1)
    {
      pseudo.cognates <- paste(pseudo.cognates,collapse=",")
    } else if (length(pseudo.cognates) == 0) {
      pseudo.cognates <- "" 
    }
    
    near.cognates <- neighbors.df %>% filter(Category == "Near-Cognate") %>% dplyr::select(AntiCodon) %>% deframe() %>% unique()
    if (length(near.cognates) > 1)
    {
      near.cognates <- paste(near.cognates,collapse=",")
    } else if (length(near.cognates) == 0) {
      near.cognates <- "" 
    }
    df[codon,"Cognate"] <- cognates
    df[codon,"Pseudo.cognate"] <- pseudo.cognates
    df[codon,"Near.cognate"] <- near.cognates
  }
  
  df <- df %>% arrange(AA,Codon)
  
  merge.df <- df
  
  
  for(i in 1:nrow(merge.df))
  {
    row = merge.df[i,]
    merge.df[i, "Rc"] <- EquationThree(row$Codon, row$Cognate, row$Pseudo.cognate, input, row$AA)
    merge.df[i, "Rn"] <- EquationFour(row$Codon, row$Near.cognate, input, row$AA)
  }
  merge.df['eN'] <- 0.003146/(merge.df$Rc + merge.df$Rn +0.003146)
  merge.df['eM'] <- merge.df$Rn/(merge.df$Rc + merge.df$Rn +0.003146)
  
  write_tsv(merge.df, file.path(output.dir,input.file))
}


tgcn.files <- list.files("../tGCN_new/Thermophile/",recursive=F,full.names = F)
rates <- lapply(tgcn.files,calculateElongationMetrics, output.dir = "../Rates_Updated/Thermophile", target.dir = "../tGCN_new/Thermophile/")

tgcn.files <- list.files("../tGCN_new/Psychrophile/",recursive=F,full.names = F)
rates <- lapply(tgcn.files,calculateElongationMetrics, output.dir = "../Rates_Updated/Psychrophile", target.dir = "../tGCN_new/Psychrophile/")

tgcn.files <- list.files("../tGCN_new/Mesophile/",recursive=F,full.names = F)
rates <- lapply(tgcn.files,calculateElongationMetrics, output.dir = "../Rates_Updated/Mesophile", target.dir = "../tGCN_new/Mesophile/")



  