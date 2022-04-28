#git add
#git commit
#git push


library(tidyverse)
#library(ggplot2)
library(stringr)
AAColumn <- c()
CodonColumn <- c()
CognateColumn <- c()
PCColumn <- c()
NearCog <- c()
RC <- c()
RN <- c()
EM <- c()
EN <- c()

#move further down.

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


## I would work with T rather than U. I think most of our files use T.
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

#wobble rules function
#parameter nucleotide
#call wobble wobble function
#paste specific spot
anticodonGenerator <- function(codon){
  CodonColumn <- c(CodonColumn, codon)
  
  string_split = strsplit(codon, split = "")
  rev_order = nchar(codon):1
  reversed_chars = string_split[[1]][rev_order]
  paste(reversed_chars, collapse = "")
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
  CognateColumn <- c(CognateColumn, s)
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

#use function below on each row and column
#make sure each row agrees with key
check.same.codons<- function(set.1,set.2)
{
  #one set is my own set of psuedo/near cognate tRNA for given codon
  #other set is key
  #first check psudeo, near cognates
  #then check values
  #get value in data frame. merge data frames
  #plot data
  #x value key elongation rate(Rc)
  #y value estimate elongation rate
  
  #do the same for all other values
  
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

#merge on codon and AA

PNCognate <- function(codon,trna){
  
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
    for(val in seq(1,length(a)))
    {
      if (input[val, 2] != 0)
      {
        if(a[[val,3]] == a[[codon,3]])
        {
          if(getWobble(codon) != val)
          {
            pcognate.append(val)
          } else{
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
    #print(reverse.split)
    #print(trna.split)
    #print(which(reverse.split != trna.split))
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
  
  return(EquationFour(codon)/(EquationThree(codon)+EquationFour(codon)+0.0003146))
  
}

EquationTwo <- function(codon){
  #return not print
  #adijust equation 1,2
  p <- EquationThree(codon)
  q <- EquationFour(codon)
  #compare cognates and values
  return(0.0003146/(p + q +0.0003146))
  print('hi')
}

EquationThree <- function(codon, cognates, pseudo.cognates, input, aa){
  rc <- 0
  #splits up string
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


makeLists <- function(codon){
  cognates <- list()
  ncognates <- list()
  pcognates <- list()
}




#removeAnticodon(input,"AAG")

#print(input)
#ggplot(data = key, aes(x = as.integer(key[,6]), y = as.integer(key[,8])) +
#geom_point()
#print(input[x,])

input <- read_tsv(file = "ecoli_tRNA_2010.tsv")
input <- input %>% rowwise() %>% 
  mutate(AntiCodon = reverseComplement(Codon)) %>% 
  mutate(AntiCodon=ifelse(Codon == "ATA","CAT",AntiCodon)) ## This is because of some weirdness with bacteria. 
key <- read_tsv(file = 'shah_gilchrist_2010_ecoli.tsv')

#### Start of Alex's additions ####

df <- data.frame(Codon = input$Codon,AA=input$AA,Cognate=NA,Pseudo.cognate=NA,Near.cognate=NA)
rownames(df) <- df$Codon


for(val in 1:nrow(input))
{
  codon <- input[val,"Codon"] %>% deframe()
  codon.aa <- input[val,"AA"] %>% deframe()
  anticodon <- input[val,"AntiCodon"] %>% deframe()
  reverse.Comp <- reverseComplement(codon)
    
  ## Get anticodon neighbors
  one.step.neighbors <- getNeighbors(reverse.Comp)
  neighbors.df <- input %>% 
    filter((AntiCodon %in% one.step.neighbors))
  anticodon.neighbors <- neighbors.df$AntiCodon
  names(anticodon.neighbors) <- anticodon.neighbors

  
  ## Find the anticodons that could match codon via wobble
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
        #return(data.frame(Wobble=TRUE))
      } 
    }
    tmp
  }) %>% bind_rows() %>% distinct()
  neighbors.df <- neighbors.df %>% left_join(wobbles,by=c("AntiCodon","AA"))
  
  ## Label codons as cognate, pseudo-cognate, or near-cognate
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

print(df)
#### End of Alex's additions ####
library (tidyverse)
merge.df <- df %>% left_join(key, by= c("AA", "Codon"))
merge.df <- merge.df %>% mutate(Cognates = ifelse(is.na(Cognates), "", Cognates))
merge.df <- merge.df %>% mutate(`Pseudo-cognates` = ifelse(is.na(`Pseudo-cognates`), "", `Pseudo-cognates`))
merge.df <- merge.df %>% mutate(`Near-cognates` = ifelse(is.na(`Near-cognates`), "", `Near-cognates`))
merge.df <- merge.df %>% mutate(`Near-cognates` = ifelse(is.na(`Near-cognates`), "", `Near-cognates`))
merge.df <- merge.df %>% mutate(`Near-cognates` = ifelse(is.na(`Near-cognates`), "", `Near-cognates`))

for(i in 1:nrow(merge.df)){
  
  row <- merge.df[i,]
  check <- check.same.codons(row$Cognate, row$Cognates)
  if(check == FALSE){
    print(row$Codon)
    
  }
}

for(i in 1:nrow(merge.df)){
  
  row <- merge.df[i,]
  check <- check.same.codons(row$Near.cognate, row$`Near-cognates`)
  if(check == FALSE){
    print(row$Codon)
    
  }
}

for(i in 1:nrow(merge.df)){
  
  row <- merge.df[i,]
  check <- check.same.codons(row$Pseudo.cognate, row$`Pseudo-cognates`)
  if(check == FALSE){
    print(row$Codon)
    
  }
}

# check to see if it is longer than empty string. see if anything is in psudeo  cognate(otherwise skip)
#merge.df <- left_join(df, input)
#merge.df$Rc <- rep(NA,nrow(merge.df))
for(i in 1:nrow(merge.df)){
  row = merge.df[i,]
  print(row$Codon)
  merge.df[i, "Rc_my"] <- EquationThree(row$Codon, row$Cognate, row$Pseudo.cognate, input, row$AA)
  merge.df[i, "Rn_my"] <- EquationFour(row$Codon, row$Near.cognate, input, row$AA)
  merge.df[i, "eM_my"] <- EquationTwo(row$Codon)
  merge.df[i, "eN_my"] <- EquationOne(row$Codon)
  
  #^used to add Rc values^
}

Scatter.Plot <- function(datas){
  #need to create column called Rc_key

  
}

p1 <- ggplot(merge.df, aes(x = Rc, y = Rc_my, label = Codon)) +
  geom_text()+
  geom_point()+
  labs(
    x = "Key",
    y = "My Estimate",
    #I need to actually create a column and add the values from equations 1-4 into the df for top part to work. Im not sure how to.
    title = "Comparing Elongation Rates to Shah And Gilchrist 2010"
  )+ 
  geom_abline(intercept = 0, slope=1, linetype = "dashed")+
  theme_bw()


p2 <- ggplot(merge.df, aes(x = Rn, y = Rn_my, label = Codon)) +
  geom_text()+
  geom_point()+
  labs(
    x = "Key",
    y = "My Estimate",
    #I need to actually create a column and add the values from equations 1-4 into the df for top part to work. Im not sure how to.
    title = "Comparing Near Cognate Elongation Rates to Shah And Gilchrist 2010"
  )+ 
  geom_abline(intercept = 0, slope=1, linetype = "dashed")+
  theme_bw()
