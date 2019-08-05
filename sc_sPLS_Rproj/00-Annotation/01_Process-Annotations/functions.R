# Function to process de BioMart results
# - change strand from -1/+1 to -/+
# - sort using genomic coordinates
processBioMart <- function(df, chrnames) {
  # chrnames: name of the chromosomes
  # colnames(df) <- coln
  df$strand <- c("-","+")[factor(df$strand)]
  # df$chr <- factor(paste("chr",df$chr,sep=""), levels=chrnames)
  df <- df[with(df, order(chr,start,end)), ]
  return(df)
}

convert_chr_format <- function(chr, to) {
  # Function to convert the chr from short to long format and viceversa
  # to: "short" or "long"
  chr <- as.character(chr)
  stopifnot(to %in% c("short","long"))
  # short_alphabet <- c(1:19,"X","Y","MT")
  short_alphabet <- c(1:22,"X","Y","MT")
  long_alphabet <- paste("chr",short_alphabet,sep="")
  if (to == "short") {
    if (all(chr %in% short_alphabet)) { 
      return(factor(chr)) 
    } else {
      stopifnot(all(chr %in% long_alphabet))
      names(short_alphabet) <- long_alphabet
      return(factor(unname(short_alphabet[chr])))
    }
  }
  if (to == "long") {
    if (all(chr %in% long_alphabet)) { 
      return(factor(chr)) 
    } else {
      stopifnot(all(chr %in% short_alphabet))
      names(long_alphabet) <- short_alphabet
      return(factor(unname(long_alphabet[chr])))
    }
  }
}

# Function to filter non-interesting genes by description
remove <- function(df, riken=T, no_symbol=F, dna_segments=T, empty_description=T, 
                   duplicated=F, pseudogenes=T, predicted=T, expressed_sequence=T,
                   olfactory=T) {
  
  # Remove the set of useless predicted genes
  if (predicted == T)
    df <- df %>% filter(!str_detect(description,"predicted "))

  # Remove genes with no symbol
  if (no_symbol == T) {
    df <- df %>% filter(!symbol=="")
  }
  
  # Remove genes with no description
  if (empty_description == T) {
    df <- df %>% filter(!description=="")
  }
  
  # Remove RIKEN cDNA
  if (riken == T) {
    df <- df %>% filter(!str_detect(description,"RIKEN cDNA"))
    df <- df %>% filter(!str_detect(description,"Riken cDNA"))
  }
  
  # Remove useless DNA segments
  if (riken == T) {
    df <- df %>% filter(!str_detect(description,"DNA segment"))
    df <- df %>% filter(!str_detect(description,"cDNA sequence"))
  }

  # Remove pseudogenes
  if (pseudogenes == T) {
    df <- df %>% filter(!str_detect(description,"pseudogene"))
  }
  
  # Remove 'expressed sequence XXX'
  if (pseudogenes == T) {
    df <- df %>% filter(!str_detect(description,"expressed sequence"))
  }
  
  # Remove duplicated elements
  if (duplicated == T) {
    dup <- df$ens_id[duplicated(df$ens_id)]
    df <- df %>% filter(!ens_id %in% dup)
  }
  
  # Remove olfactory receptors (why are there so many?????)
  if (olfactory == T) {
    df <- df %>% filter(!str_detect(description,"olfactory receptor"))
  }
  
  return(df)
}
