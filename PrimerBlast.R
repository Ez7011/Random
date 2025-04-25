
# this function takes two input, typein the DNA sequence, and the primer list file including its path
PrimerBlast <- function (DNAsequence,primerfile) {
library(Biostrings)
library(readxl)

# standarize all the inputs, get the reverse complement strand 
  DNA <- toupper(DNAsequence)
  AND <- chartr("ATGC","TACG",DNA)
  revDNA <- paste0(rev(strsplit(AND, NULL)[[1]]), collapse = "")
# clean up and standardize the primer table
  primerlist <- suppressMessages(read_excel(primerfile, col_types = "text"))
  primerlist_std <- primerlist[,c("OLIGO NAME","SEQUENCE", "NUMBER")]
  primerlist_std$SEQUENCE <- gsub("[^ACTG]", "", toupper(trimws(primerlist_std$SEQUENCE )))
  primerlist_std <- primerlist_std[!is.na(primerlist_std$SEQUENCE) & primerlist_std$SEQUENCE != "", ]
  
# blast for fwd primers
  results_fwd <- lapply(primerlist_std$SEQUENCE, function(query) {
    
    fwd <- matchPattern(DNAString(query), DNA)
      if (length(fwd) > 0) {
        data.frame(
          ID = primerlist_std$NUMBER[primerlist_std$SEQUENCE == query],
          query_sequence =  as.character(query),
          start_position =  start(fwd),
          end_position =  end(fwd),
          width =  width(fwd),
          direction = "forward"
          )
      }
    
  }
  )
# blast for reverse primers
  results_rvs <- lapply(primerlist_std$SEQUENCE, function(query) {
    
    rvs <- matchPattern(DNAString(query), revDNA)
    if (length(rvs) > 0) {
      data.frame(
        ID = primerlist_std$NUMBER[primerlist_std$SEQUENCE == query],
        query_sequence =  as.character(query),
        start_position =  start(rvs),
        end_position =  end(rvs),
        width =  width(rvs),
        direction = "reverse"
      )
    }
  }
  )
    
    
# Combine all into one dataframe
  match_df <- do.call(rbind, c(results_fwd, results_rvs))

  match_df
}
