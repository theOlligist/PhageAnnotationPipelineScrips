# Created by theOlligist: 7Mar2024
# Function to read a FASTA file line by line and produce a list containing the HEADERS and SEQUENCES in separate vectors.

MakeSeqDictionary <- function(file_path) {
  #load tidyverse for replace and concatinate functions
  library(tidyverse)
  #Create Headers and sequence as empty containers to populate while reading the input
  HEADERS <- c()
  SEQUENCES <- c()
  
  #Establish a connection to the input file (should be a fasta)
  file_connection <- file(file_path, "r")
  #While loop to go through the entire file, line by line. While there are lines to read, do the following to each line.
  while (TRUE) {
    #read the first line (n=1) do some tidying up and save the text string as an object called currline
    currline = readLines(file_connection, n = 1, warn = F) %>% 
      str_replace_all(., "(\\s)|(\\#)|(\\\t>)","|") %>% 
      str_replace_all(., "\\|\\|\\|", "|") %>% 
      str_replace_all(., "\\*", "") 
    #if the length of the currline is 0 kill the while loop. We have arrived at the end
    if(length(currline) == 0) {
      print("done")
      break()}
    #if the first character of the line is a >, do these two things, because it's a header
    if(substr(currline, 1, 1) == ">") {
      HEADERS <- c(HEADERS, currline) #add it to the headers list.
      SEQUENCES <- c(SEQUENCES, "") #create an empty entry in the SEQUENCES vector
    } else if (nchar(currline) > 1) { #Otherwise do the following, but only if there is more than one character in currline
      SEQUENCES[length(SEQUENCES)] <- str_c(SEQUENCES[length(SEQUENCES)], currline, sep = "")  #here I am opting to concatinate strings of sequences together
    }
  }
  #once out of the while loop, close the file connection and stuff HEADERS and SEQUENCES into one list where they will parallel each other.
  close(file_connection)
  return(list(HEADERS = HEADERS, SEQUENCES = SEQUENCES))
}

#usage
#provide the path to a fasta formatted file.
#seqdictionary = MakeSeqDictionary(file_path)

#Function to filter a seqdictionary for a vector of headernames.
GetSeqs = function(seqdictionary, search_vector){
  seqdictionary_filt = list(HEADERS = seqdictionary[[1]][seqdictionary[[1]] %in% search_vector], SEQUENCES = seqdictionary[[2]][seqdictionary[[1]] %in% search_vector])
  return(seqdictionary_filt)
}
#usage
#GetSeqs(seqdictionary, c(">gene_1|GeneMark.hmm|90_aa|-|2|271|>put5_pilon", ">gene_11|GeneMark.hmm|66_aa|-|9151|9351|>put5_pilon"))
