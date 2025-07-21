#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  if (!require("tidyverse")){
    install.packages("tidyverse", dependencies = TRUE)
    require("tidyverse")} 
  # Load required packages
  library(tidyverse)
})

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage (in order!): Rscript run_compile.R GeneMarks.txt Prodigal.txt MetaGene.txt Glimmer.txt")
}

# Assign to variables
genemarks_file <- args[1]
prodigal_file  <- args[2]
metagene_file  <- args[3]
glimmer_file   <- args[4]


## Four scripts that are used to make the ID information for Column 1 of the annotation spreadsheet.

# 1. GeneMarksS as input --------------------------------------------------------
RowIDproducer.GMS = function(filepath) {
  # This function produces the column 1 information using (the filepath to) an output .txt file from GeneMarks as input.
  read.table(filepath, header = TRUE, skip = 9) %>% 
    set_names("Gene","Strand","LeftEnd","RightEnd","GeneLength","Class") %>% 
    mutate(LeftEnd = str_replace_all(LeftEnd, "[><]",""),
           RightEnd = str_replace_all(RightEnd, "[><]",""),
           Column1 = str_c(">gene_",Gene,"|GeneMark.hmm|",GeneLength,"_aa|",Strand,"|",LeftEnd,"|",RightEnd, sep = "")) %>% 
    select(Column1) %>% 
    return(.)
}
#Sample usage.
#filepath1 = "~/Documents/Postdoc/PhageAnnotation/GeneCallerTables/GMSOut_put1.txt"
#RowIDproducer.GMS(filepath1)


# 2. Prodigal as input --------------------------------------------------
# Note this script should be used in combination with bash commands used to snatch the headers off of the prodigal fasta file.
# input for this is the result of the names from a fasta file. for example: grep ">" sample_prodigal.faa | sed 's/ # /,/g' > Put1Prodigalheaders.csv
# script to produce column1 id for prodigal headers file
RowIDproducer.Prod = function(filepath){
  read.csv(filepath, header = F) %>% 
    mutate(V3 = str_replace_all(V3, "<",""),
           V2 = str_replace_all(V2, "<",""),
           Seqlength = as.numeric(V3) - as.numeric(V2) + 1,
           Strand = ifelse(V4 > 0, "+", "-"),
           Column1 = paste(V1,"|Prodigal|",Seqlength,"_aa|",Strand,"|",V2,"|",V3, sep = "")) %>% 
    select(Column1) %>% 
    return(.)
}
#filepath2 = "~/Documents/Postdoc/PhageAnnotation/GeneCallerTables/ProdigalOut_Put1.csv"
#RowIDproducer.Prod(filepath2)

# 3. MetaGene as input -------------------------------------------------------
## Function to produce column1 id for MetaGene output

RowIDproducer.MetaG = function(filepath) {
  ## Function to produce column1 id for MetaGene output
  read.table(filepath, header = F, skip = 3) %>% 
    mutate(V3 = str_replace_all(V3, "<","") %>% as.numeric(),
           V2 = str_replace_all(V2, "<","") %>% as.numeric(),
           len = V3 - V2 +1,
           Column1 = paste(">",V1,"|MetaGene|",len,"_aa|",V4,"|",V2,"|",V3, sep = "")) %>% 
    select(Column1) %>% 
    return(.)
}
#filepath3 = "~/Dropbox/Annotation_R_working_directory/GeneCallerTables/MetaGeneOut_Put1.txt"
#RowIDproducer.MetaG(filepath3)

# 4. Glimmer output ----------------------------------------------------------
RowIDproducer.Glimm = function(filepath) {
  # Produce the ID column from a glimmer table output
  read.delim(filepath, skip = 1, header = F, sep = "") %>% 
    mutate(V3 = str_replace_all(V3, "<","") %>% as.numeric(),
           V2 = str_replace_all(V2, "<","") %>% as.numeric(),
           Begin = ifelse(V2 > V3, V3, V2),
           End = ifelse(V2 > V3, V2, V3),
           Strand = ifelse(V4 <0, "-", "+"),
           Len = End - Begin +1,
           Column1 = paste(">",V1,"|Glimmer|",Len,"_aa|",Strand,"|",Begin,"|",End, sep = "")) %>% 
    select(Column1) %>% 
    return(.)
}
#filepath4 = "~/Documents/Postdoc/PhageAnnotation/GeneCallerTables/GlimmerOut_Put1.txt"
#RowIDproducer.Glimm(filepath4) #%>% getCoordsfromRowID()


# Get gene coords from RowIDs (product of RowIDproducer functions--------------------------------------------------------
getCoordsfromRowID = function(RowID){
  ## Usage: Function to get the gene coordinates from any Column 1 ID output
  RowID %>% 
    separate(col = Column1, sep = "\\|", into = c("d","f","g","h","j","k")) %>% 
    unite(coords, c("j","k"), sep = "-") %>% 
    pull(coords)
}

# Combine tables; keep RowIDs with unique coordinates---------------------------------------------


Compile_uniqe_genes = function(GeneMarks, Prodigal, MetaGene, Glimmer) {
  # function takes as input one putida phage name, and a parent directory for the For example put1, put13, putCompost, etc. Returns genes with unique coordinates and scores.
  
  #Process gene caller tables  
  rbind(RowIDproducer.GMS(GeneMarks) %>% 
          mutate(Coords = getCoordsfromRowID(.), #Place the coordinates next to the gene for filtering purposes
                 Genemarks = Coords %in% (RowIDproducer.GMS(GeneMarks) %>% getCoordsfromRowID()), # Gene in GMS? All should be true.
                 Prodigal = Coords %in% (RowIDproducer.Prod(Prodigal) %>% getCoordsfromRowID()), # Gene in prodigal?
                 Glimmer = Coords %in% (RowIDproducer.Glimm(Glimmer) %>% getCoordsfromRowID()), #placeholder until I can get Glimmer up and running
                 Metagene = Coords %in% (RowIDproducer.MetaG(MetaGene) %>% getCoordsfromRowID())),
        RowIDproducer.Prod(Prodigal) %>% 
          mutate(Coords = getCoordsfromRowID(.),
                 Genemarks = Coords %in% (RowIDproducer.GMS(GeneMarks) %>% getCoordsfromRowID()),
                 Prodigal = Coords %in% (RowIDproducer.Prod(Prodigal) %>% getCoordsfromRowID()),
                 Glimmer = Coords %in% (RowIDproducer.Glimm(Glimmer) %>% getCoordsfromRowID()),
                 Metagene = Coords %in% (RowIDproducer.MetaG(MetaGene) %>% getCoordsfromRowID())) %>% 
          filter(Genemarks == FALSE),
        RowIDproducer.MetaG(MetaGene) %>% 
          mutate(Coords = getCoordsfromRowID(.),
                 Genemarks = ifelse(Coords %in% (RowIDproducer.GMS(GeneMarks) %>% getCoordsfromRowID()),1,0),
                 Prodigal = ifelse(Coords %in% (RowIDproducer.Prod(Prodigal) %>% getCoordsfromRowID()),1,0),
                 Glimmer = ifelse(Coords %in% (RowIDproducer.Glimm(Glimmer) %>% getCoordsfromRowID()),1,0),
                 Metagene = ifelse(Coords %in% (RowIDproducer.MetaG(MetaGene) %>% getCoordsfromRowID()),1,0)) %>% 
          filter(Genemarks == FALSE, Prodigal == FALSE),
        RowIDproducer.Glimm(Glimmer) %>% 
          mutate(Coords = getCoordsfromRowID(.),
                 Genemarks = ifelse(Coords %in% (RowIDproducer.GMS(GeneMarks) %>% getCoordsfromRowID()),1,0),
                 Prodigal = ifelse(Coords %in% (RowIDproducer.Prod(Prodigal) %>% getCoordsfromRowID()),1,0),
                 Glimmer = ifelse(Coords %in% (RowIDproducer.Glimm(Glimmer) %>% getCoordsfromRowID()),1,0),
                 Metagene = ifelse(Coords %in% (RowIDproducer.MetaG(MetaGene) %>% getCoordsfromRowID()),1,0)) %>% 
          filter(Genemarks == FALSE, Prodigal == FALSE, Metagene == FALSE)) %>%
    separate(Coords, into = c("Gene_range1", "Gene_range2"), sep = "-", remove = TRUE) %>% 
    mutate(Strand = ifelse(str_detect(Column1, "\\|\\+\\|"), "pos", "neg")) %>% 
    arrange(Strand, as.numeric(Gene_range1)) %>% 
    #arrange(as.numeric(Gene_range1)) %>% 
    select(Column1, Strand, Gene_range1, Gene_range2, everything()) %>% 
    group_by(Column1) %>% 
    mutate(Score = sum(Genemarks, Prodigal, Glimmer, Metagene),
           Gene_range1 = as.numeric(Gene_range1),
           Gene_range2 = as.numeric(Gene_range2)) %>%
    group_by(Strand) %>%
    arrange(Strand, ifelse(Strand == "pos", Gene_range1, -Gene_range2)) %>%
    #Calculate gaps and overlaps, (strand conscious)
    mutate(prev_end   = lag(Gene_range2),
           prev_start = lag(Gene_range1),
           Gap = if_else(!is.na(prev_end) & Gene_range2 < prev_start & Strand == "neg", prev_start - Gene_range2,
                         if_else(!is.na(prev_end) & Gene_range1 > prev_end & Strand == "pos",Gene_range1 - prev_end ,0)),
           Overlap = if_else(!is.na(prev_end) & Gene_range2 > prev_start & Strand == "neg", Gene_range2 - prev_start,
                             if_else(!is.na(prev_end) & Gene_range1 < prev_end & Strand == "pos", Gene_range1 - prev_end,0))) %>%
    ungroup() %>% 
    #Handle duplicate ORFs: prioritize orfs with highest score, secondarily, if a tie, the longest orf
    mutate(orf_length = abs(Gene_range2 - Gene_range1) + 1) %>% 
    group_by(Gene_range1) %>%
    # Step 2: For each group (i.e., same coordinates), pick the one with the highest Score
    slice_max(order_by = Score, with_ties = TRUE) %>%
    slice_max(order_by = orf_length, with_ties = FALSE) %>% 
    ungroup() %>% 
    arrange(Strand, ifelse(Strand == "pos", Gene_range1, -Gene_range2)) %>% 
    #select(-starts_with("prev")) %>%
    # Calculate Scores: operon, overlap, length
    group_by(Column1) %>% 
    #Calculate Operon Score
    mutate(Operon_score0 = ifelse(Overlap == 0, 1, NA_integer_),
           Operon_score3 = ifelse(Overlap == -3, 1, NA_integer_),
           Operon_score6 = ifelse(Overlap == -6, 1, NA_integer_),
           OPERON_SCORE = sum(Operon_score0,Operon_score3,Operon_score6, na.rm = TRUE)) %>% 
    #Calculate Length Score
    mutate(Len_200plus = ifelse(orf_length >= 200, 0, NA_integer_),
           Len_150.199 = ifelse(orf_length >= 150 & orf_length < 200, 1, NA_integer_),
           Len_120.149 = ifelse(orf_length >= 120 & orf_length < 150, 2, NA_integer_),
           Len_90.119 = ifelse(orf_length >= 90 & orf_length < 120, 3, NA_integer_),
           Len_sub90 = ifelse(orf_length < 90, 4, NA_integer_),
           Len_Score = sum(Len_200plus, Len_150.199, Len_120.149, Len_90.119, Len_sub90,na.rm = TRUE)) %>% 
    #Calculate Overlap Score
    mutate(Overlap_10.0 = ifelse(Overlap <= 10, 0, NA_integer_),
           Overlap_10.39 = ifelse(Overlap > 10 & Overlap <= 39, 0, NA_integer_),
           Overlap_40.69 = ifelse(Overlap > 40 & Overlap <= 69, 0, NA_integer_),
           Overlap_70.99 = ifelse(Overlap > 70 & Overlap <= 99, 0, NA_integer_),
           Overlap_100plus = ifelse(Overlap > 100, 4, NA_integer_),
           Overlap_Score = sum(Overlap_10.0, Overlap_10.39, Overlap_40.69, Overlap_70.99, Overlap_100plus, na.rm = TRUE)) %>% 
    ungroup() %>% 
    return(.)
}
# source("functions.R")  # <- If you'd rather keep them in a separate file

# Run the function
final_df <- Compile_uniqe_genes(
  GeneMarks = genemarks_file,
  Prodigal = prodigal_file,
  MetaGene = metagene_file,
  Glimmer = glimmer_file
)

# Save the result
write_csv(final_df, "StructureCompiled.csv")