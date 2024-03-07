
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
#filepath3 = "~/Documents/Postdoc/PhageAnnotation/GeneCallerTables/MetaGeneOut_Put1.txt"
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

Compile_uniqe_genes = function(PUT) {
  # function takes as input one putida phage name. For example put1, put13, putCompost, etc. Returns genes with unique coordinates and scores.
  filepath.a = paste("~/Dropbox/Annotation_R_working_directory/GeneCallerTables/GMSOut_",PUT,".txt", sep = "")
  filepath.b = paste("~/Dropbox/Annotation_R_working_directory/GeneCallerTables/ProdigalOut_",PUT,".csv", sep = "")
  filepath.c = paste("~/Dropbox/Annotation_R_working_directory/GeneCallerTables/MGAOut_",PUT,".txt", sep = "")
  filepath.d = paste("~/Dropbox/Annotation_R_working_directory/GeneCallerTables/GlimmerOut_",PUT,".txt", sep = "") 

  rbind(RowIDproducer.GMS(filepath.a) %>% 
          mutate(Coords = getCoordsfromRowID(.), #Place the coordinates next to the gene for filtering purposes
                 Genemarks = Coords %in% (RowIDproducer.GMS(filepath.a) %>% getCoordsfromRowID()), # Gene in GMS? All should be true.
                 Prodigal = Coords %in% (RowIDproducer.Prod(filepath.b) %>% getCoordsfromRowID()), # Gene in prodigal?
                 Glimmer = Coords %in% (RowIDproducer.Glimm(filepath.d) %>% getCoordsfromRowID()), #placeholder until I can get Glimmer up and running
                 Metagene = Coords %in% (RowIDproducer.MetaG(filepath.c) %>% getCoordsfromRowID())),
        RowIDproducer.Prod(filepath.b) %>% 
          mutate(Coords = getCoordsfromRowID(.),
                 Genemarks = Coords %in% (RowIDproducer.GMS(filepath.a) %>% getCoordsfromRowID()),
                 Prodigal = Coords %in% (RowIDproducer.Prod(filepath.b) %>% getCoordsfromRowID()),
                 Glimmer = Coords %in% (RowIDproducer.Glimm(filepath.d) %>% getCoordsfromRowID()),
                 Metagene = Coords %in% (RowIDproducer.MetaG(filepath.c) %>% getCoordsfromRowID())) %>% 
          filter(Genemarks == FALSE),
        RowIDproducer.MetaG(filepath.c) %>% 
          mutate(Coords = getCoordsfromRowID(.),
                 Genemarks = ifelse(Coords %in% (RowIDproducer.GMS(filepath.a) %>% getCoordsfromRowID()),1,0),
                 Prodigal = ifelse(Coords %in% (RowIDproducer.Prod(filepath.b) %>% getCoordsfromRowID()),1,0),
                 Glimmer = ifelse(Coords %in% (RowIDproducer.Glimm(filepath.d) %>% getCoordsfromRowID()),1,0),
                 Metagene = ifelse(Coords %in% (RowIDproducer.MetaG(filepath.c) %>% getCoordsfromRowID()),1,0)) %>% 
          filter(Genemarks == FALSE, Prodigal == FALSE),
        RowIDproducer.Glimm(filepath.d) %>% 
          mutate(Coords = getCoordsfromRowID(.),
                 Genemarks = ifelse(Coords %in% (RowIDproducer.GMS(filepath.a) %>% getCoordsfromRowID()),1,0),
                 Prodigal = ifelse(Coords %in% (RowIDproducer.Prod(filepath.b) %>% getCoordsfromRowID()),1,0),
                 Glimmer = ifelse(Coords %in% (RowIDproducer.Glimm(filepath.d) %>% getCoordsfromRowID()),1,0),
                 Metagene = ifelse(Coords %in% (RowIDproducer.MetaG(filepath.c) %>% getCoordsfromRowID()),1,0)) %>% 
          filter(Genemarks == FALSE, Prodigal == FALSE, Metagene == FALSE)) %>%
    separate(Coords, into = c("Gene_range1", "Gene_range2"), sep = "-", remove = TRUE) %>% 
    arrange(as.numeric(Gene_range1)) %>% 
    select(Column1, Gene_range1, Gene_range2, everything()) %>% 
    return(.)
}
#Compile_uniqe_genes("put1")

# Write a table to copy and paste to the lab excel spreadsheet ------------

write_compiled_table = function(PUT,PATH = "~/Desktop/"){
  Compile_uniqe_genes(PUT) %>% 
    write_delim(file = paste(PATH,"Compiled",PUT,".txt", sep = ""), delim = "\t")
}
#write_compiled_table("put1","~/Desktop/")

#End.