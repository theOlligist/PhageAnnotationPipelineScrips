#!/usr/bin/env Rscript

#Load in packages
suppressMessages({
  library(tidyverse)
  library(Biostrings)})

#Enable command line argument handling
args = commandArgs(trailingOnly = TRUE) #trailingonly makes only the defined arbuments?

#error handling to ensure two arguments are given
if(length(args) != 2) {
  stop("Rscript make_bed.R genome.fna SequenceCompiled.csv")
}

#Define the input arguments
fasta_file <- args[1]
orf_table  <- args[2]
output_bed <- "ORFs.bed"

#Load genome fasta and get the name of the genome for the contig information

genome=readDNAStringSet(fasta_file)
contig_name = names(genome)[1]

#Load in StructureCompiled for the ORFs

StructureCompiled = read.csv(orf_table)

# Build Bed format file
bed_df <- StructureCompiled %>%
  mutate(
    gene_id = str_remove(Column1, "^>"),  # use Column1 or extract short ID
    contig  = contig_name,
    start   = pmin(Gene_range1, Gene_range2) - 1,  # 0-based for BED
    end     = pmax(Gene_range1, Gene_range2),
    strand  = case_when(
      Strand %in% c("+", "pos", "plus") ~ "+",
      Strand %in% c("-", "neg", "minus") ~ "-",
      TRUE ~ "."
    ),
    score = 0
  ) %>%
  select(contig, start, end, gene_id, score, strand)

# ---- Write to BED file ----
write_tsv(bed_df, output_bed, col_names = FALSE)
message("✅ BED file written to: ", output_bed)