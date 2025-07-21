# PhageAnnotationPipelineScripts
Automate_GeneCalling.pl: Optional script that is useful for looping over multiple files and jobs.

CompileStructure.R (Structural Annotation; Gets comprehensive list of Open Reading frames)
Usage: Open reading frames are predicted from a genome using four ORF calling tools: Glimmer, GeneMarkS, Prodigal, and Metagene. From these tools a table is produced where each row represents an ORF and several columns of metadata describing important features of the ORF, such as its coordinates and lengtn. Our goal is to produce a comprehensive list of ORFs utilizing all four program's ORF calls. This is where CompileStructure.R comes in. 
CompileStructure.R is executed from the command line with four inputs as follows:
Rscript CompileStructure.R PATH/TO/GeneMarkS_table.txt PATH/TO/Prodigal_table.csv PATH/TO/MetaGene_table.txt PATH/TO/Glimmer_table.txt

The program will produce a csv called StructureCompiled.csv

/ ! \ In this version the order is important: GeneMarkS, Prodigal, Metagene, Glimmer.
/ ! \ Note that prodigal produces a csv, not a txt

### functional annotation
We need ORF sequences now.
StructureCompiled.csv serves as raw material for the sequences at the positions noted by Gene_range1 and Gene_range2, in a strand conscious way.

bedtools will be employed to produce a fasta of the ORFs from the genome, ORF coordinates and strand information. Unfortunately it requires the OFR coord and strand to be contained in a special file; a .bed file which we will create with Make_Bed.R

Usage:
Rscript Make_Bed.R PATH/TO/genome.fasta PATH/TO/StructureCompiled.csv  -->ORFs.bed

Use the .bed file to get the ORF sequences as a fasta with bedtools:
/ ! \ installing bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html.
For (mac osx): brew install bedtools

#Run Bedtools
bedtools getfasta -fi assembly.fna -bed ORFs.bed -s -name -fo orfs.fa


RunDram-V to annotate the sequences
DRAM-v.py annotate -i orfs.fa -o dramvAnnotation --min_contig_size 100

#Join the annotations to the structure table to produce a complete table with both structural and functional annotations.
left_join(StructureCompiled.csv annotations.tsv)
