# PhageAnnotationPipelineScripts
Scripts that, once finished, will combine to make a pipeline for automating single phage genome annoataion.
Automate_GeneCalling.pl:
Perl script that takes as input a set of genome prefixes and runs each gene calling software on each prefix. It saves these tables and extracts unique ORFs via Annotation_scripts_v2.R

CompileStructure.R
Usage: Open reading frames are predicted from a genome using four ORF calling tools: Glimmer, GeneMarkS, Prodigal, and Metagene. From these tools a table is produced where each row represents an ORF and several columns of metadata describing important features of the ORF, such as its coordinates and lengtn. Our goal is to produce a comprehensive list of ORFs utilizing all four program's ORF calls. This is where CompileStructure.R comes in. 
CompileStructure.R is executed from the command line with four inputs as follows:
Rscript CompileStructure.R PATH/TO/GeneMarkS_table.txt PATH/TO/Prodigal_table.csv PATH/TO/MetaGene_table.txt PATH/TO/Glimmer_table.txt

The program will produce a csv called StructureCompiled.csv

/ ! \ In this version the order is important: GeneMarkS, Prodigal, Metagene, Glimmer
/ ! \ Note that prodigal produces a csv, not a txt


### On the horizon: functional annotation
StructureCompiled.csv serves as input for GetSequence.R

MakeBED.R StructureCompiled.csv ORFs.bed
GetSequence.R(StructureCompiled.csv, genome.fasta) --> StructureFXCompiled.csv


bedtools getfasta -fi assembly.fna -bed ORFs.bed -s -name -fo orfs.fa
write_fasta(Column1, Seq) 

RunDram-V

left_join()
