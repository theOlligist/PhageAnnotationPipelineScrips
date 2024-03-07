# PhageAnnotationPipelineScripts
Scripts that, once finished, will combine to make a pipeline for automating single phage genome annoataion.

Annotation_scripts_v2.R:
Used to pull the unique open reading frames from the four gene calling software outputs. (needs to be redone. I want them to use the protein fasta headers to do this instead of the tables. The way this script handles prodigal inputs is an an example. The other three use the output tables. I want to make this change bcause MakeSeqDictionary will use the fasta files as well.)

MakeSeqDictionary.R:
These two functions take a fasta as input and stores the Header and Sequence in a list (aka dictionary) parallel to one another. This enables me to get the sequences from the name, sorta like a key-value pair.

************************* Unfinished Business
To do upstream:
-Automate the production of the four tables for use in the annotation_scripts_v2.R and their concatenation via bash and perl scripts. The concatenated version will be passed to the MakeSeqDirectionary.R which will be filtered.
Test to perfom: 
$cat GMS_output.faa Prodigal_output.faa > combined.faa
$MakeSeqDictionary(combined.faa) does it work?

To do downstream:
-Do something with the sequences from the unique ORFs. I can either place them in a column of the dataframe with the other information OR write a separate filtered and ordered fasta file. These will enable the enduser to blastp these sequences for functional annotation.
