#! /usr/bin/perl -w

#This will be used for $i for each individual sample
#open INFILE, "prefix_odd.txt";
#open INFILE, "prefix_file.txt";
open INFILE, "prefix_put1.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;

#Perl script to automate V4 tag sequencing QC process
#last updated 09-29-2017,  Sarah K. Hu
# First step is to perform an initial QC. This step will take as input FWD and REV reads, and will output Quality filtered reads with and without mates, PE vs SE respectively.
#print "source activate qiime1\n\n";
#print "perl create_map.pl\n\n";
#print "mv *_map.txt* mapfiles/\n\n";

for $i(@prefix) {
##Place print statements of tasks to execute within this loop.
# Unzip all fastq files.(done)
    print "/Users/geridollison/Documents/Bioinformatics/Useful-Scripts/./mga_osx /Users/geridollison/Dropbox/PhageSeqs/",$i,"_final.fasta > /Users/geridollison/Documents/Postdoc/PhageAnnotation/GeneCallerTables/MGAOut_",$i,".txt\n";
# Produce the prodigal table
    print "prodigal -a /Users/geridollison/Dropbox/PhageSeqs/outputs/ProdigalAA_",$i,".faa -i /Users/geridollison/Dropbox/PhageSeqs/",$i,"_final.fasta\n";
    print "grep \">\" /Users/geridollison/Dropbox/PhageSeqs/outputs/ProdigalAA_",$i,".faa | sed 's/ \# /,/g' > /Users/geridollison/Documents/Postdoc/PhageAnnotation/GeneCallerTables/ProdigalOut_",$i,".csv\n\n";

# Produce genemarks table
#placeholder

# Produce glimmer table
#placeholder
}

## Eventually connect to R automation scripts.
