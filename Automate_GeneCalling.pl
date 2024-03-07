#! /usr/bin/perl -w

#This will be used for $i for each line in the prefix file/each genome prefix
open INFILE, "prefix_odd.txt";

@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;

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
