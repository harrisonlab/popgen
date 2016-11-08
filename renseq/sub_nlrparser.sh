#!/bin/bash

fasta=$1
fasta_h=$(basename "$fasta")

### Setting variables
nlr=/home/sobczm/bin/NLR-Parser
meme=/home/sobczm/bin/meme_4.9.1/bin

chop="${fasta_h%.*}_chopped"
java -jar $nlr/ChopSequence.jar -i $fasta_h -o $chop -l 20000 -p 5000

#Then you run NLR-Parser on the choppedFasta.
#You will need meme-4.9.1 (not a higher version).
#You can download the meme.xml from github NLR-Parser
java -jar $nlr/NLR-Parser.jar -i $chop -x $nlr/meme.xml -y $meme/mast \
-c $chop.nlr.xml

#Then you run the NLR-Annotator.
#Currently I have TSV, BED, GFF output and a few undocumented experimental things.
#Let me know if you want anything else. The BED outputs have the most information.
java -jar $nlr/NLR-Annotator.jar -i $chop.nlr.xml -o "${fasta_h%.*}"_nlr.tsv \
-g "${fasta_h%.*}"_nlr.gff -b "${fasta_h%.*}"_nlr.bed -m "${fasta_h%.*}"_nlr.motifs.bed
