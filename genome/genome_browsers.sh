#!/bin/bash
input=/home/sobczm/popgen/genome
scripts=/home/sobczm/bin/popgen/genome
browsers=/home/sobczm/bin/browsers

#Install UCSC, GBrowse and Ensembl genome browsers for ananassa latest genome.

#Copy input genome. (Latest D. Swarbreck annotation)
cp -r /home/sobczm/popgen/rnaseq/ananassa_annotation $input

#USCS
cd $browsers
git clone http://genome-source.cse.ucsc.edu/kent.git
#Most users want to use the beta branch, which has tested, released versions of the browser. To create a beta tracking branch:
cd kent
git checkout -t -b beta origin/beta
#The -b creates a local branch with the name "beta", and checks it out.
#The -t makes it a tracking branch, so that "git pull" pulls in updates from origin/beta, the "real" beta branch in our public central read-only repository.

#To get the latest UCSC release:
git pull