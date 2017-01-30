## *Programs installed locally on Maria's account:* :smile:
##R packages:
Using R version 3.2.2 installed locally:
```export PATH=/home/armita/prog/R/R-3.2.2/bin:${PATH}```
and libraries stored in ```export R_LIBS=/home/sobczm/R/x86_64-pc-linux-gnu-library/3.2:$R_LIBS```

```diff
- If aiming to use those libraries, append the path to them in the following way:
R
.libPaths( c( .libPaths(), "/home/sobczm/R/x86_64-pc-linux-gnu-library/3.2") )
.libPaths( c( .libPaths(), "/home/armita/prog/R/R-3.2.2/library") )
```
##Programs:
###Source file with all dependencies for the programs below. If in doubt, load all of them into your current shell instance prior to execution of any pipeline by adding the line below to the top of your script:
```source /home/sobczm/bin/marias_profile```

###Alternatively, export each dependency to PATH individually by hand
Type ```nano ~/.profile``` to start editing your BASH profile.
Press ```Alt and /``` to navigate until the end of the file and paste the export command on a new line, for instance:
```export PATH=/home/sobczm/bin/mcl-14-137/bin:${PATH}```. Save changes and exit by pressing ```Ctr and x```followed by Return (AKA Enter). 
For the changes to take place, either type ```source ~/.profile``` or close and re-open the terminal window.

[ade4] (https://cran.r-project.org/web/packages/ade4/index.html): Analysis of Ecological Data : Exploratory and Euclidean Methods in Environmental Sciences

[adegenet] (http://adegenet.r-forge.r-project.org/): a R package for the multivariate analysis of genetic markers

[ggplot2] (http://ggplot2.org/): graphing package implemented on top of the R statistical package

[PCAdapt] (http://membres-timc.imag.fr/Michael.Blum/PCAdapt.html): pcadapt implements a genome scan for detecting genes involved in local adaptation

[pegas] (https://cran.r-project.org/web/packages/pegas/index.html): Population and Evolutionary Genetics Analysis System

[PopGenome] (https://cran.r-project.org/web/packages/PopGenome/index.html): An Efficient Swiss Army Knife for Population Genomic Analyses

[poppr] (http://grunwaldlab.cgrb.oregonstate.edu/poppr-r-package-population-genetics): Genetic Analysis of Populations with Mixed Reproduction

[SNPRelate] (http://corearray.sourceforge.net/tutorials/SNPRelate/): Parallel Computing Toolset for Relatedness and Principal Component Analysis of SNP Data

[vcfR] (https://cran.r-project.org/web/packages/vcfR/index.html): Manipulate and Visualize VCF Data

[WhopGenome] (https://cran.r-project.org/web/packages/WhopGenome/index.html): High-Speed Processing of VCF, FASTA and Alignment Data

##Standalone
[4P] (https://github.com/anbena/4p): 4P (Parallel Processing of Polymorphism Panels) is a software for computing
population genetics statistics from large SNPs dataset. `/home/sobczm/bin/4p/bin`

[ABRA ver. 0.97] (https://github.com/mozack/abra): Improved coding indel detection via assembly-based
realignment `/home/sobczm/bin/abra/bin`

[BayeScan ver. 2.1] (http://cmpg.unibe.ch/software/BayeScan/): detecting natural selection from population-based genetic data
`/home/sobczm/bin/bayescan2.1/binaries/bayescan_2.1`

[Beagle ver. 4.1] (https://faculty.washington.edu/browning/beagle/beagle.html) Beagle is a software package that performs genotype calling, genotype phasing, imputation of ungenotyped markers, and identity-by-descent segment detection. `/home/sobczm/bin/beagle`

**BEAST requires Java 8 - downloaded it to a local directory and changed path for default Java in my profile:**
```
export JAVA_HOME=/home/sobczm/bin/jre1.8.0_101
export PATH="$JAVA_HOME/bin:$PATH"
```
[BEAST ver. 1.8.3] (http://beast.bio.ed.ac.uk/) package - Bayesian analysis of molecular sequences using MCMC. Includes: BEAST, BEAUti, LogCombiner, TreeAnnotator.  `/home/sobczm/bin/beast/BEASTv1.8.3/bin`

[BEAST ver. 2.4.2] (http://beast2.org/) package. Includes: BEAST, BEAUti, LogCombiner, TreeAnnotator, DensiTree. `/home/sobczm/bin/beast/BEASTv2.4.2/bin`

[BEASTGen ver. 1.0.2] (http://beast.bio.ed.ac.uk/beastgen) Creates BEAST XML input files.
`/home/sobczm/bin/beast/BEASTGenv1.0.2/bin`

[bioawk] (https://github.com/lh3/bioawk) BWK awk modified for biological data. `/home/sobczm/bin/bioawk`

[BUSCO ver 2.0] (http://busco.ezlab.org/): Assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs
`/home/sobczm/bin/BUSCO_v2`

**Dependencies**:
```
export PATH=/home/sobczm/bin/hmmer-3.1b2/binaries:${PATH}
export PATH=/home/armita/prog/python3/Python-3.3.5/bin:${PATH}
export PATH=/home/armita/prog/ncbi-rmblastn-2.2.28/bin:${PATH}
export AUGUSTUS_CONFIG_PATH=/home/sobczm/bin/augustus-3.1/config
export PATH=/home/sobczm/bin/augustus-3.1/bin:${PATH}
export PATH=/home/sobczm/bin/augustus-3.1/scripts:${PATH}
export PATH=/home/armita/prog/emboss/EMBOSS-4.0.0/bin:${PATH}
```
[CodonW ver.1.3] (http://codonw.sourceforge.net/) Multivariate analysis (correspondence analysis) of codon and amino acid usage. `/home/sobczm/bin/codonW`

[CViT ver. 1.2.1] (https://sourceforge.net/projects/cvit/) CViT - Chromosome Viewing Tool `/home/sobczm/bin/cvit.1.2.1`

[DAGchainer] (http://dagchainer.sourceforge.net/) DAGchainer: Computing Chains of Syntenic Genes in Complete Genomes `/home/sobczm/bin/DAGCHAINER`

[DendroPy ver. 4.1.0] (https://pythonhosted.org/DendroPy/) Python library for phylogenetic computing. `/home/sobczm/bin/DendroPy` *Need to export it to PYTHONPATH

[DivStat] (http://www.portugene.com/DivStat.html): A User-Friendly Tool for Single Nucleotide Polymorphism Analysis of
Genomic Diversity `/home/sobczm/bin/DivStat`

[EIGENSOFT ver. 6.1.2] (https://genetics.med.harvard.edu/reich/Reich_Lab/Software.html) A set of population structure dectection methods. `/home/sobczm/bin/EIG/bin`

[FastTree ver. 2.1.9] (http://www.microbesonline.org/fasttree/) Approximately Maximum-Likelihood Trees for Large Alignments `/home/sobczm/bin/FastTree2.1.9`  

[FigTree ver. 1.4.2] (http://tree.bio.ed.ac.uk/software/figtree/) Viewing of phylogenetic trees and production of publication-ready figures. `/home/sobczm/bin/FigTree_v1.4.2/bin`

[freebayes ver. v1.0.2] (https://github.com/ekg/freebayes) Bayesian haplotype-based polymorphism discovery and genotyping. `/home/sobczm/bin/freebayes/bin`

[GATK ver. 3.6] (https://software.broadinstitute.org/gatk/) Genome Analysis Toolkit - Variant Discovery in High-Throughput Sequencing Data. `/home/sobczm/bin/GenomeAnalysisTK-3.6`

[GeneProteinViz ver. 1.2.8] (http://www.icbi.at/software/gpviz/gpviz.shtml) Dynamic visualization of genomic regions and variants affecting protein domains. 
`/home/sobczm/bin/GPViz`

[LMAP ver. 1.0] (http://lmapaml.sourceforge.net/) A collection of perl scripts to automate PAML use. Requires a number of dependencies `/home/sobczm/bin/LMAPv1.0.0/LMAP`

[LUMPY] (https://github.com/arq5x/lumpy-sv) A probabilistic framework for structural variant discovery `/home/sobczm/bin/lumpy-sv/bin`

[MAFFT ver. 7.222] (http://mafft.cbrc.jp/alignment/software/) Rapid multiple sequence alignment based on fast Fourier transform
`/home/sobczm/bin/mafft-7.222/bin`

[MEGA ver. 7] (http://www.megasoftware.net/) Sophisticated and user-friendly software suite for analyzing DNA and protein sequence data from species and populations. `/home/sobczm/bin/mega`

[MEME Suite ver. 4.11.2] (http://meme-suite.org/) MEME SUITE: tools for motif discovery and searching `/home/sobczm/bin/meme_4.11.2/bin`

[NLR-Parser] (https://github.com/steuernb/NLR-Parser) A tool to rapidly annotate the NLR complement from sequenced plant genomes. `home/sobczm/bin/NLR-Parser` Requires meme 4.9.1 in `/home/sobczm/bin/meme_4.9.1/bin`

[OrthoFinder v1.0.7] (https://github.com/davidemms/OrthoFinder/) OrthoFinder: Accurate inference of orthogroups, orthologues, gene trees and rooted species tree made easy `/home/sobczm/bin/OrthoFinder-1.0.7/orthofinder`

**Dependencies**:
```
export PATH=/home/sobczm/bin/mcl-14-137/bin:${PATH}
export PATH=/home/sobczm/bin/fastme-2.1.5/bin:${PATH}
export PATH=/home/sobczm/bin/dlcpar-1.0/bin:${PATH}
export PATH=/home/sobczm/bin/mafft-7.222/bin:${PATH}
export PATH=/home/sobczm/bin/FastTree2.1.9:${PATH}

```
[PAL2NAL] (http://www.bork.embl.de/pal2nal/): Robust conversion of protein sequence alignments into the corresponding codon alignments `/home/sobczm/bin/pal2nal.v14`

[PAML ver. 4.8] (http://abacus.gene.ucl.ac.uk/software/paml.html) A package of programs for phylogenetic analyses of DNA or protein sequences using maximum likelihood. `/home/sobczm/bin/paml4.8/bin`

[PartitionFinder ver. 1.1.1] (http://www.robertlanfear.com/partitionfinder/) A program to select best-fit partitioning schemes and models of molecular evolution for phylogenetic analyses `/home/sobczm/bin/PartitionFinder1.1.1` *Needs to be run with Anaconda Python distribution installed in `/home/sobczm/bin/anaconda2/bin`

[PGDSpider ver. 2.1.0.3] (http://www.cmpg.unibe.ch/software/PGDSpider/) An automated data conversion tool for connecting population genetics and genomics programs `/home/sobczm/bin/PGDSpider_2.1.0.3`

[PicardTools ver. 2.5.0] (https://broadinstitute.github.io/picard/) A set of command line tools for manipulating formats such as SAM/BAM/CRAM and VCF. `/home/sobczm/bin/picard-tools-2.5.0`

[PhyloNet ver. 3.5.5] (http://bioinfo.cs.rice.edu/phylonet) & [PhyloNetHMM ver. 0.1] (http://bioinfo.cs.rice.edu/software/phmm) Bayesian inference of reticulate phylogenies under the multispecies network coalescent `/home/sobczm/bin/phmm`

[Phyutility ver. 2.6.6] (https://code.google.com/archive/p/phyutility/) Phyutility provides a set of phyloinformatics tools for summarizing and manipulating phylogenetic trees, manipulating molecular data and retrieving data from NCBI. `/home/sobczm/bin/phyutility`

[popoolation ver. 1.2.2] (https://sourceforge.net/projects/popoolation/) PoPoolation is a pipeline for analysing pooled next generation sequencing data. `home/sobczm/bin/popoolation_1.2.2`

[RAxML ver. 8.2.9] (http://sco.h-its.org/exelixis/web/software/raxml/index.html): a ML a tool for phylogenetic analysis and post-analysis of large phylogenies `/home/sobczm/bin/RAxML8.2.9`

[RGAugury] (https://bitbucket.org/yaanlpc/rgaugury/): A pipeline for genome-wide prediction of resistance gene analogs (RGAs) in plants. **Multiple dependencies, including own copy of Phobius. Phobius output not saved correctly when running via qsub**
`/home/sobczm/bin/rgaugury`

[snpEff & snpSift ver. 4.3] (http://snpeff.sourceforge.net/) Genetic variant annotation, effect prediction, VCF filtering and manipulation toolbox. `/home/sobczm/bin/snpEff`

[SplitsTree ver. 4] (http://www.splitstree.org/) Program for computing unrooted phylogenetic networks from molecular sequence data. `/home/sobczm/bin/splitstree4`

[STRUCTURE ver. 2.3.4] (http://pritchardlab.stanford.edu/structure.html) A package for using multi-locus genotype data to investigate population structure `/home/sobczm/bin/structure`

[structure Harvester ver. 0.6.93] (http://taylor0.biology.ucla.edu/structureHarvester) Downstream processing of STRUCTURE results to calculate Evanno’s Δk value and prepares input file for CLUMPP `/home/sobczm/bin/structureHarvester`

[CLUMPP ver. 1.1.2] (https://web.stanford.edu/group/rosenberglab/clumpp.html) Permutes the clusters output by independent runs of STRUCTURE, so that they match up as closely as possible. `/home/sobczm/bin/CLUMPP_Linux64.1.1.2`

[distruct ver. 1.1] (https://web.stanford.edu/group/rosenberglab/distruct.html) A program to graphically display results produced by STRUCTURE or by other similar programs. `/home/sobczm/bin/distruct1.1`

[Stampy ver. 1.0.29] (http://www.well.ox.ac.uk/project-stampy) Sensitive mapping of Illumina reads. `/home/sobczm/bin/stampy-1.0.29`

[Tracer ver. 1.6] (http://tree.bio.ed.ac.uk/software/tracer/) A program for analysing the trace files generated by Bayesian MCMC runs.
`/home/sobczm/bin/beast/Tracer_v1.6/bin`

[transAlign ver. 1.2] (https://www.uni-oldenburg.de/ibu/systematik-evolutionsbiologie/programme/) An open-source Perl script that aligns protein-coding DNA sequences via their amino-acid translations `/home/sobczm/bin/transalign` Dependency: `/home/sobczm/bin/clustalw1.83`

[Treemix ver. 1.12] (https://bitbucket.org/nygcresearch/treemix/wiki/Home): estimation of population trees with admixture.
`/home/sobczm/bin/treemix-1.12`

[Trinity ver 2.2] (https://trinityrnaseq.github.io/): RNA-Seq assembly  `/home/sobczm/bin/trinityrnaseq-2.2.0`

[vcflib] (https://github.com/vcflib/vcflib): a simple C++ library for parsing and manipulating VCF files, + many command-line utilities. `/home/sobczm/bin/vcflib/bin`

[VCFtools] (https://vcftools.github.io): another set of C++ and Perl libraries for analysing VCF files. `/home/sobczm/bin/vcftools/bin`

```
#SET PERL PATH
export PERL5LIB=/home/sobczm/bin/vcftools/share/perl/5.14.2
```

[PyVCF] (https://github.com/jamescasbon/PyVCF) A VCF v. 4.0 and 4.1 parser for Python. `/home/sobczm/bin/PyVCF/bin`

```
#SET PYTHON PATH
export PYTHONPATH="$PYTHONPATH:/home/sobczm/bin/PyVCF/lib/python2.7/site-packages/"
```
[USEARCH v. 9.0] (http://www.drive5.com/usearch/manual/) High-throughput search and clustering `/home/sobczm/bin/usearch`

[Weeder ver. 2.0] (http://159.149.160.51/modtools/) Discovery of transcription factor binding sites in a set of sequences from co-regulated genes `/home/sobczm/bin/weeder`

##Nanopore sequencing
[poretools ver. 0.6] (https://poretools.readthedocs.io/en/latest/index.html) A toolkit for working with nanopore sequencing data from Oxford Nanopore `/home/sobczm/bin/poretools/poretools` Usage: `python ./poretools`

[marginAlign] (https://github.com/benedictpaten/marginAlign) The marginAlign package can be used to align reads to a reference genome and call single nucleotide variations (SNVs). It is specifically tailored for Oxford Nanopore Reads. `/home/sobczm/bin/marginAlign`

[NanoOK] (https://documentation.tgac.ac.uk/display/NANOOK/NanoOK) Flexible, multi-reference software for pre- and post-alignment analysis of nanopore sequencing data, quality and error profiles `/home/sobczm/bin/NanoOK/bin`

[minion] (http://www.ebi.ac.uk/research/enright/software/kraken) A small utility program to infer or test the presence of 3' adapter sequence in sequencing data. `/home/sobczm/bin/minion`
