#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace11.blacklace

# Script to prepare rna for downstream applications.
# Will filter poor quality reads, perform trimming and
# remove illumina adapters.
# To be run from the project directory. Usage:
# rna_qc_fastq-mcf <RNASeq_F.fq> <RNASeq_R.fq>

#######  Step 1	 ########
# Initialise values	#
#########################


fastq_mcf=/home/armita/all_idris/idris/prog/ea-utils.1.1.2-537/fastq-mcf

WORK_DIR=$TMPDIR/rna_qc

F_IN=$1
R_IN=$2
OUTDIR=$PWD
ILLUMINA_ADAPTERS=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa

#######  Step 2	 ########
# 	unzip reads			#
#########################

mkdir -p "$WORK_DIR"/F "$WORK_DIR"/R
cd "$WORK_DIR"

F_FILE=$(echo $F_IN | rev | cut -d "/" -f1 | rev )
R_FILE=$(echo $R_IN | rev | cut -d "/" -f1 | rev )

F_OUT=$(echo "$F_FILE" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')
R_OUT=$(echo "$R_FILE" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')


#######  Step 4	 ########
# 	Quality trim		#
#########################

$fastq_mcf $ILLUMINA_ADAPTERS $F_FILE $R_FILE -o F/"$F_OUT" -o R/"$R_OUT" -C 1000000 -u -k 20 -t 0.01 -q 30 -p 5

gzip F/"$F_OUT"
gzip R/"$R_OUT"
cp -r F/"$F_OUT".gz $OUTDIR
cp -r R/"$R_OUT".gz $OUTDIR


#######  Step 8  ########
#       Cleanup         #
#########################

rm -r $WORK_DIR/
