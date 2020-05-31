# Preprocessing RNA-seq analysis pipeline for paired ends fastq reads
# using FastQC, STAR and featureCounts

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <path/to/read1.fq.gz> <path/to/read2.fq.gz> <path/to/genome> <output_prefix_> "
    echo ""
    exit -1
fi


# Custom parameters
STHREADS=8 # STAR threads
FTHREADS=5 # featureCounts threads


# Arguments association
FQ1=${1} # read 1 - gz file
FQ2=${2} #read 2 - gz file
GEN=${3} # path to genome
OUTP=${4} #output prefix


# FastQC
mkdir -p ${OUTP}_prefastqc/ && fastqc -o ${OUTP}_prefastqc/ ${FQ1} && fastqc -o ${OUTP}_prefastqc/ ${FQ2}

##################################################################

# STAR mapping
STAR —runThreadN ${STHREADS} —genomeDir ${3} \
—readFilesIn ${1} ${2} \
—outFileNamePrefix ${4} \
—sjdbGTFfile ${3}/genes.gtf \
—outSAMstrandField intronMotif \
—outSAMtype BAM SortedByCoordinate \
—outReadsUnmapped Fastx;

#################################################################

#featureCounts
featureCounts -T ${FTHREADS} -t exon -g gene_id -a ${3}/genes.gtf -o ${4}counts.txt ${4}SortedByCoordinate.sam
