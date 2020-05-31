# script to simply control quality, remove adaptors and align and filter.
# to modify as trimming can be optionnal.


if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <path/to/read1.fq.gz> <path/to/read2.fq.gz> <path/to/genome> <adaptor> <output prefix> "
    echo ""
    exit -1
fi


# Custom parameters
THREADS=4

# CMD parameters if using bowtie
FQ1=${1} #gz file
FQ2=${2} #gz file
GEN=${3}
ADAPTOR=${4}
OUTP=${5}


# pre Fastqc
mkdir -p ${OUTP}_prefastqc/ && fastqc -o ${OUTP}_prefastqc/ ${FQ1} && fastqc -o ${OUTP}_prefastqc/ ${FQ2}


####################################################
# Adapter trimming

 # -q <int> to trim low-quality ends from reads
 # -m <int> to discard reads shorter than <int>
 # -e <int> minimum error rate
 # -O <int> minimum overlap
 # -a 3' adapter for R1 (-g for 5')
 # -A 3' adapter for R2 (-G for 5')
cutadapt -q 10 -m 15 -e 0.10 -O 5 -a ${ADAPTOR} -A ${ADAPTOR} -o ${OUTP}.1.fq.gz -p ${OUTP}.2.fq.gz ${FQ1} ${FQ2} | gzip -c > ${OUTP}.cutadapt.log.gz

###################################################

# post Fastqc
mkdir -p ${OUTP}_postfastqc/ && fastqc -o ${OUTP}_postfastqc/ ${OUTP}.1.fq.gz && fastqc -o ${OUTP}_postfastqc/ ${OUTP}.2.fq.gz

##################################################

# Bowtie2
 # --very-sensitive  : -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
bowtie2 --very-sensitive -x ${GEN} -1 ${OUTP}.1.fq.gz -2 ${OUTP}.2.fq.gz -S ${OUTP}.sam

# SAM TO BAM
samtools view -S -b ${OUTP}.sam > ${OUTP}.bam

#adjustments // remove mitochondrial Chromosome
samtools view -h ${OUTP}.bam |  grep -v chrM  |  samtools view -b -  >  ${OUTP}_filtered.bam

#non unique alignements
samtools view -b  -q 10  ${OUTP}_filtered.bam  >  ${OUTP}_filtered.bam


###### MERGE BAM FILES FOR REPLICATES AND CHECK OVERLAPPING ######


# BAM TO SAM TO BED
samtools view -h  ${OUTP}_filtered.sam  |  SAMtoBED  -i -  -o ${OUTP}.bed  -x  -v

# MACS2
macs2 callpeak  -t ${OUTP}.bed  -f BEDPE  -n macs2_res  -g ce  --keep-dup all
