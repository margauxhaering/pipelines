# merging 3 replicates and filtering to perform peak calling with HMMRATAC

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <path/to/rep1.sorted.bam> <path/to/rep2.sorted.bam> <path/to/rep3.sorted.bam> <output prefix> "
    echo ""
    exit -1
fi



# CMD parameters
REP1=${1} #bam file
REP2=${2} #bam file
REP3=${3} #bam file
OUPT=${4} #output prefix



#merge replicates prior to hmmratac :
samtools merge –n ${OUTP}_merged.bam ${REP1}.sorted.bam ${REP2}.sorted.bam ${REP3}.sorted.bam

# if uncertain of the merge to be sorted :
samtools sort -@ 4 -o ${OUTP}_merged.sorted.bam ${OUTP}_merged.bam
#########################################################


#index file for the bam : (named atac.forHMMRATAC.bam.bai)
samtools index ${OUTP}_merged.sorted.bam


#########################################################
# genome info file for chr size
mysql –-user=genome --host=genome-mysql.cse.ucsc.edu –A –e \ “select chrom, size from dm6.chromInfo” > genome.info

# if the sam file has header, bam is generated (step 1 and 2) : can extract info from bam file :
# samtools view -H ${OUTP}_merged.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info
########################################################

#hmmratac
java -jar HMMRATAC_V1.2_exe.jar -b ${OUTP}_merged.sorted.bam-i ${OUTP}_merged.sorted.bam.bai -g genome.info <options>
