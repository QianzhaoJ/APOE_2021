#!/bin/bash
if [ $# -lt 2 ]
	then
	echo 1>&2 Usage: align_local_fq.sh input_file.fq.gz assembly_version codedir [model-matrix-file]
	exit 1
fi

# set some paths for executables
BOWTIE2='bowtie2'
CUTADAPT='cutadapt'
FASTX_REVCOM='fastx_reverse_complement'
BOWTIE2_INDEXES="/work2/liugh/wuqiong/03_Database/02_bowtie2_index/hg19"
SAMT='samtools'

# echo some version info to log:
echo 'using bowtie2 version:'
echo `${BOWTIE2} --version`
echo ''
echo 'using cutadapt version:'
echo `${CUTADAPT} --version`
echo ''
echo 'using fastx_reverse_complement version:'
echo `${FASTX_REVCOM} -h | grep FASTX`
echo ''

#set some parameters
IN_FQ=$1
#adaptor sequences used in DamID-seq
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_LONG_5="CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAG"
#adaptor sequences used by Illumina
#ILLUMINA_5="GCTCTTCCGATCT"
#make reverse complement of adapter sequences 
#(${FASTX_REVCOM} expects fasta input, "awk 'NR > 1'" means print everything beyond line 1)
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
ADPTR_LONG_3=`echo -e ">\n${ADPTR_LONG_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
#ILLUMINA_3=`echo -e ">\n${ILLUMINA_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #AGATCGGAAGAGC

fq_base=${IN_FQ} # save only file name
basef=${fq_base%.fastq.gz}
stats=${fq_base%.fastq.gz}/stats
mkdir ${fq_base%.fastq.gz} $stats



# set base filename
case $1 in
*.fq.gz) 
OUT_BAM=$basef/${IN_FQ%.fq.gz}_local.bam
CAT=zcat ;;
*.fastq.gz) 
OUT_BAM=$basef/${IN_FQ%.fastq.gz}_local.bam
CAT=zcat ;;
*.fastq) 
OUT_BAM=$basef/${IN_FQ%.fastq}_local.bam
CAT=cat ;;
*.fq) 
OUT_BAM=$basef/${IN_FQ%.fq}_local.bam
CAT=cat ;;
*) 
echo "inputfile (${IN_FQ}) should either be gzipped fastq"
echo "(*.fq.gz or fastq.gz) or fastq (*.fq or *.fastq)" 
exit 1 ;;
esac

# set additional output files based on base filename
BWT_STATS=${OUT_BAM%_local.bam}_local.bowtie_stats
CLIP_STATS=${OUT_BAM%_local.bam}_local.clip_stats

# set species and assembly (for read alignment)
case $2 in
dm3)
ASSEMBLY=dm3 ;;
hg18)
ASSEMBLY=hg18 ;;
hg19)
ASSEMBLY=hg19 ;;
mm9)
ASSEMBLY=mm9 ;;
*)
echo "assembly $2 not known, exiting"
exit 1 ;;
esac

################################################
# create temporary output files
################################################
TMP_BAM=`mktemp $basef/tmp_bam.XXXXXXXXXX`
TMP_FQ=`mktemp $basef/tmp_fq.XXXXXXXXXX`

# define function that checks exit code last command
CheckExit()
{
	if [ $1 -ne 0 ]; then
		echo "$2" 1>&2
			exit 1
			fi
}

# check whether outfiles exists already
if [ -f "${OUT_BAM}" -a -f "${BWT_STATS}" -a -f "${BAM_OUT}".bai ]; then
echo "${IN_FQ} is aligned already, exiting" 1>&2
exit 0
else
# remove if only some exist
rm -f "${OUT_BAM}" "${BWT_STATS}" "${BAM_OUT}".bai
fi

##################################################################################
# trim adapter sequences from reads, sort reads by precence of the GATC site into reads 
# Those reads where were found adapters are edge1, those reads where there are GATC at 5' or 3'
# ends are edge2. Others are inner
##################################################################################
# and write statistic to file 																									 #
##################################################################################

# Main trim reads. Processing cutadapt with the following parameters: 5' & 3'
# adapters encountered 2 or more times, the overlap of the adapter 30 for long
# (AdRt), 12 for short (AdRb) or more bases. Those reads with adapters go to
# trimmed[1-2].fastq, where there were no adapters write to file untrim_out.fastq 
${CUTADAPT} -g "${ADPTR_LONG_5}" -a "${ADPTR_LONG_3}" -e 0.01 -n 2 -m 20 --overlap 30 --match-read-wildcards --untrimmed-output $basef/untrimmed_pre_out.fastq ${IN_FQ} -o $basef/trimmed_1.fastq > $stats/clip_long_${fq_base%.fastq.gz}.stats

${CUTADAPT} -g "${ADPTR_SHORT_5}" -a "${ADPTR_SHORT_3}" -e 0.01 -n 2 -m 20 --overlap 12 --match-read-wildcards --untrimmed-output $basef/untrimmed.fastq $basef/untrimmed_pre_out.fastq -o $basef/trimmed_2.fastq > $stats/clip_short_${fq_base%.fastq.gz}.stats

# Combine adapters which had long or short adapters to trimmed.fastq
cat $basef/trimmed_1.fastq $basef/trimmed_2.fastq > $basef/trimmed.fastq


# Combine all found reads to one file
	cat $basef/trimmed.fastq $basef/untrimmed.fastq > $basef/interim_gatcs.fastq

# Remove reads with inner GATC's
	cat $basef/interim_gatcs.fastq | perl -e 'while($h=<>){$s=<>;$t=<>.<>; if($s!~/.+GATC.+/){print $h.$s.$t}}' > ${TMP_FQ}


# remove intermediate files
rm -R $stats $basef/*trimmed*.fastq $basef/interim_gatcs.fastq

# exit succesful
exit 0
