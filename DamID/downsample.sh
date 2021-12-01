#!/bin/bash
#BSUB -q c_liugh2
#BSUB -o extract.out
#BSUB -o extract.err
#BSUB -J extract
#BSUB -n 24
#BSUB -R "span[ptile=24]"

function SubSample {

## Calculate the sampling factor based on the intended number of reads:
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]
  then 
  echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
fi

/work2/liugh/wuqiong/04_Softwares/sambamba-0.7.0-linux-static view -s $FACTOR -f bam -l 9 -t 20 $1

}

path=/work2/liugh/wuqiong/05_Results/08_APOE/Dam_P9/04_unique

for sample in {Dam_APOE_P9,Dam_H1_P9,EMD_APOE_P9,EMD_H1_P9}

do

inbam=$path/$sample/${sample}_merge.bam
n=123000000
outbam=$path/$sample/${sample}_extract.bam

SubSample $inbam $n > $outbam
/work2/liugh/wuqiong/04_Softwares/sambamba-0.7.0-linux-static flagstat -t 24 $outbam
samtools index $outbam
done

