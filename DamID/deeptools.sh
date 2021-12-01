#!/bin/bash
tar=/work2/liugh/wuqiong/05_Results/08_APOE/Dam
raw=/work2/liugh/wuqiong/02_Data/APOE/Dam/data

deep=$tar/05_bamCoverage
for sample in {Dam_H9_P10_1,Dam_H9_P10_2,Dam_H9_P10_3,Dam_ZK3_P10_1,Dam_ZK3_P10_2,Dam_ZK3_P10_3,EMD_H9_P10_1,EMD_H9_P10_2,EMD_H9_P10_3,EMD_ZK3_P10_1,EMD_ZK3_P10_2,EMD_ZK3_P10_3}

do

echo -e '#!/bin/bash' > $deep/scripts/${sample}_deep_s.sh
echo '#BSUB -q c_liugh2' >> $deep/scripts/${sample}_deep_s.sh
echo "#BSUB -o $sample.out" >> $deep/scripts/${sample}_deep_s.sh
echo "#BSUB -o $sample.err" >> $deep/scripts/${sample}_deep_s.sh
echo "#BSUB -J DEEP_$sample" >> $deep/scripts/${sample}_deep_s.sh
echo '#BSUB -n 24' >> $deep/scripts/${sample}_deep_s.sh
echo '#BSUB -R "span[ptile=24]"' >> $deep/scripts/${sample}_deep_s.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

uni=$tar/04_unique
deep=$tar/05_bamCoverage
echo Deeptools: $sample is Starting
srt_rmdup_bam=$uni/$sample/${sample}_no_chrY_chrM.sort.rmdup.bam
#########################
bin=10
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=10  normalization= RPKM   Application= Track 
#bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
#bamCoverage -p 24 -b $srt_rmdup_bam --outFileFormat bedgraph -o $bg --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done
########################
bin=100
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads=NO   Application= TSS
#bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
#bamCoverage -p 24 -b $srt_rmdup_bam --outFileFormat bedgraph -o $bg --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done
#######################
bin=2000
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads= no  Application= replication
#bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
#bamCoverage -p 24 -b $srt_rmdup_bam --outFileFormat bedgraph -o $bg --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done

#######################
bin=25000
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads= no  Application= replication
bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
bamCoverage -p 24 -b $srt_rmdup_bam --outFileFormat bedgraph -o $bg --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done




echo RMduplicates: $sample is Done'>> $deep/scripts/${sample}_deep_s.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *deep_s.sh`

do 
   bsub < $i &

done' >$deep/run_deep.sh

