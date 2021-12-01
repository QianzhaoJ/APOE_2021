#!/bin/bash
tar=/work2/liugh/wuqiong/05_Results/08_APOE/Dam
raw=/work2/liugh/wuqiong/02_Data/APOE/Dam/data

deep=$tar/06_bamCompare
for sample in {H1_P3,APOE_P3}

do

echo -e '#!/bin/bash' > $deep/scripts/${sample}_deep.sh
echo '#BSUB -q c_liugh2' >> $deep/scripts/${sample}_deep.sh
echo "#BSUB -o $sample.out" >> $deep/scripts/${sample}_deep.sh
echo "#BSUB -o $sample.err" >> $deep/scripts/${sample}_deep.sh
echo "#BSUB -J compare_$sample" >> $deep/scripts/${sample}_deep.sh
echo '#BSUB -n 24' >> $deep/scripts/${sample}_deep.sh
echo '#BSUB -R "span[ptile=24]"' >> $deep/scripts/${sample}_deep.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

uni=$tar/04_unique
deep=$tar/06_bamCompare
echo Deeptools: $sample is Starting
bam1=/work2/liugh/wuqiong/05_Results/08_APOE/Dam/07_HMM/Bams/EMD_${sample}_extract.bam
bam2=/work2/liugh/wuqiong/05_Results/08_APOE/Dam/07_HMM/Bams/Dam_${sample}_extract.bam
#########################
bin=10
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_EMDvsDam_bin${bin}.bw
bg=$result/${sample}_EMDvsDam_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=10  normalization= RPKM   Application= Track 
#bamCompare -p 24 -b1 $bam1 -b2 $bam2 -o $bw --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
#bamCompare -p 24 -b1 $bam1 -b2 $bam2 --outFileFormat bedgraph -o $bg --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done
########################
bin=100
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_EMDvsDam_bin${bin}.bw
bg=$result/${sample}_EMDvsDam_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads=NO   Application= TSS
#bamCompare -p 24 -b1 $bam1 -b2 $bam2 -o $bw  --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
#bamCompare -p 24 -b1 $bam1 -b2 $bam2 --outFileFormat bedgraph -o $bg --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done
#######################
bin=2000
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_EMDvsDam_bin${bin}.bw
bg=$result/${sample}_EMDvsDam_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads= no  Application= replication
#bamCompare -p 24 -b1 $bam1 -b2 $bam2 -o $bw --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
#bamCompare -p 24 -b1 $bam1 -b2 $bam2 --outFileFormat bedgraph -o $bg --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done

#######################
bin=25000

result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_EMDvsDam_bin${bin}.bw
bg=$result/${sample}_EMDvsDam_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads= no  Application= replication
bamCompare -p 24 -b1 $bam1 -b2 $bam2 -o $bw --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>$log &&
bamCompare -p 24 -b1 $bam1 -b2 $bam2 --outFileFormat bedgraph -o $bg --operation log2 --binSize $bin --ignoreForNormalization chrM chrY 2>>$log
echo binSize=$bin is Done


echo RMduplicates: $sample is Done'>> $deep/scripts/${sample}_deep.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *deep.sh`

do 
   bsub < $i &

done' >$deep/run_deep.sh

