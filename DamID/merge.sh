#!/bin/bash
tar=/work2/liugh/wuqiong/05_Results/08_APOE/Dam_P9
raw=/work2/liugh/wuqiong/02_Data/APOE/Dam_P9/data

uni=$tar/04_unique

for sample in {Dam_APOE_P9,Dam_H1_P9,EMD_APOE_P9,EMD_H1_P9}
do

echo -e '#!/bin/bash' > $uni/scripts/${sample}_m.sh
echo '#BSUB -q c_liugh2' >> $uni/scripts/${sample}_m.sh
echo "#BSUB -o $sample.out" >> $uni/scripts/${sample}_m.sh
echo "#BSUB -o $sample.err" >> $uni/scripts/${sample}_m.sh
echo "#BSUB -J QC_$sample" >> $uni/scripts/${sample}_m.sh
echo '#BSUB -n 24' >> $uni/scripts/${sample}_m.sh
echo '#BSUB -R "span[ptile=24]"' >> $uni/scripts/${sample}_m.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Merge: $sample is Starting

uni=$tar/04_unique

result=$uni/$sample
mkdir -p $result
log=$uni/logs

bam1=$uni/${sample}_1/${sample}_1_no_chrY_chrM.sort.rmdup.bam
bam2=$uni/${sample}_2/${sample}_2_no_chrY_chrM.sort.rmdup.bam
bam3=$uni/${sample}_3/${sample}_3_no_chrY_chrM.sort.rmdup.bam

samtools merge -@ 24 $result/${sample}_merge.bam \
$bam1 \
$bam2 \
$bam3

samtools index -@ 24 $result/${sample}_merge.bam
samtools flagstat $result/${sample}_merge.bam > $result/${sample}_merge.flagstat

echo Merge: $sample is Done'>> $uni/scripts/${sample}_m.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *_m.sh`

do 
   bsub < $i &

done' >$uni/run_m.sh
