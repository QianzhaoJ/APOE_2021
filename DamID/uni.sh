#!/bin/bash
tar=/work2/liugh/wuqiong/05_Results/08_APOE/Dam_P9
raw=/work2/liugh/wuqiong/02_Data/APOE/Dam_P9/data

uni=$tar/04_unique
for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $uni/scripts/${sample}_uni.sh
echo '#BSUB -q c_liugh2' >> $uni/scripts/${sample}_uni.sh
echo "#BSUB -o $sample.out" >> $uni/scripts/${sample}_uni.sh
echo "#BSUB -o $sample.err" >> $uni/scripts/${sample}_uni.sh
echo "#BSUB -J uni_$sample" >> $uni/scripts/${sample}_uni.sh
echo '#BSUB -n 24' >> $uni/scripts/${sample}_uni.sh
echo '#BSUB -R "span[ptile=24]"' >> $uni/scripts/${sample}_uni.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Unique: $sample is Starting

mapping=$tar/03_mapping
uni=$tar/04_unique
picard=/work2/liugh/wuqiong/software/picard.jar
sam=$mapping/$sample/${sample}.sam

result=$uni/$sample
mkdir -p $result
log=$uni/logs

bam=$result/${sample}_no_chrY_chrM.bam
unique_bam=$result/${sample}_no_chrY_chrM_unique.bam
srt_bam=$result/${sample}_no_chrY_chrM.sort.bam
srt_rmdup_bam=$result/${sample}_no_chrY_chrM.sort.rmdup.bam
METRICS_FILE=$result/${sample}.picard.matrix
idxstats_txt=$result/${sample}.sort.rmdup.idxstats.txt
flagstat_txt=$result/${sample}.sort.rmdup.flagstat.txt

echo RMduplicates: $sample is Starting
awk '"'"'$3!="chrM"'"'"' $sam | samtools view -@ 24 -S -b -q 20 > $bam 2>$result/${sample}.log &&
awk '"'"'$3!="chrM"'"'"' $sam | grep -v 'XS:i:' | samtools view -@ 24 -bS -q 10 > $unique_bam 2>>$result/${sample}.log &&
samtools sort -@ 24 -l 9 $bam -o $srt_bam 2>>$result/${sample}.log &&
java -jar $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT=$srt_bam METRICS_FILE=$METRICS_FILE OUTPUT=$srt_rmdup_bam 2>>$result/${sample}.log &&
samtools index -@ 24 $srt_rmdup_bam 2>>$result/${sample}.log &&
samtools flagstat -@ 24 $srt_rmdup_bam > $flagstat_txt 2>>$result/${sample}.log &&
rm $sam

echo RMduplicates: $sample is Done'>> $uni/scripts/${sample}_uni.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$uni/run_uni.sh

