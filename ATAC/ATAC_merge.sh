#!/bin/bash


################## 01 merge
tar=/work2/liugh/wuqiong/05_Results/08_APOE/ATAC_new
uni=$tar/04_unique

for sample in {WT,APOE}
do

echo -e '#!/bin/bash' > $uni/scripts/${sample}_m.sh
echo '#BSUB -q c_liugh2' >> $uni/scripts/${sample}_m.sh
echo "#BSUB -o $sample.out" >> $uni/scripts/${sample}_m.sh
echo "#BSUB -o $sample.err" >> $uni/scripts/${sample}_m.sh
echo "#BSUB -J QC_$sample" >> $uni/scripts/${sample}_m.sh
echo '#BSUB -n 24' >> $uni/scripts/${sample}_m.sh
echo '#BSUB -R "span[ptile=24]"' >> $uni/scripts/${sample}_m.sh
echo -e sample=${sample}\\ntar=$tar\\n'

echo Merge: $sample is Starting

uni=$tar/04_unique

result=$uni/$sample
mkdir -p $result
log=$uni/logs

bam1=$uni/${sample}_1/${sample}_1_no_chrY_chrM.sort.rmdup.bam
bam2=$uni/${sample}_2/${sample}_2_no_chrY_chrM.sort.rmdup.bam
bam3=$uni/${sample}_3/${sample}_3_no_chrY_chrM.sort.rmdup.bam
bam4=$uni/${sample}_4/${sample}_4_no_chrY_chrM.sort.rmdup.bam

samtools merge -@ 24 $result/${sample}_merge.bam \
$bam1 \
$bam2 \
$bam3 \
$bam4 

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

################################# 02 bamtobed

uni=$tar/05_bamtobed
for sample in {WT,APOE}
do

echo -e '#!/bin/bash' > $uni/scripts/${sample}_m.sh
echo '#BSUB -q c_liugh2' >> $uni/scripts/${sample}_m.sh
echo "#BSUB -o $sample.out" >> $uni/scripts/${sample}_m.sh
echo "#BSUB -o $sample.err" >> $uni/scripts/${sample}_m.sh
echo "#BSUB -J uni_$sample" >> $uni/scripts/${sample}_m.sh
echo '#BSUB -n 2' >> $uni/scripts/${sample}_m.sh
echo '#BSUB -R "span[ptile=2]"' >> $uni/scripts/${sample}_m.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Bam to Bed: $sample is Starting

uni=$tar/04_unique/${sample}
srt_rmdup_bam=$uni/${sample}_extract.bam
bl=/work2/liugh/wuqiong/03_Database/anno_bed/kundaje_EncodeHg19ConsensusSignalArtifactRegions.bed

bed=$tar/05_bamtobed/${sample}
log=$tar/05_bamtobed/logs
mkdir -p $bed

srt_rmdup_bed=$bed/${sample}.merge.sort.rmdup.bed
interdect_BL_bed=$bed/${sample}.merge.interdect_BL.bed
remove_BL_bed=$bed/${sample}.merge.remove_BL.bed

bamToBed -i $srt_rmdup_bam > $srt_rmdup_bed 2>$log/${sample}.log &&
bedtools  intersect -u -a $srt_rmdup_bed -b $bl > $interdect_BL_bed 2>>$log/${sample}.log &&
bedtools  intersect -v -a $srt_rmdup_bed -b $bl > $remove_BL_bed 2>>$log/${sample}.log

echo RMduplicates: $sample is Done'>> $uni/scripts/${sample}_m.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *m.sh`

do 
   bsub < $i &

done' >$uni/run_m.sh

########################## 03 MACS2
macs=$tar/07_MACS2

for sample in {WT,APOE}
do
echo -e '#!/bin/bash' > $macs/scripts/${sample}_m.sh
echo '#BSUB -q c_liugh2' >> $macs/scripts/${sample}_m.sh
echo "#BSUB -o $sample.out" >> $macs/scripts/${sample}_m.sh
echo "#BSUB -o $sample.err" >> $macs/scripts/${sample}_m.sh
echo "#BSUB -J MACS_$sample" >> $macs/scripts/${sample}_m.sh
echo '#BSUB -n 4' >> $macs/scripts/${sample}_m.sh
echo '#BSUB -R "span[ptile=4]"' >> $macs/scripts/${sample}_m.sh
echo -e tar=$tar\\n >>$macs/scripts/${sample}_m.sh

echo sample=$sample'

python=/work2/liugh/liuzunpeng/04_Softwares/python/local/python2/bin/python
macs2=/work2/liugh/liuzunpeng/04_Softwares/python/local/python2/bin/macs2
annotatePeaks=/work2/liugh/liuzunpeng/04_Softwares/homer/bin/annotatePeaks.pl
macs_stat=/work2/liugh/liuzunpeng/04_Softwares/ATAC_seq_src/macs_stat.pl

uni=$tar/04_unique
peak=$tar/07_MACS2/$sample

bed=$tar/05_bamtobed
remove_BL_bed=$bed/$sample/${sample}.merge.sort.rmdup.bed

log=$tar/07_MACS2/logs

mkdir -p $peak
echo MACS2 call peak

peaks_xls=$peak/${sample}_peaks.xls
peaks_bed=$peak/${sample}.macs2.peaks.bed
peak_stat=$peak/${sample}.MACS2_peak_stat.txt

$python $macs2 callpeak -t $remove_BL_bed -f BED -g hs -B --nomodel --shift 0 --extsize 250 --call-summits -q 0.01 -n $sample --outdir $peak 2>$log/${sample}.log
echo Calling peak: $sample is Done &&
echo Stat: $sample is Starting &&
grep '"'"'^chr\S'"'"' $peaks_xls | awk '"'"'{print $1"\t"$2"\t"$3"\t"$10"\t"$8"\t""+"}'"'"' >$peaks_bed &&

anno=$peak/anno 
mkdir -p $anno 
echo Homer annotion
/work2/liugh/liuzunpeng/04_Softwares/homer/bin/annotatePeaks.pl $peaks_bed hg19 -genomeOntology $anno -annStats $anno/${sample}_stat.txt >$anno/${sample}_peak_anno.txt 2>>$log/${sample}.log
echo annotion finished '>>$macs/scripts/${sample}_m.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *_m.sh`

do 
   bsub < $i &

done' >$macs/run_m.sh
############################# 04 deeptools

deep=$tar/08_deeptools

for sample in {WT,APOE}
do

echo -e '#!/bin/bash' > $deep/scripts/${sample}_m.sh
echo '#BSUB -q c_liugh2' >> $deep/scripts/${sample}_m.sh
echo "#BSUB -o $sample.out" >> $deep/scripts/${sample}_m.sh
echo "#BSUB -o $sample.err" >> $deep/scripts/${sample}_m.sh
echo "#BSUB -J DEEP_$sample" >> $deep/scripts/${sample}_m.sh
echo '#BSUB -n 24' >> $deep/scripts/${sample}_m.sh
echo '#BSUB -R "span[ptile=24]"' >> $deep/scripts/${sample}_m.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

uni=$tar/04_unique
deep=$tar/08_deeptools
echo Deeptools: $sample is Starting
srt_rmdup_bam=$uni/${sample}_extract.bam
#########################
bin=10
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=10  normalization= RPKM  extendReads=250   Application= Track 
bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --extendReads 250 --ignoreForNormalization chrX chrM chrY 2>$log &&
echo binSize=$bin is Done
########################
bin=100
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}_extend_250.bw
bg=$result/${sample}_bin${bin}_extend_250.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads=250   Application= TSS
bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --extendReads 250 --ignoreForNormalization chrX chrM chrY 2>$log &&
echo binSize=$bin is Done
#######################
bin=2000
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads= no  Application= replication
bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrX chrM chrY 2>$log &&
echo binSize=$bin is Done

echo RMduplicates: $sample is Done'>> $deep/scripts/${sample}_m.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *_m.sh`

do 
   bsub < $i &

done' >$deep/run_m.sh

