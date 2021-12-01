#!/bin/bash
raw=/work2/liugh/wuqiong/02_Data/APOE/ATAC_new
tar=/work2/liugh/wuqiong/05_Results/08_APOE/ATAC_new

index=/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19

################## 01 mkdir 
cd $tar
mkdir -p 01_fastqc 02_trim 03_mapping 04_unique 05_bamtobed 06_bamtobw 07_MACS2 08_deeptools 09_plot
for path in `ls $tar`
do
mkdir -p  $path/scripts
mkdir -p $path/logs
done

rm -r $tar/scripts/logs
rm -r $tar/scripts/scripts

################ 02 trim

trim=$tar/02_trim

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $trim/scripts/${sample}_trim.sh
echo '#BSUB -q c_liugh2' >> $trim/scripts/${sample}_trim.sh
echo "#BSUB -o $sample.log" >> $trim/scripts/${sample}_trim.sh
echo "#BSUB -e $sample.err" >> $trim/scripts/${sample}_trim.sh
echo "#BSUB -J QC_$sample" >> $trim/scripts/${sample}_trim.sh
echo '#BSUB -n 2' >> $trim/scripts/${sample}_trim.sh
echo '#BSUB -R "span[ptile=2]"' >> $trim/scripts/${sample}_trim.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Trimming: $sample is Starting

trim=$tar/02_trim
log=$trim/logs

fq1=$raw/$sample/${sample}_R1.fq.gz
fq2=$raw/$sample/${sample}_R2.fq.gz

result=$trim/$sample
mkdir $result

trim_galore --fastqc --path_to_cutadapt ~/software/anaconda2/bin/cutadapt --stringency 3 --paired --output_dir $result $fq1 $fq2 2>$log/${sample}.log

echo Trimming has been Done' >> $trim/scripts/${sample}_trim.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$trim/run_trim.sh

################ 03 mapping

map=$tar/03_mapping

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $map/scripts/${sample}_map.sh
echo '#BSUB -q c_liugh2' >> $map/scripts/${sample}_map.sh
echo "#BSUB -o $sample.log" >> $map/scripts/${sample}_map.sh
echo "#BSUB -e $sample.err" >> $map/scripts/${sample}_map.sh
echo "#BSUB -J Mapping_$sample" >> $map/scripts/${sample}_map.sh
echo '#BSUB -n 24' >> $map/scripts/${sample}_map.sh
echo '#BSUB -R "span[ptile=24]"' >> $map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw\\nindex=${index}\\n'

echo Mapping: $sample is Starting
trim=$tar/02_trim
mapping=$tar/03_mapping

clean1=$trim/$sample/${sample}*1.fq.gz
clean2=$trim/$sample/${sample}*2.fq.gz

result=$mapping/$sample
log=$mapping/logs

mkdir -p $result

bowtie2 -p 12 -x $index --no-mixed --no-discordant -t -1 $clean1 -2 $clean2 -S $result/${sample}.sam 2>$log/${sample}.log

echo Mapping: $sample is Done

echo Trimming has been Done' >> $map/scripts/${sample}_map.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$map/run_map.sh

################## 04 unique

uni=$tar/04_unique

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $uni/scripts/${sample}_uni.sh
echo '#BSUB -q c_liugh2' >> $uni/scripts/${sample}_uni.sh
echo "#BSUB -o $sample.log" >> $uni/scripts/${sample}_uni.sh
echo "#BSUB -e $sample.err" >> $uni/scripts/${sample}_uni.sh
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

bamraw=$result/${sample}.bam
bam=$result/${sample}_no_chrY.bam
unique_bam=$result/${sample}_no_chrY_chrM_unique.bam
srt_bam=$result/${sample}_no_chrY_chrM.sort.bam
srt_rmdup_bam=$result/${sample}_no_chrY_chrM.sort.rmdup.bam
METRICS_FILE=$result/${sample}.picard.matrix
idxstats_txt=$result/${sample}.sort.rmdup.idxstats.txt
flagstat_txt=$result/${sample}.sort.rmdup.flagstat.txt

echo RMduplicates: $sample is Starting
awk '"'"'$3!="chrM"'"'"' $sam | samtools view -@ 24 -S -b -q 10 > $bam 2>$result/${sample}.log &&
samtools sort -@ 24 -l 9 $bam -o $srt_bam 2>>$result/${sample}.log &&
java -jar $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT=$srt_bam METRICS_FILE=$METRICS_FILE OUTPUT=$srt_rmdup_bam 2>>$result/${sample}.log &&
samtools index -@ 24 $srt_rmdup_bam 2>>$result/${sample}.log &&
samtools flagstat -@ 24 $srt_rmdup_bam > $flagstat_txt 2>>$result/${sample}.log &&
#rm $sam
echo RMduplicates: $sample is Done'>> $uni/scripts/${sample}_uni.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$uni/run_uni.sh

################### 05_bam_to_bed

bed=$tar/05_bamtobed
for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $bed/scripts/${sample}_bed.sh
echo '#BSUB -q c_liugh2' >> $bed/scripts/${sample}_bed.sh
echo "#BSUB -o $sample.log" >> $bed/scripts/${sample}_bed.sh
echo "#BSUB -e $sample.err" >> $bed/scripts/${sample}_bed.sh
echo "#BSUB -J BED_$sample" >> $bed/scripts/${sample}_bed.sh
echo '#BSUB -n 24' >> $bed/scripts/${sample}_bed.sh
echo '#BSUB -R "span[ptile=24]"' >> $bed/scripts/${sample}_bed.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Bam to Bed: $sample is Starting

uni=$tar/04_unique/$sample
srt_rmdup_bam=$uni/${sample}_no_chrY_chrM.sort.rmdup.bam
bl=/work2/liugh/wuqiong/03_Database/anno_bed/kundaje_EncodeHg19ConsensusSignalArtifactRegions.bed

bed=$tar/05_bamtobed/$sample
log=$tar/05_bamtobed/logs
mkdir -p $bed

srt_rmdup_bed=$bed/${sample}_no_chrM_chrY_unique.sort.rmdup.bed
interdect_BL_bed=$bed/${sample}_no_chrM_chrY_unique.sort.rmdup.interdect_BL.bed
remove_BL_bed=$bed/${sample}_no_chrM_chrY_unique.sort.rmdup.remove_BL.bed

bamToBed -i $srt_rmdup_bam > $srt_rmdup_bed 2>$log/${sample}.log &&
bedtools  intersect -u -a $srt_rmdup_bed -b $bl > $interdect_BL_bed 2>>$log/${sample}.log &&
bedtools  intersect -v -a $srt_rmdup_bed -b $bl > $remove_BL_bed 2>>$log/${sample}.log

echo Bamtobed: $sample is Done'>> $bed/scripts/${sample}_bed.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$bed/run_bed.sh

############################## MACS2

macs=$tar/07_MACS2

for sample in `ls $raw`
do
echo -e '#!/bin/bash' > $macs/scripts/${sample}_macs.sh
echo '#BSUB -q c_liugh2' >> $macs/scripts/${sample}_macs.sh
echo "#BSUB -o $sample.log" >> $macs/scripts/${sample}_macs.sh
echo "#BSUB -e $sample.err" >> $macs/scripts/${sample}_macs.sh
echo "#BSUB -J MACS_$sample" >> $macs/scripts/${sample}_macs.sh
echo '#BSUB -n 4' >> $macs/scripts/${sample}_macs.sh
echo '#BSUB -R "span[ptile=4]"' >> $macs/scripts/${sample}_macs.sh
echo -e tar=$tar\\nraw=$raw\\n >>$macs/scripts/${sample}_macs.sh

echo sample=$sample'

perl=/work2/liugh/liuzunpeng/04_Softwares/perl/local/perl/bin/perl
python=/work2/liugh/liuzunpeng/04_Softwares/python/local/python2/bin/python
macs2=/work2/liugh/liuzunpeng/04_Softwares/python/local/python2/bin/macs2

annotatePeaks=/work2/liugh/liuzunpeng/04_Softwares/homer/bin/annotatePeaks.pl
macs_stat=/work2/liugh/liuzunpeng/04_Softwares/ATAC_seq_src/macs_stat.pl

uni=$tar/04_unique
peak=$tar/07_MACS2/$sample

bed=$tar/05_bamtobed
remove_BL_bed=$bed/$sample/${sample}_no_chrM_chrY_unique.sort.rmdup.remove_BL.bed
macs_stat=/data2/liuzunpeng/03_Database/ATAC/macs_stat.pl

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
$perl $macs_stat $peaks_bed > $peak_stat 

anno=$peak/anno 
mkdir -p $anno 
echo Homer annotion

$annotatePeaks $peaks_bed hg19 -genomeOntology $anno -annStats $anno/${sample}_stat.txt >$anno/${sample}_peak_anno.txt 2>>$log/${sample}.log
echo annotion finished '>>$macs/scripts/${sample}_macs.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$macs/run_macs2.sh

###################  deeptools convert
deep=$tar/08_deeptools
for sample in `ls $raw`

do

echo -e '#!/bin/bash' > $deep/scripts/${sample}_deep.sh
echo '#BSUB -q c_liugh2' >> $deep/scripts/${sample}_deep.sh
echo "#BSUB -o $sample.log" >> $deep/scripts/${sample}_deep.sh
echo "#BSUB -e $sample.err" >> $deep/scripts/${sample}_deep.sh
echo "#BSUB -J DEEP_$sample" >> $deep/scripts/${sample}_deep.sh
echo '#BSUB -n 24' >> $deep/scripts/${sample}_deep.sh
echo '#BSUB -R "span[ptile=24]"' >> $deep/scripts/${sample}_deep.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

uni=$tar/04_unique
deep=$tar/08_deeptools
echo Deeptools: $sample is Starting
srt_rmdup_bam=$uni/$sample/${sample}_extract.rmdup.bam

###################
bin=2000
result=$deep/${sample}/bin${bin}
mkdir -p $result
bw=$result/${sample}_bin${bin}.bw
bg=$result/${sample}_bin${bin}.bedGraph
log=$deep/logs/${sample}_bin${bin}.log
echo binSize=$bin  normalization= RPKM  extendReads= no  Application= replication
bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrX chrM chrY 2>$log &&
echo binSize=$bin is Done

echo Deeptools: $sample is Done'>> $deep/scripts/${sample}_deep.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *deep.sh`

do 
   bsub < $i &

done' >$deep/run_deep.sh

########################### final merge

