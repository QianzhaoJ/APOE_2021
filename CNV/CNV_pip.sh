#!/bin/bash

tar=xxxx
raw=xxxx

cd $tar
mkdir -p 01_fastqc 02_trim 03_mapping 04_unique 05_count
for path in `ls $tar`
do
mkdir -p  $path/scripts
mkdir -p $path/logs
done

rm -r $tar/scripts/logs
rm -r $tar/scripts/scripts

######################### 02 trim

trim=$tar/02_trim

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $trim/scripts/${sample}_trim.sh
echo '#BSUB -q c_liugh2' >> $trim/scripts/${sample}_trim.sh
echo "#BSUB -o $sample.out" >> $trim/scripts/${sample}_trim.sh
echo "#BSUB -o $sample.err" >> $trim/scripts/${sample}_trim.sh
echo "#BSUB -J Trim_$sample" >> $trim/scripts/${sample}_trim.sh
echo '#BSUB -n 2' >> $trim/scripts/${sample}_trim.sh
echo '#BSUB -R "span[ptile=2]"' >> $trim/scripts/${sample}_trim.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Trimming: $sample is Starting

trim=$tar/02_trim
log=$trim/logs

fq1=$raw/$sample/${sample}_1.fq.gz
fq2=$raw/$sample/${sample}_2.fq.gz

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

############################################ 03 mapping

map=$tar/03_mapping

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $map/scripts/${sample}_map.sh
echo '#BSUB -q c_liugh2' >> $map/scripts/${sample}_map.sh
echo "#BSUB -o $sample.out" >> $map/scripts/${sample}_map.sh
echo "#BSUB -o $sample.err" >> $map/scripts/${sample}_map.sh
echo "#BSUB -J MAP_$sample" >> $map/scripts/${sample}_map.sh
echo '#BSUB -n 12' >> $map/scripts/${sample}_map.sh
echo '#BSUB -R "span[ptile=12]"' >> $map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Mapping: $sample is Starting
trim=$tar/02_trim
mapping=$tar/03_mapping

clean1=$trim/$sample/${sample}*1.fq.gz
clean2=$trim/$sample/${sample}*2.fq.gz

index=/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19

result=$mapping/$sample
log=$mapping/logs

mkdir -p $result

bowtie2 -p 24 -x $index --no-mixed --no-discordant -t -1 $clean1 -2 $clean2 -S $result/${sample}.sam 2>$log/${sample}.log

echo Mapping: $sample is Done' >> $map/scripts/${sample}_map.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$map/run_map.sh

############################################ 04 unique

uni=$tar/04_unique
mkdir -p $uni
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
awk '"'"'$3!="chrM" && $3!="chrY"'"'"' $sam | samtools view -@ 24 -S -b -q 10 > $bam 2>$result/${sample}.log &&
samtools sort -@ 24 -l 9 $bam -o $srt_bam 2>>$result/${sample}.log &&
java -jar $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT=$srt_bam METRICS_FILE=$METRICS_FILE OUTPUT=$srt_rmdup_bam 2>>$result/${sample}.log &&
samtools index -@ 24 $srt_rmdup_bam 2>>$result/${sample}.log &&
samtools flagstat -@ 24 $srt_rmdup_bam > $flagstat_txt 2>>$result/${sample}.log &&

echo RMduplicates: $sample is Done'>> $uni/scripts/${sample}_uni.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$uni/run_uni.sh

######################### 05 count

count=$tar/05_count
mkdir -p $count/logs
mkdir -p $count/scripts

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $count/scripts/${sample}_count.sh
echo '#BSUB -q c_liugh2' >> $count/scripts/${sample}_count.sh
echo "#BSUB -o $sample.out" >> $count/scripts/${sample}_count.sh
echo "#BSUB -o $sample.err" >> $count/scripts/${sample}_count.sh
echo "#BSUB -J Count_$sample" >> $count/scripts/${sample}_count.sh
echo '#BSUB -n 12' >> $count/scripts/${sample}_count.sh
echo '#BSUB -R "span[ptile=12]"' >> $count/scripts/${sample}_count.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

uni=$tar/04_unique
count=$tar/05_count/$sample
bam=$uni/$sample/${sample}_no_chrY_chrM.sort.rmdup.bam
log=$tar/05_count/logs

readCounter=/work2/liugh/liuzunpeng/04_Softwares/hmmcopy_utils/hmmcopy_utils/bin/readCounter

echo ReadCount: $sample is Starting
mkdir -p $count
wig=$count/${sample}_read.wig
$readCounter -w 500000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX $bam > $wig
echo ReadCount: $sample is Done ' >>$count/scripts/${sample}_count.sh
done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`

do 
   bsub < $i &

done' >$count/run_count.sh


