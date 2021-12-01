#!/bin/bash

tar=/work2/liugh/wuqiong/05_Results/08_APOE/Dam_P9
raw=/work2/liugh/wuqiong/02_Data/APOE/Dam_P9/data

map=$tar/03_mapping

for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $map/scripts/${sample}_map.sh
echo '#BSUB -q c_liugh2' >> $map/scripts/${sample}_map.sh
echo "#BSUB -o $sample.out" >> $map/scripts/${sample}_map.sh
echo "#BSUB -o $sample.err" >> $map/scripts/${sample}_map.sh
echo "#BSUB -J MAP_$sample" >> $map/scripts/${sample}_map.sh
echo '#BSUB -n 24' >> $map/scripts/${sample}_map.sh
echo '#BSUB -R "span[ptile=24]"' >> $map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Mapping: $sample is Starting
trim=$tar/02_trim
mapping=$tar/03_mapping

clean1=$trim/$sample/alignedReads/${sample}_1.p.fastq.gz
clean2=$trim/$sample/alignedReads/${sample}_2.p.fastq.gz

index=/work2/liugh/liuzunpeng/03_database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19

result=$mapping/$sample
log=$mapping/logs

mkdir -p $result

bowtie2 -p 24 -x $index --no-mixed --no-discordant -t -1 $clean1 -2 $clean2 -S $result/${sample}.sam 2>$log/${sample}.log

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

