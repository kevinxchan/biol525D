ref=/mnt/chan/Topic_4/ref/reference.fa
ngm=/home/biol525d/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm
fastq=/home/biol525d/Topic_4/fastq
proj=biol525d
bam=/mnt/chan/Topic_4/bam
ls $fastq | grep R1.fastq.gz | sed 's/.R1.fastq.gz//g' > $bam/samplelist.txt

while read name
do
	$ngm \
	-r $ref \
	-1 $fastq/${name}.R1.fastq.gz \
	-2 $fastq/${name}.R2.fastq.gz \
	-o $bam/${name}.ngm.sam \
	-t 2 \
	--rg-id ${name} \
	--rg-sm ${name} \
	--rg-pl illumina \
	--rg-pu $proj \
	--rg-lb ${name}_lib 

	samtools view -bh $bam/${name}.ngm.sam | samtools sort > $bam/${name}.ngm.bam 
	samtools index $bam/${name}.ngm.bam
done < $bam/samplelist.txt

