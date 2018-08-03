#!usr/bin/env bash

PERL=$(which perl)
R=$(which R)
JAVA=$(which java)
PYTHON=$(which python)
PROGRAMHOME=/home/ubuntu/bin
PRINSEQLITE=$PROGRAMHOME/prinsite-lite.pl
BWA=$PROGRAMHOME/bwa
SAMTOOLS=$PROGRAMHOME/samtools
GATK=$PROGRAMHOME/gatk
PICARD=$PROGRAMHOME/picard
PLINK=$PROGRAMHOME/plink
STRUCTURE=$PROGRAMHOME/structure.py
CHOOSEK=$PROGRAMHOME/chooseK.py
data_dir=$PWD/data
ref_genome=$PWD/ref
PROJECT=biol525d

usage() {
echo "
Written by Kevin Chan

Description: Final coding project for BIOL 525D.

Completes the following steps:
1. Trims fastq files for quality.
2. Aligns sequences to a reference genome.
3. Calls SNPs.
4. Filters vcf output.
5. Runs a principal component analysis (PCA).
6. Outputs a figure in PDF format of the PCA.

Assumptions:
1. All programs and libraries are in /home/ubuntu/bin.
2. All starting fastq files are in the same directory as this script, in a directory called 'data'.
3. A reference genome is in 'ref/reference.fa' in this directory.
4. All programs necessary are installed (obviously).
5. Rscripts used in this script are also in this directory (pca.R)
"
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi


echo starting pipeline...

echo ==============
echo fetching data
echo ==============

ls $data_dir/ | sed "s/.R[0-9].fastq.gz//g" | uniq > sample_names_list.txt

echo ================================================
echo quality trimming and filtering data with prinseq
echo ================================================

mkdir -p prinseqlite_log
mkdir -p data_cleaned

while read name; do
	$PERL $PRINSEQLITE \
	-fastq $name.R1.fastq.gz \
	-fastq2 $name.R2.fastq.gz \
	-log prinseqlite_log/log_$name \
	-out_good data_cleaned/$name.filtered \
	-min_len 70 \
	--trim_qual_left 10 \
	-trim_qual_right 10
done < sample_names_list.txt

echo ========================
echo aligned reads using BWA
echo ========================

mkdir -p bam

while read name; do
	r=$(echo -e "@RG\tID:Sample_$name\tSM:Sample_$name\tPL:illumina\tPU:$PROJECT\tLB:Sample_$name_lib")
	$BWA mem \
	$ref_genome/reference.fa \
	$data_cleaned/$name.filtered_1.fq.gz \
	$data_cleaned/$name.filtered_2.fq.gz \
	-t 2 \
	-R $r |\
	$SAMTOOLS view -bh \
	$SAMTOOLS sort | \
	$SAMTOOLS index > bam/$name.bwa.bam
done < sample_names_list.txt

echo ==================================
echo marking PCR duplicates with picard
echo ==================================

mkdir -p mark_duplicates_log

while read name; do
	$JAVA -jar $PICARD MarkDuplicates \
	I=bam/$name.bwa.bam \
	O=bam/$name.bwa.dedup.bam \
	M=mark_duplicates_log/$name.duplicateinfo.txt \
	samtools index bam/$name.bwa.dedup.bam
done < sample_names_list.txt

echo ========================================
echo creating sequence dictionary with picard
echo ========================================

$JAVA -jar $PICARD CreateSequenceDictionary R= $ref_genome/reference.fa O= $ref_genome/reference.dict
samtools faidx $ref_genome/reference.fa 

echo ===================================================================
echo calling SNPs with gatk HaplotypeCaller. this step will take a while
echo ===================================================================

mkdir -p gvcf gatk_log

while read name; do
	$GATK --java-options "-Xmx4g" HaplotypeCaller \
   -R $ref_genome \
   -I bam/$name.bwa.dedup.bam \
   --native-pair-hmm-threads 1 \
   -ERC GVCF \
   -O gvcf/$name.bwa.dedup.g.vcf
done < sample_names_list.txt

echo ================================================
echo importing into genomics DB with GenomicsDBImport
echo ================================================

for i in `ls gvcf | grep "vcf" | grep -v ".idx"| sed 's/.bwa.dedup.g.vcf//g'`; do
	echo -e "$i\tgvcf/$i.bwa.dedup.g.vcf"
done > $PROJECT.sample_map

$GATK --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $PWD/db \
       --batch-size 50 \
       -L Chr1 \
       --sample-name-map $PWD/$PROJECT.sample_map \
       --reader-threads 1

echo ===================================
echo calling variants with GenotypeGVCFs
echo ===================================

$gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $ref_genome/reference.fa \
   -V gendb://$PWD/db \
   -O $PWD/$PROJECT.vcf.gz

echo ==============================================================
echo filtering SNPs with SelectVariants.
echo keeping biallelic alleles, with minor allele frequency of 10%.
echo subsetting to 10%.
echo ==============================================================

$GATK SelectVariants \
-R $ref_genome/reference.fa \
-V $PWD/$PROJECT.vcf.gz \
-O $PWD/$PROJECT.snps.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
-fraction 0.1 \
-select "AF > 0.1" 

echo ========================================
echo converting vcf file to tab-delimited txt
echo ========================================

$GATK VariantsToTable \
-R $ref_genome/reference.fa \
-V $PWD/$PROJECT.snps.vcf.gz \
-F CHROM \
-F POS \
-GF GT \
-O $PWD/$PROJECT.snps.tab

cat $PWD/$PROJECT.snps.tab | sed 's/.GT	/	/g' | sed 's/.GT$//g' | sed 's|/||g' | sed 's/\.\./NN/g' | grep -v '*' > $PWD/$PROJECT.snps.formatted.tab
mv $PWD/$PROJECT.snps.formatted.tab $PWD/PROJECT.snps.tab

echo ================================================
echo looking at clustering in data with FastStructure
echo ================================================

mkdir -p faststructure

$PLINK \
--make-bed \
--vcf $PROJECT.snps.vcf.gz \
--out $PROJECT.snps \
--set-missing-var-ids @:# \
--double-id \
--allow-extra-chr

cd faststructure

for k in `seq 4`; do
	$PYTHON $STRUCTURE -K $k --input=../$PROJECT.snps --output=$PROJECT
done

echo getting best K value
$PYTHON $CHOOSEK --input=$PROJECT

echo ====================================
echo running PCA, plot will be in pca.pdf
echo ====================================

$R CMD BATCH pca.R

echo
echo finished pipeline!
echo
