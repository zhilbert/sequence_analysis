#!/bin/bash


# prepare reference
mkdir reference
mkdir reference/WS245
curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS245.genomic.fa.gz > reference/WS245/WS245.fa.gz
gunzip -f reference/WS245/WS245.fa.gz
bgzip --stdout reference/WS245/WS245.fa > reference/WS245/WS245.fa.gz

# Make bwa index
bwa index reference/WS245/WS245.fa.gz

# faidx
samtools faidx reference/WS245/WS245.fa.gz

# qc
fastqc --threads 4 `ls data/fastq/*`

# bwa
bwa mem -t 4 reference/WS245/WS245.fa.gz data/fastq/140314_I315_FCC3NFRACXX_L2_WHAIPI002516-13_1.fq.gz data/fastq/140314_I315_FCC3NFRACXX_L2_WHAIPI002516-13_2.fq.gz | \
samtools sort -O bam -T MY10 -@ 4 - > data/bam/MY10.bam
samtools index data/bam/MY10.bam

# QC bam file
fastqc data/bam/MY10.bam

# Variant calling v2
REF=reference/WS245/WS245.fa.gz
BAM=data/bam/MY10.bam
STRAIN=MY10
chrom_regions=`samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | sed "s/$/.${STRAIN}.vcf.gz/g" | sed 's/^/data\/vcf\//g'`
samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 6 sh -c "samtools mpileup -BQ0 -uf $REF -r \"{}\" $BAM | bcftools call -O z -mv > data/vcf/\"{}\".$STRAIN.vcf.gz"
bcftools concat -O z ${chrom_regions} > data/vcf/$STRAIN.vcf.gz
bcftools index data/vcf/$STRAIN.vcf.gz
rm ${chrom_regions}
