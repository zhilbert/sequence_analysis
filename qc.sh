#!/bin/bash


# prepare reference
mkdir reference
mkdir reference/WS245
cd reference/WS245
curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS245.genomic.fa.gz > WS245.fa.gz
gunzip -f WS245.fa.gz
bgzip --stdout WS245.fa > WS245.fa.gz

# Make bwa index
bwa index WS245.fa.gz

# qc
fastqc --threads 4 `ls data/fastq/*`

# bwa
bwa mem -t 4 reference/WS245/WS245.fa.gz data/fastq/140314_I315_FCC3NFRACXX_L2_WHAIPI002516-13_1.fq.gz data/fastq/140314_I315_FCC3NFRACXX_L2_WHAIPI002516-13_2.fq.gz | \
samtools sort -O bam -T MY10 -@ 4 - > data/bam/MY10.bam
samtools index data/bam/MY10.bam
