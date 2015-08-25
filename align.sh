#copy and rename MY18 files
 cp data/fastq/150625_I188_FCC5V31ANXX_L6_wHAIPI019442-89_2.fq.gz data/fastq/MY18_2.fq.gz

# qc
STRAIN=MY18
fastqc --threads 4 `ls data/fastq/$STRAIN*`
# bwa
bwa mem -t 4 reference/WS245/WS245.fa.gz data/fastq/$STRAIN_1.fq.gz data/fastq/$STRAIN_2.fq.gz | \
samtools sort -O bam -T $STRAIN -@ 4 - > data/bam/$STRAIN.bam
samtools index data/bam/$STRAIN.bam

# QC bam file
fastqc data/bam/$STRAIN.bam

# Variant calling v2
REF=reference/WS245/WS245.fa.gz
BAM=data/bam/MY10.bam
STRAIN=MY10
chrom_regions=`samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | sed "s/$/.${STRAIN}.vcf.gz/g" | sed 's/^/data\/vcf\//g'`
samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 6 sh -c "samtools mpileup -BQ0 -uf $REF -r \"{}\" $BAM | bcftools call -O z -mv > data/vcf/\"{}\".$STRAIN.vcf.gz"
bcftools concat -O z ${chrom_regions} > data/vcf/$STRAIN.vcf.gz
bcftools index data/vcf/$STRAIN.vcf.gz
rm ${chrom_regions}
