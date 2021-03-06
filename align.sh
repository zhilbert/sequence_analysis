#!/usr/bin/env


#copy and rename MY18 files
#cp data/fastq/150625_I188_FCC5V31ANXX_L6_wHAIPI019442-89_2.fq.gz data/fastq/MY18_2.fq.gz

# qc
STRAIN=${1}
# fastqc --threads 4 `ls data/fastq/$STRAIN*fq.gz`
# bwa
# bwa mem -t 4 reference/WS245/WS245.fa.gz data/fastq/${STRAIN}_1.fq.gz data/fastq/${STRAIN}_2.fq.gz | \
# samtools sort -O bam -T $STRAIN -@ 4 - > data/bam/$STRAIN.bam
# samtools index data/bam/$STRAIN.bam

#Remove duplicates
picard MarkDuplicates \
        I=data/bam/${STRAIN}.bam \
        O=data/bam/${STRAIN}.new.bam \
        M=data/bam/${STRAIN}.duplicate_report.txt \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false

#Find structural variants
delly -g reference/WS245/WS245.fa -o data/sv/${STRAIN}.sv.vcf data/bam/${STRAIN}.bam
bcftools view -O z data/sv/${STRAIN}.sv.vcf > data/sv/${STRAIN}.sv.vcf.gz
bcftools index data/sv/${STRAIN}.sv.vcf.gz
bcftools filter -O z --include '%FILTER=="PASS"' data/sv/${STRAIN}.sv.vcf.gz > data/sv/${STRAIN}.filtered.sv.vcf.gz

# Variant calling v2
REF=reference/WS245/WS245.fa.gz
BAM=data/bam/${STRAIN}.bam
chrom_regions=`samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | sed "s/$/.${STRAIN}.vcf.gz/g" | sed 's/^/data\/vcf\//g'`
samtools view -H $BAM | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 6 sh -c "samtools mpileup -BQ0 -uf $REF -r \"{}\" $BAM | bcftools call -O z -mv > data/vcf/\"{}\".$STRAIN.vcf.gz"
bcftools concat -O z ${chrom_regions} | \
bcftools filter --mode "+" --soft-filter 'Low coverage' --include 'DP>8' | \
bcftools filter --mode "+" --soft-filter 'High coverage' --include 'DP<2000'| \
bcftools filter -O z --mode "+" --soft-filter 'HET' --include 'FORMAT/GT=="1/1"' > data/vcf/$STRAIN.vcf.gz
bcftools index data/vcf/$STRAIN.vcf.gz
rm ${chrom_regions}

#stats
bcftools stats data/vcf/${STRAIN}.vcf.gz > data/vcf/${STRAIN}.stats.txt
plot-vcfstats -p data/vcf/${STRAIN}/ data/vcf/${STRAIN}.stats.txt

#snpeff
#snpeff databases > databases.txt
#snpeff download WS241
snpeff eff -s data/snpeff/${STRAIN}_snpEff_summary.html -no-downstream -no-intergenic -no-intron -no-upstream WS241 data/vcf/${STRAIN}.vcf.gz > data/snpeff/${STRAIN}.eff.vcf

#run snpeff on region of interest
region_of_interest=IV:5464516-5561950
#region_of_interest=X:2144602-8455640
VCF=data/vcf/${STRAIN}.vcf.gz
bcftools view -r ${region_of_interest} ${VCF} | snpeff eff -s data/snpeff/${STRAIN}_roi_snpEff_summary.html WS241 > data/snpeff/${STRAIN}_roi.eff.vcf

#convert vcf to tsv
#python snpeff2tsv.py {region} {vcf_file} {severity} > ${STRAIN}_severity.tsv
#cat data/snpeff/MY10.tsv <(cat data/snpeff/MY10.tsv | grep -v 'feature_id') | sort -k1,1 -k2,2n

#merge MY10 with WI vcf and find unique MY10 SNPs in region
#bcftools merge -O z data/vcf/MY10.vcf.gz data/vcf/MY18.vcf.gz data/WI/20150731_WI_PASS.vcf.gz > merged.vcf.gz
#bcftools view --types snps data/WI/merged.vcf.gz IV:5464516-5561950  | \
 #awk '$10 ~ "1\/1.*" { print } $0 ~ "#" { print }' | \
   #bcftools filter --include 'INFO/AC == 4 && AN > 4' |\
    #bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]'

bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\n]' data/WI/merged.vcf.gz | grep '1/1' > data/WI/merged_variants.txt
