library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(grid)
library(gridExtra)

#Read in data and fix header names
command <- "bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%SAMPLE=%GT]\n' data/vcf/MY10.vcf.gz"
df_my10 <- read_tsv(pipe(command))
names(df_my10) <-  gsub(".bam:GT","",gsub("SAMPLE=data/bam/","", gsub("(# )?\\[[0-9]+\\]", "", names(df_my10))))

command <- "bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%SAMPLE=%GT]\n' data/vcf/MY18.vcf.gz"
df_my18 <- read_tsv(pipe(command))
names(df_my18) <-  gsub(".bam:GT","",gsub("SAMPLE=data/bam/","", gsub("(# )?\\[[0-9]+\\]", "", names(df_my18))))

#Plot all variant calls, show those that are het or have low/high coverage
p1 <- ggplot(df_my10) + 
  geom_histogram(aes(x = POS, fill=FILTER)) +
  facet_grid(.~CHROM)+
  scale_fill_brewer(palette= "Spectral")+
  labs(x= "Position", y="Number of Variants", title="MY10")+
  scale_x_continuous(labels = NULL)+
  theme_bw()

p2 <-ggplot(df_my18) + 
  geom_histogram(aes(x = POS, fill=FILTER)) +
  facet_grid(.~CHROM)+
  scale_fill_brewer(palette= "Spectral")+
  labs(x= "Position", y="Number of Variants", title="MY18")+
  scale_x_continuous(labels = NULL)+
  theme_bw()

#Show MY18 and MY10 on the same plot
grid.arrange(p1,p2, ncol=1)


#Read in variants with predicted effects (only moderate and high effect variants)
df_my18_snpeff <- read_tsv("data/snpeff/MY18.tsv")
df_my10_snpeff <- read_tsv("data/snpeff/MY10.tsv")

#split effect column into the primary effect and secondary effect to cut down on effect categories
df_my18_split <- df_my18_snpeff %>%
  select(CHROM, POS, effect) %>%
  separate(effect,
           into = c("main_effect", "secondary_effect"),
           sep = "&",
           extra = "drop")

df_my10_split <- df_my10_snpeff %>%
  select(CHROM, POS, effect) %>%
  separate(effect,
           into = c("main_effect", "secondary_effect"),
           sep = "&",
           extra = "drop")

#plot variants by the main effect type  
p3 <-ggplot(df_my18_split) +
  geom_histogram(aes(x = POS, fill = main_effect)) +
  facet_grid(.~CHROM)+
  labs(x= "Position", y="Number of Variants", title="MY18")+
  scale_x_continuous(labels = NULL)+
  scale_fill_discrete(name = "Effect of SNP", labels = c("Disruptive In Frame Deletion","Disruptive In Frame Insertion", "Frameshift", "In Frame Deletion","In Frame Insertion", "Missense", "Start Lost", "Stop Gained"))
  
p4 <-ggplot(df_my10_split) +
  geom_histogram(aes(x = POS, fill = main_effect)) +
  facet_grid(.~CHROM)+
  labs(x= "Position", y="Number of Variants", title="MY10")+
  scale_x_continuous(labels = NULL)+
  theme(legend.position="none")

#Show MY18 and MY10 on the same plot
grid.arrange(p3,p4, ncol=1)

genotypes <- tbl_df(fread("~/Documents/git/sequence_analysis/data/WI/merged_variants.txt")) %>%
  rename(CHROM=V1,
         POS=V2,
         STRAIN=V3,
         GT=V4)

MY10_POS <- filter(genotypes, STRAIN %in% c("MY10","data/bam/MY10.bam"), CHROM=="IV", POS>5464516, POS<5561950)$POS

filtered_genotypes <-   filter(genotypes, POS %in% MY10_POS) %>%
  group_by(CHROM, POS) %>%
  mutate(n = n())

filtered_genotypes <-   filter(genotypes,CHROM == "X", POS %in% MY18_POS) %>%
  group_by(CHROM, POS) %>%
  mutate(n = n())
