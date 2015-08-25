#!usr/local/bin/python
from subprocess import Popen, PIPE
import re
from collections import namedtuple
from collections import defaultdict
import os
script_dir = os.path.dirname(os.path.realpath(__file__))



#=================================#
# Load Wormbase gene identifiers. #
#=================================#
import sys
import cPickle
"""
import urllib2
import StringIO
import gzip
URL = "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.current.geneIDs.txt.gz"
response = urllib2.urlopen(URL)
compressedFile = StringIO.StringIO()
compressedFile.write(response.read())
compressedFile.seek(0)
decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb').read()
wb_gene = dict([x.split(",")[1:3] for x in decompressedFile.splitlines()])

cPickle.dump(wb_gene, file = open("wb_genes.data", "w"))

# Write gene names file for R.
with open("wb_gene.txt", "w") as f:
    wb_gene_items = wb_gene.items()
    wb_gene_items = ["\t".join(x) for x in wb_gene_items]
    wb_gene_items = ["ID\tname"] + wb_gene_items
    f.write("\n".join(wb_gene_items))

"""

if len(sys.argv) == 1:
  region = "I:1-100000"
  vcf_file = "20150608_WI_PASS.snpeff.vcf.gz"
  severity = ["HIGH","MODERATE","LOW"]
else:
  region = sys.argv[1]
  vcf_file = sys.argv[2]
  severity = sys.argv[3].split(",")
if region == "ALL":
    comm = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%ANN[\t%TGT]\n", vcf_file]
else:
    comm = ["bcftools", "query", "--regions", region, "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%ANN[\t%TGT]\n", vcf_file]

proc = Popen(comm, stdout = PIPE)
last_POS = 0

base_header = ["CHROM",
               "POS",
               "REF",
               "ALT",
               "FILTER"]

ANN_header = ["allele",
              "effect",
              "impact",
              "gene_name",
              "gene_id",
              "feature_type",
              "feature_id",
              "transcript_biotype",
              "exon_intron_rank",
              "nt_change",
              "aa_change",
              "cDNA_position/cDNA_len",
              "protein_position",
              "distance_to_feature",
              "error"]

sample_q = ["bcftools", "query", "--list-samples", vcf_file]
sample_names = Popen(sample_q, stdout = PIPE).communicate()[0].strip().split("\n")

header = base_header + ANN_header + sample_names
print '\t'.join(header)

wb_gene = cPickle.load(open(script_dir + "/wb_genes.data","r"))


for line in proc.stdout:
  line = line.strip().split("\t")
  CHROM  = line[0]
  POS    = line[1]
  REF    = line[2]
  ALT    = line[3]
  FILTER = line[4]
  ANN    = line[5]
  GT     = line[6:]
  if last_POS != POS:
    for ANN_COL in [x.split("|") for x in ANN.strip().split(",")]:
      if ANN_COL == ['.']:
        ANN_COL = [""] * len(ANN_header)
      if ANN_COL[2] in severity:
        ANN_COL = ANN_COL[:-1]
        ANN_COL[3] = wb_gene[ANN_COL[3]]
        newline = [CHROM, POS, REF, ALT, FILTER] + ANN_COL + GT
        print '\t'.join(newline)
