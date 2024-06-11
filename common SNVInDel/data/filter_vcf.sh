#!/usr/bin/bash

sample=/home/data/sample.list
OUT_DIR=/home/data/snp/vcfs
# get the basename of vcf
OUT_PREFIX=$(basename $VCF)
# remove ".vcf.gz" at the end of prefix
OUT_PREFIX=${OUT_PREFIX%%.vcf.gz}

THREADS=2
# remove everything except GT
bcftools annotate --threads $THREADS -x "FORMAT" -Oz -o $OUT_DIR/tmp/$OUT_PREFIX.annot.vcf.gz $VCF
# only keep PASS SNPs
bcftools view --threads $THREADS -S $sample -f PASS -Oz -o $OUT_DIR/tmp/$OUT_PREFIX.subset.vcf.gz $OUT_DIR/tmp/$OUT_PREFIX.annot.vcf.gz
# use VFLAGS and ABHet to further filter variants
bcftools view --threads $THREADS -i "VFLAGS_one_subgroup=0 && ((ABHet_one_subgroup > 0.25 && ABHet_one_subgroup < 0.75) || ABHet_one_subgroup = '.') && AC > 1" -Oz -o $OUT_DIR/$OUT_PREFIX.vcf.gz $OUT_DIR/tmp/$OUT_PREFIX.subset.vcf.gz
bcftools index --threads $THREADS $OUT_DIR/$OUT_PREFIX.vcf.gz
