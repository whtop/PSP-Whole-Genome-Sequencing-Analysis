#!/bin/bash
#$ -S /bin/bash
#$ -N normann
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pa/2211/data/vcfs/bi/norm/id.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=2G,h_rt=72:00:00
#$ -tc 10
#$ -hold_jid bcf_pass
#$ -t 1:22

#psp con samples only
#keep min-ac =1
#filter for qual = pass only
#normalize indel only. no multiallelic to split as this is biallelic
#annotate ID's


chr=($(seq 1 1 22))

echo ${chr[${SGE_TASK_ID}-1]}

indir=/mnt/analysis-psp/users/wanping.lee/qc_vcf/bi/
infilehead=gcad.qc.r3.wgs.16905.GATK.2021.08.24.biallelic.genotypes.chr
infiletail=.ALL.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
samplefile=/mnt/analysis-psp/users/timchang/v3/pa/2211/data/clin/mod/pspcon.1.samp
reffile=/mnt/analysis-psp/users/timchang/v3/pa/2211/data/ref/fasta/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
outdir=/mnt/analysis-psp/users/timchang/v3/pa/2211/data/vcfs/bi/norm/
outfilehead=r3.bi.chr
outfiletail=.PspNhwCon.QC.f1.vcf.gz
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}

bcftools view -O u -S $samplefile $infile | bcftools norm -m -both --check-ref w -f $reffile -Ou  |  bcftools annotate --set-id '%CHROM\_%POS\_%REF\/%FIRST_ALT' -O z -o $outfile; 
tabix  -p vcf $outfile;
bcftools index -n $infile;
bcftools index -n $outfile;

