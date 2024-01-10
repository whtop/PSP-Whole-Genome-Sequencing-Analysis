#!/bin/bash
#$ -S /bin/bash
#$ -N bcf_norm
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/multi/norm/multi.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=6G,h_rt=48:00:00
#$ -tc 16
#$ -hold_jid bcf_multi
#$ -t 1;22


#norm left align
#remove spanning deletion
#annotate ID's

chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}
indir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/multi/samp/
infilehead=r3.multi.chr
infiletail=.PspNhwCon.samp.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
reffile=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/ref/fasta/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/multi/norm/
outfilehead=r3.multi.chr
outfiletail=.PspNhwCon.f1.vcf.gz
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}

bcftools norm --check-ref w -f $reffile -Ou $infile | bcftools view -e 'ALT="*"' -Ou | bcftools annotate --set-id '%CHROM\_%POS\_%REF\/%FIRST_ALT' -Oz -o $outfile

tabix -p vcf $outfile;
bcftools index -n $infile;
bcftools index -n $outfile;

