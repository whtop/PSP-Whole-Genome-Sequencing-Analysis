#!/bin/bash
#$ -S /bin/bash
#$ -N remVar
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/remVar.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=2G,h_rt=48:00:00
#$ -hold_jid hwe
#$ -t 1:22

#remove Variants from list

chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}

indir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/
infilehead=r3.bm.chr
infiletail=.PspNhwCon.eur.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
varexdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/
varexhead=r3.bm.chr
varextail=.PspNhwCon.eur.id.hwe.rem
varexfile=${varexdir}${varexhead}${chr[${SGE_TASK_ID}-1]}${varextail}
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/
outfilehead=r3.bm.chr
outfiletail=.PspNhwCon.hwe.vcf.gz
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}

bcftools view --exclude ID=@$varexfile -O z -o $outfile $infile;
tabix  -p vcf $outfile;
bcftools index -n $infile;
bcftools index -n $outfile;

