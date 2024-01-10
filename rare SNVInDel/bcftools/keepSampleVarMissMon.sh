#!/bin/bash
#$ -S /bin/bash
#$ -N keep
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/keep.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=8G,h_rt=48:00:00
#$ -t 1:22

#keep sample 
#keep <10% missing variants
#rmove monomorphic
#reannotate INFO

chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}

indir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/
infilehead=r3.bm.chr
infiletail=.PspNhwCon.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
samplefile=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/m1000G_20181203/maf.01/plink/merge/r3.bm.chrALLauto.PspNhwCon.1KG.maf.01.s.namb.g.projpca.6.8.keep.txt
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/
outfilehead=r3.bm.chr
outfiletail=.PspNhwCon.eur.vcf.gz
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}

bcftools view -O u -S $samplefile $infile | bcftools view --min-ac=1 --include 'F_MISSING<0.1' -Ou |  bcftools +fill-tags -Ou -- -t AC,AF,AN,MAF,NS | bcftools filter -e 'AC==0 || AC==AN' -Oz  -o $outfile;

tabix  -p vcf $outfile;
bcftools index -n $infile;
bcftools index -n $outfile;

