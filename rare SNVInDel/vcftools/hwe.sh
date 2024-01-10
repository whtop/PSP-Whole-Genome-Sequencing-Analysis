#!/bin/bash
#$ -S /bin/bash
#$ -N hwe
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=7G,h_rt=48:00:00
#$ -hold_jid keep
#$ -t 1:22

#find hwe less than 1E-7

chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}

indir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/
infilehead=r3.bm.chr
infiletail=.PspNhwCon.eur.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
sample_names=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/con.2.sampnames
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/
outfilehead=r3.bm.chr
outfiletail=.PspNhwCon.eur
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}


vcftools --gzvcf $infile --keep $sample_names --hardy --out $outfile; 
bcftools query -f '%ID\n' $infile > ${outfile}.id
sed -i '1s/^/ID\n/' ${outfile}.id
paste ${outfile}.id ${outfile}.hwe > ${outfile}.id.hwe
awk -v OFS='\t' '$7<1E-7 {print $1}' ${outfile}.id.hwe > ${outfile}.id.hwe.rem; 


