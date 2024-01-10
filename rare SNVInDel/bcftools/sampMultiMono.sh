#!/bin/bash
#$ -S /bin/bash
#$ -N bcf_mult
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/multi/samp/multi.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=6G,h_rt=48:00:00
#$ -tc 16
#$ -t 1;22


#psp con samples
#remove monomorphic

chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}
indir=/mnt/analysis-psp/users/wanping.lee/qc_vcf/multi/subset/
infilehead=gcad.preview.r3.wgs.16906.GATK.2020.05.26.multiallelic.genotypes.chr
infiletail=.ALL.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
samplefile=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/pspcon.1.samp
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/multi/samp/
outfilehead=r3.multi.chr
outfiletail=.PspNhwCon.samp.vcf.gz
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}

bcftools view -S $samplefile -Ou $infile | bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,NS | bcftools filter -e 'AC==0 || AC==AN' -O z -o $outfile 
tabix -p vcf $outfile;
bcftools index -n $infile;
bcftools index -n $outfile;

