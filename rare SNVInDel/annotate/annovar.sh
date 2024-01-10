#!/bin/bash
#$ -S /bin/bash
#$ -N ann
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/none/annovar/ann.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=12G,h_rt=48:00:00
#$ -hold_jid remSamp
#$ -t 1:22

#remove Variants from list

chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}

indir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/none/
infilehead=r3.bm.chr
infiletail=.PspNhwCon.none.vcf
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/none/annovar/
outfilehead=r3.bm.chr
outfiletail=.PspNhwCon.none.ann
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}

annovar_call=/mnt/analysis-psp/users/timchang/programs/annovar/annovar/table_annovar.pl
annovar_annotation_tables=/mnt/analysis-psp/users/timchang/programs/annovar/annovar/humandb/

$annovar_call --vcfinput $infile $annovar_annotation_tables --buildver hg38 --nocheckfile --remove --nastring . --outfile $outfile --protocol cytoBand,gnomad30_genome,gnomad211_exome,dbnsfp42c --operation r,f,f,f 
bgzip -c ${outfile}.hg38_multianno.vcf > ${outfile}.hg38_multianno.vcf.gz
tabix  -p vcf ${outfile}.hg38_multianno.vcf.gz;


