#!/bin/bash
#$ -S /bin/bash
#$ -N vep
#$ -cwd
#$ -o /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/none/vep/vep.$TASK_ID.log
#$ -j y
#$ -m bea
#$ -l h_data=8G,h_rt=48:00:00
#$ -hold_jid ann
#$ -t 1:22


chr=($(seq 1 1 22) X Y M)

echo ${chr[${SGE_TASK_ID}-1]}

indir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/none/annovar/
infilehead=r3.bm.chr
infiletail=.PspNhwCon.none.ann.hg38_multianno.vcf.gz
infile=${indir}${infilehead}${chr[${SGE_TASK_ID}-1]}${infiletail} 
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/none/vep/
outfilehead=r3.bm.chr
outfiletail=.PspNhwCon.none.vep.vcf.gz
outfile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${outfiletail}
donefiletail=.PspNhwCon.none.vep.done
donefile=${outdir}${outfilehead}${chr[${SGE_TASK_ID}-1]}${donefiletail}

#get conda in path. to activate use source
export PATH="/mnt/analysis-psp/users/timchang/programs/miniconda/4.9.2/bin:$PATH"
echo $PATH
if [[ -v CONDA_PREFIX ]]; then
    if [ $CONDA_PREFIX != "/mnt/analysis-psp/users/timchang/programs/miniconda/vep" ]; then
        source deactive;
        source activate /mnt/analysis-psp/users/timchang/programs/miniconda/vep 
    fi;
else
    source activate /mnt/analysis-psp/users/timchang/programs/miniconda/vep
fi;

VEP_annotation_dbs=/mnt/analysis-psp/users/timchang/programs/vep/104/cache/
VEP_plugins=/mnt/analysis-psp/users/timchang/programs/vep/104/plugins/
export PERL5LIB=/mnt/analysis-psp/users/timchang/programs/vep/104/plugins/loftee



vep --force_overwrite --format vcf --assembly GRCh38 --verbose --polyphen p --humdiv --ccds --hgvs --symbol --numbers --canonical --protein --biotype --uniprot --var_synonyms --variant_class --fork 1 --buffer_size 10000 --offline -i $infile -o $outfile --cache --dir $VEP_annotation_dbs --plugin LoF,loftee_path:/mnt/analysis-psp/users/timchang/programs/vep/104/plugins/loftee --plugin CADD,/mnt/analysis-psp/users/timchang/programs/vep/104/cadd/1.6/whole_genome_SNVs.tsv.gz --dir_plugins /mnt/analysis-psp/users/timchang/programs/vep/104/plugins/ --pick --vcf --compress_output bgzip --gencode_basic --no_stats && touch $done_file; 
tabix  -p vcf ${outfile};


