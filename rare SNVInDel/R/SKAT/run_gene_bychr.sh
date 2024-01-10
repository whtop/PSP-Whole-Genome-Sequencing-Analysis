#run SKATO for genes

#variant group directory name
groupdir=PTV.01NA_I.1;

R_file=/mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/SKAT/run_gene.Rscript;

bed_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/id/plink/;
bed_head=r3.bm.;
bed_tail=.PspNhwCon.id;
dat_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/;
fam_name=pspcon.2.fam;

covset=com3 
cov_name=pspcon.hwe.${covset}.cov;

set_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/vcffilter/${groupdir}/; #ptn and non

set_head=PSP.; 
set_tail=.c.vcffilter.out;
wt_dir=${set_dir};
wt_head=PSP.chr; 
wt_tail=ALL.BETA.weight.out; 

output_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/SKAT/${groupdir}/${covset}/weightbeta/; 

mkdir -p $output_dir;
output_head=SKATO.BETA.;
chr_list="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22";
numsplit=3;

for chr in `echo $chr_list | tr ',' '\n'`;do    
    echo $chr;
    runcallfile=${output_dir}/${chr}.call.sh
cat << EOF > $runcallfile
#!/usr/bin/env bash
thisfile=${output_dir}${output_head}chr${chr}.\$SGE_TASK_ID.res; if [ ! -f \"\${thisfile}\" ];
then
    export PATH="/mnt/analysis-psp/users/timchang/programs/miniconda/4.9.2/bin:\$PATH"
    source activate /mnt/analysis-psp/users/timchang/programs/miniconda/R3.6.3
	Rscript $R_file chr${chr} $bed_dir $bed_head $bed_tail $dat_dir $fam_name $cov_name $set_dir $set_head $set_tail $output_dir $output_head $wt_dir $wt_head $wt_tail \$SGE_TASK_ID > ${output_dir}${chr}.\$SGE_TASK_ID.wt.out 2>&1; fi;
fi;
EOF
    #echo $runcall
    chmod 775 $runcallfile;     
    cat $runcallfile
	qsub -cwd -N skW${chr} -t 1:$numsplit -l h_data=7000M,h_rt=48:59:00 -o ${output_dir}skat${chr}.\$TASK_ID.log -j y $runcallfile -wd /mnt/analysis-psp/users/timchang/v3/pc/2211/data/workingdir ; #enhHiC
 
done;                                


