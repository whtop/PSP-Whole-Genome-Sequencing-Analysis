#run SKAT by module

#choose a group
groupdir=PTV.01NA_I.1

R_file=/mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/SKAT/run_mod.Rscript;

bed_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/id/plink/${groupdir}/;
bed_head=r3.bm.chr;
chr='ALLauto'; #used for weight, nonbeta
bed_tail=.PspNhwCon.id.ex;
dat_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/;
fam_name=pspcon.2.fam;

#choose covariate file
covset=com3noh 
cov_name=pspcon.hwe.${covset}.cov;

#choose a set_dir
set_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/vcffilter/${groupdir}/kme/; #MPRA
set_head=PSP.ALL.super.; 
set_tail=.c.vcffilter.out;

#weight file
wt_dir=$set_dir; 
wt_head=PSP.;  
wt_tail=.BETAMAF.weight.out;

#choose an output directory
output_dir=/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/SKAT/${groupdir}/${covset}/super/;
mkdir -p $output_dir;
#choose an output_head
output_head=SKATO.super.betamaf;

super_list="Consensus_BannerProteomics";

for super in `echo $super_list | tr ',' '\n'`;do
  echo $super;
  #choose a full weight file
  wtfile=${wt_dir}${wt_head}${super}${wt_tail}; 

  runcallfile=${chr}.${super}.nowt_TMP.sh;
cat << EOF > $runcallfile
#!/usr/bin/env bash
thisfile=${output_dir}${output_head}chr${chr}.\$SGE_TASK_ID.res; if [ ! -f \"\${thisfile}\" ];
then
    export PATH="/mnt/analysis-psp/users/timchang/programs/miniconda/4.9.2/bin:\$PATH"  
    source activate /mnt/analysis-psp/users/timchang/programs/miniconda/R3.6.3
    Rscript $R_file $chr $bed_dir $bed_head $bed_tail $dat_dir $fam_name $cov_name $set_dir $set_head $set_tail $output_dir $output_head $super $wtfile > ${output_dir}${output_head}.${super}.CMD.out 2>&1
fi;
EOF

  echo $runcall;
  chmod 775 $runcallfile;
  qsub -cwd -N sk${super} -l h_data=30G,h_rt=48:59:00 -o ${output_dir}skat${output_head}.${super}.log -j y $runcallfile -wd /mnt/analysis-psp/users/timchang/v3/pc/2211/data/workingdir/;
  rm $runcallfile;
done;

