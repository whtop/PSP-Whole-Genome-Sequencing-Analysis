#Main pipeline steps for rare SNV/InDel analysis
#see subdirectories for code used

##Processing##
#==biallelic
#norm indel, annotate. (note set missing DP<10 or GQ<20 by Wan-Ping)
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/bcftools/sampNormAnn.sh

#==multiallelic
#psp/con, remove mono
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/bcftools/sampMultiMono.sh

#from multi. norm split multi then left align indel, rm *, annotate
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/bcftools/normAnnMulti.sh

#==combine biallelic and multiallelic
#cat splitmulti with bi. concat, sort, remove exact duplicates (keeps first), fill tags, rem monomorphic 
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/bcftools/catSMultiBi.sh

#run PCA
/mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/bigsnpr_pca.R 
#analyze
/mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/bigsnpr_pca_anal.R

#keep eur sample, missing>10%, fill tags, rem monomorphic 
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/bcftools/keepSampleVarMissMon.sh
#remove hwe
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/vcftools/hwe.sh
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/bcftools/removeVar_hwe.sh

#pca eur com 
Rscript /mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/bigsnpr_pca.r2.R 

#annotation: annovar
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/annotate/annovar.sh
#annotation: vep
qsub /mnt/analysis-psp/users/timchang/v3/pc/2211/code/annotate/vep.sh

##GCTA##
#GCTA LDMS
/mnt/analysis-psp/users/timchang/v3/pc/2211/code/gcta/run_gcta.sh

##SKATO##
#create weights based on MAF
Rscript /mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/ptn_filterWeight.R

#run gene by chr
/mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/SKAT/run_gene_bychr.sh

#run mod modules
/mnt/analysis-psp/users/timchang/v3/pc/2211/code/R/SKAT/run_mod.sh



