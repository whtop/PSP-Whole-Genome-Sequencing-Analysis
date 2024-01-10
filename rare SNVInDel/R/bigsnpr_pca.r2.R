#bigsnpr to calculate pca
#conda activate /mnt/analysis-psp/users/timchang/programs/miniconda/R3.6.3
.libPaths('/mnt/analysis-psp/users/timchang/programs/R/3.6.3/')
library(bigsnpr);library(bigreadr);library(data.table)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)


geno_path='/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/id/plink/'
genocom_tail='r3.bm.chrALLauto.PspNhwCon.id.com'
genocom_file=paste0(geno_path,genocom_tail)
genorare_tail='r3.bm.chrALLauto.PspNhwCon.id.rare'
genorare_file=paste0(geno_path,genorare_tail)

#pca out
pca_outfile=paste0(geno_path,'r3.bm.chrALLauto.PspNhwCon.id.comrare','.projpca','.rds')

NCORES <- 12
numPC=20


#read bed do once
snp_readBed(paste0(genocom_file,'.bed'))
# now attach the genotype object
obj.bigSNP.com <- snp_attach(paste0(genocom_file,".rds"))
#read in bed file
obj.bed.com=bed(paste0(genocom_file,'.bed'))
# extract the SNP information from the genotype
map.com <- obj.bigSNP.com$map[-3]
names(map.com) <- c("chr", "rsid", "pos", "a1", "a0")
#renaming
genotype.com <- obj.bigSNP.com$genotypes
# Rename the data structures
CHR.com <- map.com$chr
POS.com <- map.com$pos
fam.com = obj.bigSNP.com$fam


#SVD for only onekg
obj.bed.autoSVD.com = bed_autoSVD(obj.bed=obj.bed.com, k=20,ncores=NCORES)

#project everyone onto onekg space
scores.com = predict(obj.bed.autoSVD.com)
#add ID
scores.com.id=as.data.frame(cbind(scores.com,fam.com[,c('sample.ID','family.ID','affection')]))
names(scores.com.id)=c(paste0('PC',1:numPC),'sample.ID','family.ID','diagnosis')

save('scores.com.id','scores.com','obj.bed.autoSVD.com',file=pca_outfile)

