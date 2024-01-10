#bigsnpr to calculate pca
#conda activate /mnt/analysis-psp/users/timchang/programs/miniconda/R3.6.3
#conda install r-codetools r-matrix
.libPaths('/mnt/analysis-psp/users/timchang/programs/R/3.6.3/')
#also install bigsnpra, RhpcBLASctl
library(bigsnpr);library(bigreadr);library(data.table)

geno_path='/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/m1000G_20181203/maf.01/plink/merge/'
geno_file=paste0(geno_path,'r3.bm.chrALLauto.PspNhwCon.1KG.maf.01.s.namb.g')
#1000G fam
onekgfam_file='/mnt/analysis-psp/users/timchang/v3/pc/2211/data/ref/1000G/20181203/norm.filt.id/plink/ALL.chrALLauto.maf.01.s.namb.1.fam'
#pca out
pca_outfile=paste0(geno_path,'r3.bm.chrALLauto.PspNhwCon.1KG.maf.01.s.namb.g.projpca','.rds')

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
NCORES <- 4 

#read onekg fam file
onekgfam=read.table(onekgfam_file);
names(onekgfam)=c('FID','IID','mID','pID','sex','diagnosis')

#read bed do once
snp_readBed(paste0(geno_file,'.bed'))
# now attach the genotype object
obj.bigSNP <- snp_attach(paste0(geno_file,".rds"))
#read in bed file
obj.bed=bed(paste0(geno_file,'.bed'))
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
#renaming
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
fam = obj.bigSNP$fam

#default threshold r2 0.2, size is 100/thr.r2 = 200. also removes outlier snps
#SVD for only onekg
indonekg = which(fam$sample.ID %in% onekgfam$IID)
obj.bed.autoSVD = bed_autoSVD(obj.bed=obj.bed, ind.row=indonekg,k=20,ncores=NCORES)

#project everyone onto onekg space
scores = bed_projectSelfPCA(obj.bed.autoSVD,obj.bed=obj.bed,ind.row=1:dim(obj.bed$fam)[1])

save('scores','obj.bed.autoSVD','fam',file=pca_outfile)

