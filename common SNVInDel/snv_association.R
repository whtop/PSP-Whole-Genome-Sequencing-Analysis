library(SeqArray)
library(SeqVarTools)
library(SNPRelate)
library(Biobase)
library(GENESIS)

#### Inputs #####
vcffile <- paste0("gcad.qc.r3.wgs.16905.GATK.2021.08.24.biallelic.genotypes.chr", 1:22, ".ALL.vcf.gz")
convert vcf to gds
seqVCF2GDS(vcffile, 'biallelic.gds', fmt.import="GT", storage.option="LZMA_RA", parallel=12, verbose=TRUE)


gds <- seqOpen('biallelic.gds')
# read samples
phe <- read.csv("pheno.csv")
phe$sample.id <- phe$Sample_ID
phe <- phe[match(seqGetData(gds, "sample.id"), phe$Sample_ID),]
phe$Age <- as.integer(stringr::str_replace(phe$Age, '\\+', ""))
phe$Age[is.na(phe$Age)] <- median(phe$Age, na.rm=T)
phe$Age <- (phe$Age - mean(phe$Age))/sd(phe$Age)
phe$Sex <- plyr::mapvalues(phe$Sex, c(1,2), c('male', 'female'))
phe$Outcome <- as.integer(plyr::mapvalues(phe$Diagnosis, c('Control', 'PSP'), c(0, 1)))
metadata <- data.frame(labelDescription=names(phe), row.names=names(phe))
annot <- AnnotatedDataFrame(phe, metadata)
all.equal(annot$sample.id, seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)
samples <- read.csv("pspcon.2.fam", sep=" ", header=F)


# use KING for infer relationships
king <- snpgdsIBDKING(gds, sample.id=phe$Sample_ID, missing.rate=0.1, maf=0.001, num.thread = 16, verbose=FALSE)
saveRDS (object = king, file = 'data/king.RDS')


# PC-AiR, for PC calculation
king <- readRDS("data/king.RDS")
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
snpset <- snpgdsLDpruning(gds, method="corr", sample.id=phe$Sample_ID, missing.rate=0.1, maf=0.05, num.thread = 16, verbose=FALSE, slide.max.bp=10e6, ld.threshold=sqrt(0.1))
saveRDS (object = snpset, file = "data/ld_pruned.rds")
pruned <- unlist(snpset, use.names=FALSE)
pcs <- pcair(seqData, sample.include=phe$Sample_ID, kinobj=kingMat, kin.thresh=2^(-9/2), divobj=kingMat, div.thresh=-2^(-9/2), snp.include=pruned, num.cores=16)
# output PCs
pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
write.csv(pc.df, "data/pcs.csv")
saveRDS (object = pcs, file = "data/pcair.rds")

# pcrelate
snpset <- readRDS("./data/ld_pruned.rds")
pcs <- readRDS("./data/pcair.rds")
pruned <- unlist(snpset, use.names=FALSE)
seqSetFilter(seqData, sample.id=phe$Sample_ID, variant.id=pruned)
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)
pcrel <- pcrelate(iterator, pcs=pcs$vectors[,1:2], training.set=pcs$unrels, BPPARAM=BiocParallel::MulticoreParam(9))
saveRDS(object = pcrel, file = "data/pcrel.rds")
seqResetFilter(seqData, verbose=FALSE)

# association
pc.df <- read.csv("./data/pcs.csv")
pc.df <- dplyr::right_join(pc.df, pData(annot), by="sample.id")
pc.df <- pc.df[match(seqGetData(seqData, "sample.id"), pc.df$sample.id),]
annot <- AnnotatedDataFrame(pc.df)
all.equal(annot$sample.id, seqGetData(gds, "sample.id"))
sampleData(seqData) <- annot
pcrel <- readRDS("data/pcrel.rds")
# covariance matrix from pcrelate output
grm <- pcrelateToMatrix(pcrel, scaleKin=2)
# fit the null model
seqSetFilter(seqData, sample.id=samples$V2)
nullmod <- fitNullModel(seqData, outcome="Outcome", family='binomial',
                        covars=c("Sex", paste0("PC", 1:5)),
                        cov.mat=grm, verbose=FALSE)
seqSetFilter(seqData, sample.id=samples$V2)
seqSetFilterCond(seqData, maf = 0.01, missing.rate = 0.1, parallel=5)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE, BPPARAM=BiocParallel::MulticoreParam(6))

# add SNV informatio
seqResetFilter(gds)
seqSetFilter(seqData, sample.id=samples$V2)
seqSetFilterCond(gds, maf = 0.01, missing.rate = 0.1, parallel=6)
alleles <- data.frame(row.names = seqGetData(gds, "variant.id"),seqGetData(gds, 'allele'), seqGetData(gds, 'annotation/id'), seqGetData(gds, 'annotation/info/VariantType'))
colnames(alleles) <- c("V1", 'V2', 'V3')
all(assoc$variant.id %in% rownames(alleles))
# must change to character
assoc$alleles <- alleles[as.character(assoc$variant.id), "V1"]
assoc$ref <- sapply(strsplit(assoc$alleles, ","), `[`, 1)
assoc$alt <- sapply(strsplit(assoc$alleles, ","), `[`, 2)
assoc$rsid <- alleles[as.character(assoc$variant.id), "V2"]
assoc$variant.type <- alleles[as.character(assoc$variant.id), "V3"]
write.csv(assoc, file = "data/associations.csv")
