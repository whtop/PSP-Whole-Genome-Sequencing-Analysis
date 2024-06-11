library(SeqArray)
library(SeqVarTools)
library(SNPRelate)
library(Biobase)
library(GENESIS)

##### Inputs #####

seqBED2GDS(bed.fn = "psp.bed", bim.fn = "psp.bim", fam.fn = "psp.fam", parallel=12, out.gdsfn = "psp.gds")

gds <- seqOpen('psp.gds')
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

# association
pc.df <- read.csv("./sv_pcs/psp.eigenvec", sep="", header=F)
pc.snp <- read.csv("./snv_pcs/pcs.csv")
colnames(pc.df) <- c(c('IID', 'sample.id'), paste0('SV_PC', 1:20))
pc.df <- dplyr::right_join(pc.df, pc.snp, by="sample.id")
pc.df <- dplyr::right_join(pc.df, pData(annot), by="sample.id")
pc.df <- pc.df[match(seqGetData(seqData, "sample.id"), pc.df$sample.id),]
annot <- AnnotatedDataFrame(pc.df)
all.equal(annot$sample.id, seqGetData(gds, "sample.id"))
sampleData(seqData) <- annot
samples_used <- read.csv("pspcon.2.fam", sep=" ", header=F)$V2
pcrel <- readRDS("/pcrel.rds")
# covariance matrix from pcrelate output
grm <- pcrelateToMatrix(pcrel, scaleKin=2, sample.include=samples_used)
seqSetFilter(seqData, sample.id=samples_used)
nullmod <- fitNullModel(seqData, outcome="Outcome", family='binomial',
                        covars=c("Sex", "PCR_info", paste0("PC", 1:5), paste0("SV_PC", 1:5)),
                        cov.mat=grm, verbose=FALSE)

seqResetFilter(seqData)
seqSetFilter(seqData, sample.id=samples_used)
seqSetFilterCond(seqData, maf = 0.01, missing.rate = 0.5, parallel=3)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE, BPPARAM=BiocParallel::MulticoreParam(5))
# add SV information
seqResetFilter(gds)
seqSetFilter(gds, sample.id=samples_used)
seqSetFilterCond(gds, maf = 0.01, missing.rate = 0.5, parallel=2)
alleles <- data.frame(row.names = seqGetData(gds, "variant.id"),seqGetData(gds, 'allele'), seqGetData(gds, 'annotation/id'))
colnames(alleles) <- c("V1", 'V2')
all(assoc$variant.id %in% rownames(alleles))
# must change to character
assoc$alleles <- alleles[as.character(assoc$variant.id), "V1"]
assoc$ref <- sapply(strsplit(assoc$alleles, ","), `[`, 1)
assoc$alt <- sapply(strsplit(assoc$alleles, ","), `[`, 2)
assoc$sv <- alleles[as.character(assoc$variant.id), "V2"]
write.csv(assoc, file = "SV_associations.csv")
