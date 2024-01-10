#analyze data
#conda activate /mnt/analysis-psp/users/timchang/programs/miniconda/R3.6.3
#conda install r-ggforce rcolorbrewer

.libPaths('/mnt/analysis-psp/users/timchang/programs/R/3.6.3/')
library(ggforce);library(ggplot2);library(RColorBrewer)
options(stringsAsFactors=F)
geno_path='/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/m1000G_20181203/maf.01/plink/merge/'
geno_file=paste0(geno_path,'r3.bm.chrALLauto.PspNhwCon.1KG.maf.01.s.namb.g')
pca_infile=paste0(geno_path,'r3.bm.chrALLauto.PspNhwCon.1KG.maf.01.s.namb.g.projpca','.rds')
load(pca_infile)
#'scores','obj.bed.autoSVD','fam'
#get finer population
eurpop_file='/mnt/analysis-psp/users/timchang/v3/pc/2211/data/ref/1000G/20181203/integrated_call_samples_v3.20200731.ALL.20181203.eurpop'
eurpop=read.table(eurpop_file,header=T)

totalpc=20
casename='Case'
conname='Control'
eurname='EUR'
numpc=6
numsd=8
mypal='Set3'
outOrigPlot=paste0(geno_file,'.projpca.',numpc,'.',numsd,'.orig.pdf')
outOutList=paste0(geno_file,'.projpca.',numpc,'.',numsd,'.outlier.txt')
outpedfile=paste0(geno_file,'.projpca.',numpc,'.',numsd,'.out.pedind')
outpopfile=paste0(geno_file,'.projpca.',numpc,'.',numsd,'.out.poplist')
outkeepfile=paste0(geno_file,'.projpca.',numpc,'.',numsd,'.keep.txt')

#convert scores to df. add columns (in same order)
pcadat=as.data.frame(scores$simple_proj)
names(pcadat)=paste0('PC',1:dim(pcadat)[2])
pcadat$sample.ID=fam$sample.ID
pcadat$diagnosis=fam$affection
pcadat$pop=fam$family.ID
pcadat$pop[pcadat$diagnosis==1]='Control'
pcadat$pop[pcadat$diagnosis==2]='Case'

#add fine grain eur
pcadat$eurpop=pcadat$pop
indeur=match(eurpop$Individual.ID,pcadat$sample.ID)
pcadat[indeur,'eurpop'] = eurpop$Population2

table(pcadat$pop)
table(pcadat$eurpop)

cat('pdf name: ', outOrigPlot,'\n')
pdf(outOrigPlot)

for (i in c(1,3,5,7,9)){
   iplot = ggplot(pcadat, aes_string(x=paste0('PC',i),y=paste0('PC',i+1))) + geom_point(aes(color=pop),alpha=0.6,size=1) + scale_color_brewer(palette=mypal)
   print(iplot)
}
for (i in c(1,3,5,7,9)){
   iplot = ggplot(pcadat, aes_string(x=paste0('PC',i),y=paste0('PC',i+1))) + geom_point(aes(color=eurpop),alpha=0.6,size=1) + scale_color_brewer(palette=mypal)
   print(iplot)
}

tfcasecon=pcadat$pop==casename | pcadat$pop==conname
tfeur = pcadat$pop==eurname
pcacol=paste0('PC',1:numpc)
eurmean=colMeans(pcadat[tfeur,pcacol])
eursd=apply(pcadat[tfeur,pcacol],2,sd)

pcadatd = pcadat[,pcacol]

#within numsd ellipsoid. (x-h)^2/r^2 <=1
pcadatdell=sweep(sweep(pcadatd,2,eurmean,'-')^2,2,(numsd*eursd)^2,'/')
pcadatdelld=rowSums(pcadatdell)
tfout=pcadatdelld>1
write.table(file=outOutList,pcadat[tfout & tfcasecon,'sample.IDp'],row.names=F,col.names=F,quote=F)

#make outlier category
pcadat$popout=as.character(pcadat$pop)
pcadat[tfout & tfcasecon, 'popout'] = paste0(pcadat[tfout & tfcasecon,'pop'],'-out')
pcadat$popout=as.factor(pcadat$popout)
print(table(pcadat$popout))
indeurpal=which(levels(pcadat$popout)==eurname)
coleurpal=brewer.pal(length(levels(pcadat$popout)),mypal)[indeurpal]
for (i in c(1,3,5,7,9)){
   iplot = ggplot(pcadat, aes_string(x=paste0('PC',i),y=paste0('PC',i+1))) + geom_point(aes(color=popout),alpha=0.75,size=1) + scale_color_brewer(palette=mypal)+geom_ellipse(aes(x0 = eurmean[i], y0=eurmean[i+1], a=numsd*eursd[i],b=numsd*eursd[i+1],angle=0),color=coleurpal)
   print(iplot)
}

#plot outlier without ellipse
for (i in c(1,3,5,7,9)){
   iplot = ggplot(pcadat, aes_string(x=paste0('PC',i),y=paste0('PC',i+1))) + geom_point(aes(color=popout),alpha=0.75,size=1) + scale_color_brewer(palette=mypal)
   print(iplot)
}

dev.off()

#write sample.ID to keep
write.table(file=outkeepfile, pcadat[tfcasecon & !tfout, 'sample.ID'],row.names=F,col.names=F,quote=F)
