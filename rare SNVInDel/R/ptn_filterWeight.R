#create weight of variants based on MAF
MAFcut=0.0001
rareCut=0.001

dirtail='PTV.misCADD15.Poly.01NA_I.1'

dir=paste0('/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/vcffilter/',dirtail,'/')
inhead='PSP.chr'
intail='.filter.out'
outtail='.c.vcffilter.out'
wtouttail='.weight.out'
assocdir='/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/plink/assoc.cn/'
assochead='log.chr'
assoctail='.assoc'
famfile='/mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/pspcon.2.fam'

famdat=read.table(famfile,stringsAsFactors=FALSE) 
famtab=table(famdat$V6)
caseprob=famtab['2']/(famtab['2']+famtab['1'])
chr_list=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
alldat=data.frame()
for (chr in chr_list){
	infile=paste0(dir,inhead,as.character(chr),intail);
	associnfile=paste0(assocdir,assochead,as.character(chr),assoctail,'.rda');

	chrdat=read.table(infile,stringsAsFactors=FALSE,col.names=c('gene','var','MAF','VEP'))
	
	chrdatnodup = chrdat[!duplicated(paste(chrdat$gene,chrdat$var)),]
	chrdatnodup$normed=chrdatnodup$MAF/MAFcut
	chrdatnodup$weight=1-chrdatnodup$normed
	chrdatnodup$weightbeta=dbeta(chrdatnodup$normed,1,4)

	load(associnfile)
	chrassdat=chrassoc 
	chrassfdat=merge(chrdatnodup,chrassdat,by.x='var',by.y='SNP',keep.x=TRUE)
	chrassfdat$CCAF=(chrassfdat$C_A+chrassfdat$C_U)/(2*nrow(famdat))
	chrassfdat=chrassfdat[chrassfdat$CCAF!=0 & chrassfdat$CCAF<rareCut,] #must be present in case/control, and case + control allele frequency must be less than rareCut. Don't have enough samples to actually get below MAFcut
	alldat=rbind(alldat,chrassfdat)
	
	filteroutfile=paste0(dir,inhead,as.character(chr),outtail)
	write.table(chrassfdat[,c('gene','var')],file=filteroutfile,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

	weightoutfile=paste0(dir,inhead,as.character(chr),wtouttail)	
	write.table(chrassfdat[,c('var','weightbeta')],file=weightoutfile,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
}

filteroutfile=paste0(dir,inhead,'ALL',outtail)
write.table(alldat[,c('gene','var')],file=filteroutfile,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

filteroutfile=paste0(dir,inhead,'ALL.MAF',outtail)
write.table(alldat[,c('gene','var','MAF')],file=filteroutfile,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

wtout=alldat[!duplicated(alldat$var),]
weightoutfile=paste0(dir,inhead,'ALL','.BETA',wtouttail)
write.table(wtout[,c('var','weightbeta')],file=weightoutfile,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

outfile=paste0(dir,inhead,'ALL','.BETA',wtouttail,'.rda')
save('alldat',file=outfile)



