#run gcta. 

gcta64=/mnt/analysis-psp/users/ukumar/programs/gcta_1.93.2/gcta64

#setup fam file
outdir=/mnt/analysis-psp/users/timchang/v3/pc/2211/results/hwe/gcta/
mkdir $outdir
cat /mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/pspcon.2.fam | awk '{print $1, $2, $6}' | awk '{ $3 = ($3 == "1" ? 0 : 1) } 1' > ${outdir}pheno.txt

#1. Segment based LD score
chrs=($(seq 1 1 22));
for ichr in "${chrs[@]}";
do
	$gcta64 --bfile /mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/id/plink/r3.bm.chr${ichr}.PspNhwCon.id --ld-score-region 200 --ld-wind 200 --thread-num 12 --out ${outdir}r3.bm.chr${ichr}.PspNhwCon.id.out;
done;
#combine #outfiles. remove headers from 2-22
ldallout=${outdir}r3.bm.chrALLauto.PspNhwCon.id.out.score.ld
cat ${outdir}r3.bm.chr1.PspNhwCon.id.out.score.ld > ${outdir}r3.bm.chrALLauto.PspNhwCon.id.out.score.ld
chrs=($(seq 2 1 22))
for ichr in "${chrs[@]}"; do             
	tail -n +2 ${outdir}r3.bm.chr${ichr}.PspNhwCon.id.out.score.ld >> ${outdir}r3.bm.chrALLauto.PspNhwCon.id.out.score.ld;
done

#2. finds 4 quartiles of SNPs based on LD score
Rscript /mnt/analysis-psp/users/timchang/v3/pc/2211/code/gcta/GCTA_LDMS_chrALL_step2.R

#3. making GRMs using SNPs stratified into different groups with MAF<0.01, 0.1-0.2, 0.2-0.5
plink_file=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/vcfs/bm/eur/hwe/id/plink/r3.bm.chrALLauto.PspNhwCon.id;
maxmafs=( 0.01 0.1 0.2 0.5 );
minmafs=(0 0.01 0.1 0.2);
length=${#maxmafs[@]};

#write file names to common+rare and common grm
crmultigrm=${outdir}common_rare_grm.txt;
cmultigrm=${outdir}common_grm.txt;
>${crmultigrm};
>${cmultigrm};
for (( i=0; i<${length}; i++ ));
do
	echo ${maxmafs[$i]};
	#echo ${minmafs[$i]}

	$gcta64 --bfile $plink_file --maf ${minmafs[$i]} --max-maf ${maxmafs[$i]} --extract ${outdir}snp_group1.txt --make-grm --thread-num 10 --out ${outdir}test_group1_${maxmafs[$i]}; 
	$gcta64 --bfile $plink_file --maf ${minmafs[$i]} --max-maf ${maxmafs[$i]} --extract ${outdir}snp_group2.txt --make-grm --thread-num 10 --out ${outdir}test_group2_${maxmafs[$i]};
	$gcta64 --bfile $plink_file --maf ${minmafs[$i]} --max-maf ${maxmafs[$i]} --extract ${outdir}snp_group3.txt --make-grm --thread-num 10 --out ${outdir}test_group3_${maxmafs[$i]};
	$gcta64 --bfile $plink_file --maf ${minmafs[$i]} --max-maf ${maxmafs[$i]} --extract ${outdir}snp_group4.txt --make-grm --thread-num 10 --out ${outdir}test_group4_${maxmafs[$i]};

	echo "${outdir}test_group1_${maxmafs[$i]}" >> ${crmultigrm}; 
	echo "${outdir}test_group2_${maxmafs[$i]}" >> ${crmultigrm};
	echo "${outdir}test_group3_${maxmafs[$i]}" >> ${crmultigrm}; 
	echo "${outdir}test_group4_${maxmafs[$i]}" >> ${crmultigrm}; 
	
	echo "${outdir}test_group2_${maxmafs[$i]}" >> ${cmultigrm}; 
	echo "${outdir}test_group3_${maxmafs[$i]}" >> ${cmultigrm}; 
	echo "${outdir}test_group4_${maxmafs[$i]}" >> ${cmultigrm}; 

done

#4. REML analysis with multiple GRM
covfile=/mnt/analysis-psp/users/timchang/v3/pc/2211/data/clin/mod/pspcon.hwe.sex.cov #covariate.txt
$gcta64 --mgrm ${outdir}common_grm.txt --pheno ${outdir}pheno.txt --covar $covfile --prevalence 0.00005 --reml --reml-no-constrain --thread-num 12 --out ${outdir}REML_common;
$gcta64 --mgrm ${outdir}common_rare_grm.txt --pheno ${outdir}pheno.txt --covar $covfile --prevalence 0.00005 --reml --reml-no-constrain --thread-num 12 --out ${outdir}REML_common_rare;


