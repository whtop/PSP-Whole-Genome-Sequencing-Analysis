import os
for i in range(22):
    os.system("sbatch --job-name=chr{} -o chr{}.out -e chr{}.err --ntasks=1 --qos=large --mem=4g --export=ALL,VCF=gcad.qc.r3.wgs.16905.GATK.2021.08.24.biallelic.genotypes.chr{}.ALL.vcf.gz filter_vcf.sh".format(
        i+1, i+1, i+1, i+1
    ))
