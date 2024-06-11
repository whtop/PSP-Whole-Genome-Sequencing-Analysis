import os
import re
import pandas as pd

from io import StringIO
from pandas_plink import read_plink1_bin

def samplot(names, cram, fasta, chro, start, end, ty, outname, expand=False):
    if expand:  # whether to expand range
        length = int(end) - int(start)
        start = str(int(start) - length)
        end = str(int(end) + length)
    os.system("samplot plot -n {} -b {} -o {} -c {} -s {} -e {} -t {} -r {}".format(
            " ".join(names),
            " ".join(cram),
            outname + ".png",
            chro, start, end,
            ty,
            fasta
        ))
    return

def get_samples(bed, var, p):
    mapper = {":".join(i.split(sep=":")[:-2]):j for i,j in zip(bed.snp.values, bed.variant.values)}
    geno = bed.sel(variant=mapper[var])
    assert geno.a1=="A"
    df = pd.DataFrame(index=geno.sample.values)
    df['geno'] = geno.values
    return df.merge(p, left_index=True, right_on='Sample_ID')


svs = pd.read_csv("./sv_associations.csv")
bed1 = read_plink1_bin("./psp.bed", ref='a0')
df = pd.read_csv("./pheno.csv")

fasta = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
f = open("10264-sample-S3-locations.csv")
next(f)
crams = {re.search(r'/(.*?)_vcpa', i).group(1).split(sep="/")[-1]:i.strip() for i in f if i.strip()}
direc = "bams"
if not os.path.exists(direc): os.mkdir(direc)
for sv in svs.SV:
    chro, pos, ty = sv.split(sep=":")
    st, ed = map(int, pos.split(sep="-"))
    if ty == "INS": continue
    samples = get_samples(bed1, sv, df)

    used = {}
    rg = ed - st
    for i in samples.geno[~samples.geno.isna()].unique():
        ss = samples[samples.geno==i]
        ss = ss.iloc[:min(10, ss.shape[0])]
        used[i] = ss

    # download bam
    for i in used:
        for s in used[i].Sample_ID:
            cram = crams[s]
            crai = os.path.join(direc, os.path.basename(cram) + ".crai")
            outname = "{}_{}.cram".format(crai[:-10], sv)
            if not os.path.exists(crai):
                os.system("aws s3 cp --request-payer requester {}.crai {}".format(cram, crai))
            if (not os.path.exists(outname)) or (not os.path.exists(outname+".crai")):
                url=os.popen("python3 ./s3-presigned-samtools.py {}".format(cram)).read().strip()
                os.system("samtools view -X '{}' {} {}:{}-{} -o {}".format(url, crai, chro, st - max(rg, 10000), ed + max(rg, 10000), outname))
                os.system("samtools index {}".format(outname))

    # plot
    for i in used:
        names = used[i].Sample_ID.values
        files = ["{}_{}.cram".format(os.path.join(direc, os.path.basename(crams[j]))[:-5], sv) for j in names]
        samplot(names, files, fasta, chro, st, ed, ty, "plots/{}_{}".format(sv, i))
