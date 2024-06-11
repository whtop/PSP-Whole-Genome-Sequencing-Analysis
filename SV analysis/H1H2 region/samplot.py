import os
import re
import pandas as pd

from io import StringIO

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

svs = pd.read_csv("SVs.txt", sep="\t", index_col=0)
fasta = "~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
f = open("./10264-sample-S3-locations.csv")
next(f)
crams = {re.search(r'/(.*?)_vcpa', i).group(1).split(sep="/")[-1]:i.strip() for i in f if i.strip()}
direc = "bams"
if not os.path.exists(direc): os.mkdir(direc)
for variant, sv in svs.iterrows():
    alt = re.search(r'(<.*>)', variant).group(1)
    variant = variant.replace(">", "").replace("<","")  # to avoild error in using "<>" in commands
    if "INS" in variant: continue

    output  = os.popen("bcftools view {}.reheader.vcf.gz -r {}".format(
        variant.split(sep=":")[0][3:], ":".join(variant.split(sep=":")[:2])
        )).read()
    output = "\n".join([i for i in output.split(sep="\n") if not i.startswith("##")])
    df = pd.read_csv(StringIO(output), sep="\t")
    df = df[df.POS == int(variant.split(sep=":")[1].split(sep="-")[0])]
    print(df.iloc[:, :9])
    df = df[[alt==i for i in df.ALT]]  # choose the unique one by MODEL (AGG, BP, CO)
    
    chro, pos, ty = variant.split(sep=":")[:3]
    st, ed = map(int, pos.split(sep="-"))
    rg = ed - st

    used = []
    k = 0   # plot for homo reference
    for i, j in df.iloc[0,9:].iteritems():
        if j.startswith('0/0') and ("FAIL" not in j):
            k += 1
            cram = crams[i]
            crai = os.path.join(direc, os.path.basename(cram) + ".crai")
            outname = "{}_{}.cram".format(crai[:-10], variant)
            if not os.path.exists(crai):
                os.system("aws s3 cp --request-payer requester {}.crai {}".format(cram, crai))
            if (not os.path.exists(outname)) or (not os.path.exists(outname+".crai")):
                url=os.popen("python3 ./s3-presigned-samtools.py {}".format(cram)).read().strip()
                os.system("samtools view -X '{}' {} {}:{}-{} -o {}".format(url, crai, chro, st - max(rg, 10000), ed + max(rg, 10000), outname))
                os.system("samtools index {}".format(outname))
            used.append(i)
            if k >= 5: break
    if used:
        files = ["{}_{}.cram".format(os.path.join(direc, os.path.basename(crams[j]))[:-5], variant) for j in used]
        if not os.path.exists("plots/{}_0.png".format(variant)):
            samplot(used, files, fasta, chro, st, ed, ty, "plots/{}_0".format(variant))
        else:
            print("plots/{}_0.png".format(variant))

    used = []
    k = 0   # plot for heter
    for i, j in df.iloc[0,9:].iteritems():
        if (j.startswith('1/0') or j.startswith('0/1')) and ("FAIL" not in j):
            k += 1
            cram = crams[i]
            crai = os.path.join(direc, os.path.basename(cram) + ".crai")
            outname = "{}_{}.cram".format(crai[:-10], variant)
            if not os.path.exists(crai):
                os.system("aws s3 cp --request-payer requester {}.crai {}".format(cram, crai))
            if (not os.path.exists(outname)) or (not os.path.exists(outname+".crai")):
                url=os.popen("python3 ./s3-presigned-samtools.py {}".format(cram)).read().strip()
                os.system("samtools view -X '{}' {} {}:{}-{} -o {}".format(url, crai, chro, st - max(rg, 10000), ed + max(rg, 10000), outname))
                os.system("samtools index {}".format(outname))
            used.append(i)
            if k >= 5: break
    if used:
        files = ["{}_{}.cram".format(os.path.join(direc, os.path.basename(crams[j]))[:-5], variant) for j in used]
        if not os.path.exists("plots/{}_1.png".format(variant)):
            samplot(used, files, fasta, chro, st, ed, ty, "plots/{}_1".format(variant))
        else:
            print("plots/{}_1.png".format(variant))

    used = []
    k = 0   # plot for homo SV
    for i, j in df.iloc[0,9:].iteritems():
        if j.startswith('1/1') and ("FAIL" not in j):
            k += 1
            cram = crams[i]
            crai = os.path.join(direc, os.path.basename(cram) + ".crai")
            outname = "{}_{}.cram".format(crai[:-10], variant)
            if not os.path.exists(crai):
                os.system("aws s3 cp --request-payer requester {}.crai {}".format(cram, crai))
            if (not os.path.exists(outname)) or (not os.path.exists(outname+".crai")):
                url=os.popen("python3 ./s3-presigned-samtools.py {}".format(cram)).read().strip()
                os.system("samtools view -X '{}' {} {}:{}-{} -o {}".format(url, crai, chro, st - max(rg, 10000), ed + max(rg, 10000), outname))
                os.system("samtools index {}".format(outname))
            used.append(i)
            if k >= 5: break
    if used:
        files = ["{}_{}.cram".format(os.path.join(direc, os.path.basename(crams[j]))[:-5], variant) for j in used]
        if not os.path.exists("plots/{}_2.png".format(variant)):
            samplot(used, files, fasta, chro, st, ed, ty, "plots/{}_2".format(variant))
        else:
            print("plots/{}_2.png".format(variant))
