# chr17:45309498-46903773 (about 1.5Mb) refer to alternative contig in the UCSC genome browser
# bcftools view bcftools view 17.reheader.vcf.gz -r chr17:45309498-46903773 -Oz -o h1.vcf.gz

import re
import gzip
import pandas as pd

from io import StringIO

df_s = ""
for i in gzip.open("./h1.vcf.gz"):
    i = i.decode()
    if i.startswith("##"):
        continue
    df_s += i

df = pd.read_csv(StringIO(df_s), sep="\t")
df.index = df.ID.values
print(df.iloc[:, :9])

w = open('SVs.txt', 'w')
w.write("SV")
for i in df.columns[9:]:
    w.write("\t{}".format(i))
w.write("\n")

for _, row in  df.iterrows():
    ed = re.search(r'END=(.*?);', row.INFO)
    if ed:
        ed = ed.group(1)
        w.write("{}:{}-{}:{}".format(row['#CHROM'], row['POS'],
        ed, row.ALT))
        for i in row.values[9:]:
            geno, passed = i.split(sep=":")[:2]
            if passed == "PASS":
                w.write("\t{}".format(int(geno[0])+int(geno[-1])))
            else:
                w.write("\tNA")
        w.write("\n")
w.close()
