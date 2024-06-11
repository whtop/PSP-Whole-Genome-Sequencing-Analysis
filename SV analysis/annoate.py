import os
import numpy as np
import pandas as pd
import pyranges as pr

def annotate(df, genes):
    gr = pr.PyRanges(df)
    return gr.nearest(genes)

def fdr(p, n=None):
    """Give a numpy array, convert P values to FDRs"""

    if not n:
        n = len(p)
    p = np.array(p)
    # do not allow nan values
    assert pd.isna(p).any()==False
    lp = len(p)
    assert n >= lp

    if n <= 1 : 
        return(p)
    else:
        i = np.arange(lp, 0, -1)
        o = np.argsort(-p)   # using -p to sort in decreasing order
        ro = np.argsort(o)
        cummin = pd.Series(n/i * p[o]).cummin().values
        return np.minimum(1, cummin)[ro]



infile = "./sv_associations.csv"
outprefix = os.path.basename(infile)
df = pd.read_csv(infile, index_col=0)
df.rename(columns={"chr":"Chromosome", "sv":"ID"}, inplace=True)
df["Start"] = [int(i.split(sep=":")[1].split(sep="-")[0]) for i in df.ID]
df["End"] = [int(i.split(sep=":")[1].split(sep="-")[1]) for i in df.ID]
df["SVType"] = [i.split(sep=":")[-3] for i in df.ID]
df["Source"] = [i.split(sep=":")[-2] for i in df.ID]
df["AF"] = [i.split(sep=":")[-1] for i in df.ID]
df["Name"] = ["chr{}:{}-{}:{}".format(i,j,m,n) for i,j,m,n in zip(df.Chromosome, df.Start, df.End, df.SVType)]
df['FDR'] = fdr(df['Score.pval'].fillna(1).values)
df = df.sort_values(by="Score.pval")
print(df.shape)
df.to_csv(outprefix + ".annotated", sep=",", index=False)
df = df[df['FDR']<0.05]

print(df.shape)
genes = pd.read_csv("ensemble_hg38_biomart_all_genes_chr+MT_20201003.txt", sep="\t")
genes = genes[["Gene stable ID", "Gene name", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"]]
genes = genes[~genes.duplicated()]  # id is unique, while name is not
genes.index = genes["Gene stable ID"].values
genes.columns = ["ENSG", "SYMBOL", "Chromosome", "Start", "End"]
genes = pr.PyRanges(genes)
dfg = annotate(df, genes)
dfg = dfg.as_df()
dfg = dfg.sort_values(by="Score.pval")

# add annotation
annosv = pd.read_csv("./annoSV/adsp.tsv", sep="\t")
annosv = annosv[annosv.Annotation_mode=="full"]
annosv = annosv[['name', 'AnnotSV_ranking_score', 'AnnotSV_ranking_criteria', 'ACMG_class']]
annosv.index = annosv.name.values
vep = pd.read_csv("./vep/SVs.csv", sep="|", index_col=0)
dfg['vep'] = vep.loc[dfg.Name.values, "most_severe_consequence"].values
dfg['AnnotSV_ranking_score'] = annosv.loc[dfg.Name.values, "AnnotSV_ranking_score"].values
dfg['ACMG_class'] = annosv.loc[dfg.Name.values, "ACMG_class"].values
dfg.to_csv(outprefix+".sig", index=False)
