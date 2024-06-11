import numpy as np
import pandas as pd
import pyranges as pr


def annotate(df, genes):
    gr = pr.PyRanges(df)
    return gr.nearest(genes)

df = pd.read_csv("associations.csv", index_col=0)
df['Chromosome'] = df.chr.values
df["Start"] = df.pos.values
df["End"] = df.pos.values

genes = pd.read_csv("ensemble_hg38_biomart_all_genes_chr+MT_20201003.txt", sep="\t")
genes = genes[["Gene stable ID", "Gene name", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"]]
genes = genes[~genes.duplicated()]  # id is unique, while name is not
genes.index = genes["Gene stable ID"].values
genes.columns = ["ENSG", "SYMBOL", "Chromosome", "Start", "End"]
genes = pr.PyRanges(genes)
dfg = annotate(df, genes)
dfg = dfg.as_df()
dfg.to_csv("snvs.csv", index=False)
