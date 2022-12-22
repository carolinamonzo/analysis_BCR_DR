# ---
# jupyter:
#   jupytext:
#     cell_markers: '{{{,}}}'
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# {{{
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import matplotlib
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm
import matplotlib.patches as mpatches
from scipy.stats import norm
import scipy.stats as stats
import numpy as np
import re
import scikit_posthocs as sp

new_day = datetime.datetime.now().strftime("%Y%m%d")
path = "../../analysis/VJusage/"

#run_type = "dry"
run_type = "wet"

organ = "ILE"
#organ = "SPL"

gene = "V"
#gene = "J"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", "AL_DR16M":"teal", "AL_DR20M":"gold"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]

# Reading metadata file
metadata = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";", usecols = ["Diet", "Age", "File", "Mouse"])
metadata = metadata.replace(" ", "_", regex = True)
metadata["Biogroup"] = metadata["Age"].astype(str) + metadata["Diet"]

# Read file with frequencies
df = pd.read_csv(f"../../analysis/results_tables/{gene}gene_usage_{organ}.tsv", sep = " ")

# {{{
# Unifing subtypes of gene into original one

## DONT DO FOR ALLELES
df["gene"] = [x.split("S")[0] for x in df["gene"]]
if organ == "ILE":
    df["gene"] = [x.split("-")[0] for x in df["gene"]]
# }}}

def make_biorep_heatmap(freq, figsize):
    
    freq = freq.reindex(index=sorted(freq.index, key = natsort))
    fig, ax = plt.subplots(figsize = figsize)
    sns.heatmap(freq.T, cmap="Blues", cbar = False, linewidths=.5)
    matplotlib.rcParams['pdf.fonttype'] = 42 
    plt.tight_layout()


## Remove duplicate entries to the same gene
df = df.groupby(["sample_id", "gene"]).aggregate({"seq_count":"sum", "seq_freq":"sum"}).reset_index()

df_matrix = df.loc[:, ["sample_id", "gene", "seq_freq"]].fillna(0).pivot(index = "sample_id", columns = "gene")
df_matrix.columns = df_matrix.columns.get_level_values(1)
df_matrix.fillna(0, inplace = True)
df_matrix = df_matrix[sorted(df_matrix.columns, key = natsort)]

# Keep only metadata of the mice of our tissue
metadata = metadata[metadata["File"].isin(df_matrix.index)]

fig, ax = plt.subplots(figsize = (20, 20))
sns.heatmap(df_matrix, cmap = "Blues")

for e in gr:

    samps = metadata[metadata["Biogroup"] == e]["File"]
    if gene == "V":
        fig, ax = plt.subplots(figsize = (15, 3))
    else:
        fig, ax = plt.subplots(figsize = (5, 3))
    sns.heatmap(df_matrix.loc[samps, :], cmap = "Blues", cbar = False, linewidths = .5)
    plt.title(e)
    matplotlib.rcParams['pdf.fonttype'] = 42 
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig(f"{path}/{organ}_{gene}_{e}_biorep.pdf")
    else:
        plt.show()

mer = df_matrix.merge(metadata.loc[:, ["File", "Biogroup"]], right_on = "File", left_index = True)
mer.drop(columns = "File", inplace = True)

# {{{
biog = mer.groupby("Biogroup").median()
biog = biog.loc[gr, :]

if gene == "V":
    g = sns.clustermap(biog, cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 4))
else:
    g = sns.clustermap(biog, cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 4))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_ageDiet_all.pdf")
else:
    plt.show()
# }}}

if gene == "V":
    g = sns.clustermap(biog.loc[["5AL_lifelong", "20AL_lifelong", "24AL_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 2.5))
else:
    g = sns.clustermap(biog.loc[["5AL_lifelong", "20AL_lifelong", "24AL_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 2.5))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_AL.pdf")
else:
    plt.show()

if gene == "V":
    g = sns.clustermap(biog.loc[["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 2.5))
else:
    g = sns.clustermap(biog.loc[["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 2.5))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_DR.pdf")
else:
    plt.show()

if gene == "V":
    g = sns.clustermap(biog.loc[["5AL_lifelong", "20AL_DR16M", "24AL_DR16M"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 2.5))
else:
    g = sns.clustermap(biog.loc[["5AL_lifelong", "20AL_DR16M", "24AL_DR16M"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 2.5))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_ALDR16.pdf")
else:
    plt.show()

if gene == "V":
    g = sns.clustermap(biog.loc[["5AL_lifelong", "20AL_lifelong", "24AL_DR20M"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 2.5))
else:
    g = sns.clustermap(biog.loc[["5AL_lifelong", "20AL_lifelong", "24AL_DR20M"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 2.5))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_ALDR20.pdf")
else:
    plt.show()

if gene == "V":
    g = sns.clustermap(biog.loc[["5AL_lifelong", "5DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 1.5))
else:
    g = sns.clustermap(biog.loc[["5AL_lifelong", "5DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 1.5))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_5M.pdf")
else:
    plt.show()

if gene == "V":
    g = sns.clustermap(biog.loc[["20AL_lifelong", "20AL_DR16M", "20DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 2))
else:
    g = sns.clustermap(biog.loc[["20AL_lifelong", "20AL_DR16M", "20DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 2))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_20M.pdf")
else:
    plt.show()

if gene == "V":
    g = sns.clustermap(biog.loc[["24AL_lifelong", "24AL_DR20M", "24AL_DR16M", "24DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (15, 2.5))
else:
    g = sns.clustermap(biog.loc[["24AL_lifelong", "24AL_DR20M", "24AL_DR16M", "24DR_lifelong"], :], cmap = "Blues", linewidths = .5, 
               row_cluster = False, figsize = (5, 2.5))
g = g.cax.set_visible(False)
if run_type != "dry":
    plt.savefig(f"{path}/{organ}_{gene}_24M.pdf")
else:
    plt.show()

merall = df_matrix.merge(metadata, left_index = True, right_on = "File")
merall.drop(columns = ["Mouse", "File", "Biogroup"], inplace = True)

print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
for e in [5, 20, 24]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = merall[merall["Age"] == e]
    
    if gene == "V":          
        for j in range(0, 16):
            print(merall.columns.to_list()[j])
            print("Kruskall")
            print(stats.kruskal(test[merall.columns.to_list()[j]].to_list(), test["Diet"].to_list()))
            print("Mann-Whitney")
            print(sp.posthoc_dunn(test, val_col = merall.columns.to_list()[j], group_col = "Diet", p_adjust = 'bonferroni'))
            print("\\n\\")
    else:
        for j in range(0, 4):
            print(merall.columns.to_list()[j])
            print("Kruskall")
            print(stats.kruskal(test[merall.columns.to_list()[j]].to_list(), test["Diet"].to_list()))
            print("Mann-Whitney")
            print(sp.posthoc_dunn(test, val_col = merall.columns.to_list()[j], group_col = "Diet", p_adjust = 'bonferroni'))
            print("\\n\\")


def boxplots_cond(df_prov, d, method):
    
    # Subset the info we need
    df = df_prov[df_prov["Diet"] == d]
    if d == "AL_lifelong":
        df = df[df["Age"] != 16]
    
    #Linear regression
    x = df.Age.to_list()
    y = df[method].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\np-value: {p_value}\n")


# {{{
print("\nAL")

if gene == "V":
    for j in range(0, 16):
        print(merall.columns.to_list()[j])
        boxplots_cond(merall, "AL_lifelong", merall.columns.to_list()[j])
else:
    for j in range(0, 4):
        print(merall.columns.to_list()[j])
        boxplots_cond(merall, "AL_lifelong", merall.columns.to_list()[j])
# }}}

# {{{
print("\nDR")

if gene == "V":
    for j in range(0, 16):
        print(merall.columns.to_list()[j])
        boxplots_cond(merall, "DR_lifelong", merall.columns.to_list()[j])
else:
    for j in range(0, 4):
        print(merall.columns.to_list()[j])
        boxplots_cond(merall, "DR_lifelong", merall.columns.to_list()[j])
# }}}
# {{{
merall = merall[merall["Age"] != 16]

for e in merall.columns[:-2]:
    print(e)
    dfs_small = merall[merall["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

    model = ols('{} ~ C(Age) + C(Diet) + C(Age):C(Diet)'.format(e), data = dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))
# }}}















