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
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.patches as mpatches
from scipy.stats import norm
import scipy.stats as stats
import numpy as np
import re
import scikit_posthocs as sp

new_day = datetime.datetime.now().strftime("%Y%m%d")
path = "../../analysis/plots/dnds_ratio/"

#run_type = "dry"
run_type = "wet"

organ = "ILE"
#organ = "SPL"

palette2 = {"DR_lifelong":"#D55E00", "AL_lifelong":"#0072B2", "AL_DR16M":"#009E73", "AL_DR20M":"#F0E442"}

gr = ['5AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

# {{{
if organ == "SPL":
    df = pd.read_csv("../../analysis/results_tables/SHM_absolute_SPL.tsv", sep = " ")
else:
    df = pd.read_csv("../../analysis/results_tables/SHM_absolute_ILE.tsv", sep = " ")

#df["Diet"] = df["Diet"].replace({0:"DR_lifelong", 3:"AL_lifelong", 1:"AL_DR16M", 2:"AL_DR20M"})
df = df[df["biogroup"] != "16AL_lifelong"]

# Reading metadata file
metadata = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";", usecols = ["Diet", "Age", "File", "Mouse"])
metadata = metadata.replace(" ", "_", regex = True)
metadata["biogroup"] = metadata["Age"].astype(str) + metadata["Diet"]

df = df.merge(metadata.loc[:, ["Diet", "Age", "biogroup"]], on = "biogroup", how = "inner")
# }}}

df["dnds_ratio"] = df["mu_count_seq_r"] / df["mu_count_seq_s"]
df["dnds_ratio"].replace([np.inf, -np.inf], np.nan, inplace = True)
df.dropna(inplace = True)

df1 = pd.DataFrame(df.groupby(["sample_id", "biogroup", "Age", "Diet"])["dnds_ratio"].mean()).reset_index()

fig, ax = plt.subplots(figsize = (10, 5))
sns.boxplot(data = df1, x = "biogroup", y = "dnds_ratio", palette = palette3, order = gr)

# {{{
# This has standard error of the mean
fig, ax = plt.subplots(figsize = (3, 2.5))
sns.boxplot(data = df1, x = "Age", y = "dnds_ratio", hue = "Diet", palette = palette2, showfliers = False, 
           hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)
sns.swarmplot(data = df1, x = "Age", y = "dnds_ratio", hue = "Diet", dodge = True, color = ".25", 
           hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)

ax.tick_params(axis = "x", labelsize=8)
ax.tick_params(axis = "y", labelsize=8)
ax.set_xlabel("Age [Months]", fontsize = 9)
ax.set_ylabel("dN/dS ratio", fontsize = 9)

ax.get_legend().remove()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

plt.savefig("{}{}_dNdS_ratio_{}.pdf".format(path, organ, new_day))
# }}}

# {{{
dfs_small = df1[df1["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

model = ols('dnds_ratio ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
for e in [5, 20, 24]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = df1[df1["Age"] == e]

    print(stats.kruskal(test["dnds_ratio"].to_list(), test["Diet"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, "dnds_ratio"))
    print(sp.posthoc_mannwhitney(test, val_col = "dnds_ratio", group_col = "Diet", p_adjust = 'fdr_bh'))

print("//K-W BETWEEN timepoints WITHIN treatments//")
for e in ["AL_lifelong", "DR_lifelong", "AL_DR16M"]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = df1[df1["Diet"] == e]

    print(stats.kruskal(test["dnds_ratio"].to_list(), test["Age"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, "dnds_ratio"))
    print(sp.posthoc_mannwhitney(test, val_col = "dnds_ratio", group_col = "Age", p_adjust = 'fdr_bh'))


def boxplots_cond(df1, order, method):
    df2 = df1[df1["biogroup"].isin(order)]

    #Linear regression
    x = df2.Age.to_list()
    y = df2[method].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")


print("Synonymous: AL")
boxplots_cond(df1, ["5AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
             "dnds_ratio")
print("\nSynonymous: DR")
boxplots_cond(df1, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
             "dnds_ratio")
print("\nSynonymous: AL_DR16M")
boxplots_cond(df1, ["20AL_DR16M", "24AL_DR16M"], 
             "dnds_ratio")
print("\nSynonymous: AL_DR20M")
boxplots_cond(df1, ["20AL_lifelong", "24AL_DR20M"], 
             "dnds_ratio")




