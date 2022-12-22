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
path = "../../analysis/plots/naiveClassSwitch/"

run_type = "dry"
#run_type = "wet"

organ = "ILE"
#organ = "SPL"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", "AL_DR16M":"teal", "AL_DR20M":"gold"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

# Read metadata so we can get age and diet from there
metadata = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8"])
metadata.columns = ["Mouse", "Diet", "Age", "Illumina", "Illumina2", "Barcode", "sample_id"]

# Read full file to get counts
df = pd.read_csv("../../analysis/created_germlines/reannotated_db/merged_changeo_{}_isosum.tsv".format(organ), sep = "\t")

df.head()

shm = pd.read_csv("../../analysis/results_tables/SHM_frequency_{}.tsv".format(organ), sep = " ")

dfshm = pd.merge(df, shm.loc[:, ["sequence_id", "sample_id", "mu_freq_seq_r", "mu_freq_seq_s"]], on = ["sequence_id", "sample_id"], how = "inner")

# {{{
sample = []
IgMnIgDn = []
IgMpIgDpSHMp = []
IgMpIgDpSHMn = []

for samp in dfshm["sample_id"].unique():
    sample.append(samp)

    test = dfshm[dfshm["sample_id"] == samp]
    test = test.loc[:, ["sample_id", "sum_iso", "clone_id", "mu_freq_seq_r", "mu_freq_seq_s"]]
    total = test.shape[0]

    # To extract full switched, get clone id of the clones that have a naive, and then remove all those clones.
    list_CSclones = list(test[test["sum_iso"].isin(["IGD", "IGM"])]["clone_id"].unique())
    test_CS = test[~test["clone_id"].isin(list_CSclones)]
    post_antigen = test_CS.shape[0]
    post_antigen = float(post_antigen) / float(total)
    IgMnIgDn.append(post_antigen)

    #From the ones that are not full switched, get naive and antigen exposed from SHM frequency
    rest = test[test["clone_id"].isin(list_CSclones)]
    rest["total_mu_freq"] = rest["mu_freq_seq_r"] + rest["mu_freq_seq_r"]

    naive = rest[rest["total_mu_freq"] < 0.01]
    antigen = rest[rest["total_mu_freq"] >= 0.01]

    naive = naive.shape[0]
    naive = float(naive)/float(total)
    IgMpIgDpSHMn.append(naive)
    antigen = antigen.shape[0]
    antigen = float(antigen)/float(total)
    IgMpIgDpSHMp.append(antigen)
# }}}

CS_df = pd.DataFrame()
CS_df["sample_id"] = sample
CS_df["IgM-IgD-"] = IgMnIgDn
CS_df["IgM+IgD+SHM+"] = IgMpIgDpSHMp
CS_df["IgM+IgD+SHM-"] = IgMpIgDpSHMn

CS_df_met = pd.merge(CS_df, metadata.loc[:, ["Mouse", "Diet", "Age", "sample_id"]], on = "sample_id", how = "inner")
CS_df_met["biogroup"] = CS_df_met["Age"].astype(str) + "M " + CS_df_met["Diet"].astype(str)
CS_df_met.drop(columns = ["Age", "Diet", "sample_id"], inplace = True)

CS_df_met_mean = CS_df_met.drop(columns = ["Mouse"]).groupby("biogroup").mean()
CS_df_met_mean.drop("16M AL lifelong", axis = 0, inplace = True)

# {{{
CS_df_met_mean = CS_df_met_mean.reindex(['5M AL lifelong', '5M DR lifelong', '20M AL lifelong', '20M DR lifelong', '20M AL_DR16M', '24M AL lifelong',
       '24M DR lifelong', '24M AL_DR16M', '24M AL_DR20M'])
fig, ax = plt.subplots(figsize = (6, 4))
ax = sns.heatmap(CS_df_met_mean, vmin = 0, vmax = 1, cmap = "Greens", linewidths = .5)
ax.set_ylabel('') 
ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 10)
ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize = 11)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig(f"{path}{organ}_all_Heat_{new_day}.pdf")
else:
    plt.plot()
# }}}

CS_df_met_mean


def plot_heats(df, l, s, siz):
    df = df[df.index.isin(l)]
    fig, ax = plt.subplots(figsize = siz)
    ax = sns.heatmap(df, vmin = 0, vmax = 1, cmap = "Greens", linewidths = .5)
    ax.set_ylabel('') 
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 10)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize = 11)
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig(f"{path}{organ}_{s}_Heat_{new_day}.pdf")
    else:
        plt.plot()


plot_heats(CS_df_met_mean, ["5M AL lifelong", "5M DR lifelong"], "5MALDR", (6, 1.4))
plot_heats(CS_df_met_mean, ["20M AL lifelong", "20M AL_DR16M", "20M DR lifelong"], "20MALDR16", (6, 1.8))
plot_heats(CS_df_met_mean, ["24M AL lifelong", "24M AL_DR20M", "24M AL_DR16M", "24M DR lifelong"], "24MALDR1620", (6, 2.2))

plot_heats(CS_df_met_mean, ["5M AL lifelong", "20M AL lifelong", "24M AL lifelong"], "AgeAL", (6, 1.8))
plot_heats(CS_df_met_mean, ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"], "AgeDR", (6, 1.8))
plot_heats(CS_df_met_mean, ["20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"], "AgeSwitch", (6, 1.8))

CS_df_met.to_csv(f"{path}{organ}_NaiveClassw_{new_day}.csv", sep = ";", index = False)

CS_df_met2 = CS_df_met.merge(metadata.loc[:, ["Mouse", "Diet", "Age"]], on = "Mouse", how = "inner").drop_duplicates()
CS_df_met2.head()


def linreg_age(df, order, co):
    x = df[df["biogroup"].isin(order)].sort_values(by = ["Age"], ascending = True)["Age"].to_list()
    y = df[df["biogroup"].isin(order)].sort_values(by = ["Age"], ascending = True)[co].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print(co)
    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")


linreg_age(CS_df_met2, ["5M AL lifelong", "20M AL lifelong", "24M AL lifelong"], "IgM-IgD-")
linreg_age(CS_df_met2, ["5M AL lifelong", "20M AL lifelong", "24M AL lifelong"], "IgM+IgD+SHM+")
linreg_age(CS_df_met2, ["5M AL lifelong", "20M AL lifelong", "24M AL lifelong"], "IgM+IgD+SHM-")

linreg_age(CS_df_met2, ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"], "IgM-IgD-")
linreg_age(CS_df_met2, ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"], "IgM+IgD+SHM+")
linreg_age(CS_df_met2, ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"], "IgM+IgD+SHM-")

# {{{
from scipy import stats

def KW_bwTreat(CS_df_met2, method):
    print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS: {}//".format(method))
    for e in [5, 20, 24]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = CS_df_met2[CS_df_met2["Age"] == e]

        print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))
        
def KW_bwTime(d, method):
    print("//K-W BETWEEN timepoints WITHIN treatments: {}//".format(method))
    for e in ["AL lifelong", "DR lifelong", "AL_DR16M"]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Diet"] == e]

        print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))
        
def twoway(d, method):
    dfs_small = d[d["Diet"].isin(["AL lifelong", "DR lifelong"])]

    model = ols('{} ~ C(Age) + C(Diet) + C(Age):C(Diet)'.format(method), data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))
# }}}


KW_bwTreat(CS_df_met2, "IgM-IgD-")
KW_bwTreat(CS_df_met2, "IgM+IgD+SHM+")
KW_bwTreat(CS_df_met2, "IgM+IgD+SHM-")

KW_bwTime(CS_df_met2, "IgM-IgD-")
KW_bwTime(CS_df_met2, "IgM+IgD+SHM+")
KW_bwTime(CS_df_met2, "IgM+IgD+SHM-")

CS_df_met3 = CS_df_met2.copy()
CS_df_met3.columns = ['antibody', 'antigen', 'naive', 'Mouse', 'biogroup', 'Diet', 'Age']

CS_df_met3.head()

print("IgM-IgD-")
twoway(CS_df_met3, 'antibody')
print("IgM+IgD+SHM+")
twoway(CS_df_met3, 'antigen')
print("IgM+IgD+SHM-")
twoway(CS_df_met3, 'naive')






