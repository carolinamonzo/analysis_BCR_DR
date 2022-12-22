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
path = "../../analysis/SHM/"

run_type = "dry"
#run_type = "wet"

organ = "ILE"
#organ = "SPL"

#shm_type = "absolute"
shm_type = "frequency"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", "AL_DR16M":"teal", "AL_DR20M":"gold"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

# Reading metadata file
metadata = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";", usecols = ["Diet", "Age", "File", "Mouse"])
metadata = metadata.replace(" ", "_", regex = True)
metadata["biogroup"] = metadata["Age"].astype(str) + metadata["Diet"]

metadata.head()

# Reading data SHM
shm = pd.read_csv("../../analysis/results_tables/SHM_{}_{}.tsv".format(shm_type, organ), sep = " ")
df = pd.read_csv("../../analysis/created_germlines/reannotated_db/merged_changeo_{}.tsv".format(organ), sep = "\t")

shm.head()

# {{{
if shm_type == "frequency":
    shm.columns = ["sequence_id", "sample_id", "biogroup", "mu_count_seq_r", "mu_count_seq_s"]

# add shm columns to df
if organ == "ILE":
    df["mu_count_seq_r"] = shm["mu_count_seq_r"].to_list()
    df["mu_count_seq_s"] = shm["mu_count_seq_s"].to_list()
    df = df.merge(metadata.loc[:, ["Diet", "Age", "biogroup"]], on = "biogroup")
elif organ == "SPL":
    if shm_type != "frequency":
        df["mu_count_seq_r"] = shm["mu_count_seq_r"].to_list()
        df["mu_count_seq_s"] = shm["mu_count_seq_s"].to_list()
        df = df.merge(metadata.loc[:, ["Diet", "Age", "biogroup"]], on = "biogroup")
    else:
        df = df.merge(shm.loc[:, ["sequence_id", "sample_id", "mu_count_seq_r", "mu_count_seq_s"]], on = ["sequence_id", "sample_id"], how = "inner")
        df = df.merge(metadata.loc[:, ["Diet", "Age", "biogroup"]], on = "biogroup")
# }}}

df = df.merge(metadata.loc[:, ["Diet", "Age", "biogroup"]], on = "biogroup")

# {{{
fig, ax = plt.subplots(figsize = (10, 5))
sns.boxplot(data = df, x = "biogroup", y = "mu_count_seq_r", palette = palette3, order = gr, showfliers = False)

ax.tick_params(axis = "x", labelsize=14)
ax.tick_params(axis = "y", labelsize=14)
ax.set_xlabel("")
if shm_type != "frequency":
    ax.set_ylabel("Mu count [Non-synonymous]", fontsize = 20)
else:
    ax.set_ylabel("Mu freq [Non-synonymous]", fontsize = 20)


ax.set_xticklabels(gr, rotation = 20)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
# }}}

dfs = df.loc[:, ["sample_id", "Diet_x", "Age_x", "mu_count_seq_r", "mu_count_seq_s"]]
dfs.columns = ["sample_id", "Diet", "Age", "mu_count_seq_r", "mu_count_seq_s"]

dfs_grouped = dfs.groupby(by = "sample_id").mean()
dfs_grouped = dfs_grouped.merge(metadata.loc[:, ["File", "Diet", "biogroup"]], left_index = True, right_on = "File")

# {{{
# This has standard error of the mean
fig, ax = plt.subplots(figsize = (9, 5))
sns.boxplot(data = dfs_grouped, x = "Age", y = "mu_count_seq_s", hue = "Diet", palette = palette2, showfliers = False, 
           hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)
sns.swarmplot(data = dfs_grouped, x = "Age", y = "mu_count_seq_s", hue = "Diet", dodge = True, color = ".25", 
           hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)

ax.tick_params(axis = "x", labelsize=14)
ax.tick_params(axis = "y", labelsize=14)
ax.set_xlabel("Age [Months]", fontsize = 18)
if shm_type != "frequency":
    ax.set_ylabel("Mu count [Synonymous]", fontsize = 18)
else:
    ax.set_ylabel("Mu freq [Synonymous]", fontsize = 18)

ax.get_legend().remove()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

plt.savefig("{}{}_SHM_synonym_{}_{}.pdf".format(path, organ, shm_type, new_day))
# }}}

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

model = ols('mu_count_seq_s ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
dfs_small = dfs_small[dfs_small["Age"] > 17]

model = ols('mu_count_seq_s ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
dfs_small = dfs_small[dfs_small["Age"] > 17]

model = ols('mu_count_seq_s ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

method = "mu_count_seq_s"
print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
for e in [5, 20, 24]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = dfs_grouped[dfs_grouped["Age"] == e]

    print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
    print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))

dfs_grouped.head()

method = "mu_count_seq_s"
print("//K-W BETWEEN timepoints WITHIN treatments//")
for e in ["AL_lifelong", "DR_lifelong"]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = dfs_grouped[dfs_grouped["Diet"] == e]

    print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
    print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))

# {{{
fig, ax = plt.subplots(figsize = (10, 5))
sns.boxplot(data = dfs_grouped, x = "Age", y = "mu_count_seq_r", hue = "Diet", palette = palette2, showfliers = False, 
           hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)
sns.swarmplot(data = dfs_grouped, x = "Age", y = "mu_count_seq_r", hue = "Diet", dodge = True, color = ".25", 
           hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)

ax.tick_params(axis = "x", labelsize=14)
ax.tick_params(axis = "y", labelsize=14)
ax.set_xlabel("Age [Months]", fontsize = 18)
if shm_type != "frequency":
    ax.set_ylabel("Mu count [Non-synonymous]", fontsize = 18)
else:
    ax.set_ylabel("Mu freq [Non-synonymous]", fontsize = 18)

ax.get_legend().remove()


ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.savefig("{}{}_SHM_NOTsynonym_{}_{}.pdf".format(path, organ, shm_type, new_day))
# }}}

method = "mu_count_seq_r"
print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
for e in [5, 20, 24]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = dfs_grouped[dfs_grouped["Age"] == e]

    print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
    print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))

method = "mu_count_seq_r"
print("//K-W BETWEEN timepoints WITHIN treatments//")
for e in ["AL_lifelong", "DR_lifelong"]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = dfs_grouped[dfs_grouped["Diet"] == e]

    print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
    print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

def boxplots_cond(df_prov, order, method):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup"].isin(order)]
    
    #Linear regression
    x = df.Age.to_list()
    y = df[method].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")


print("NON-synonymous: AL")
boxplots_cond(dfs_grouped, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
             "mu_count_seq_r")
print("\nNON-synonymous: DR")
boxplots_cond(dfs_grouped, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
             "mu_count_seq_r")
print("\nNON-synonymous: AL_DR16M")
boxplots_cond(dfs_grouped, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
             "mu_count_seq_r")
print("\nNON-synonymous: AL_DR20M")
boxplots_cond(dfs_grouped, ["20AL_lifelong", "24AL_DR20M"], 
             "mu_count_seq_r")

print("Synonymous: AL")
boxplots_cond(dfs_grouped, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
             "mu_count_seq_s")
print("\nSynonymous: DR")
boxplots_cond(dfs_grouped, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
             "mu_count_seq_s")
print("\nSynonymous: AL_DR16M")
boxplots_cond(dfs_grouped, ["16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
             "mu_count_seq_s")
print("\nSynonymous: AL_DR20M")
boxplots_cond(dfs_grouped, ["20AL_lifelong", "24AL_DR20M"], 
             "mu_count_seq_s")

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
dfs_small = dfs_small[dfs_small["Age"] > 17]

model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

# {{{
dfs_small = dfs_grouped[dfs_grouped["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
dfs_small = dfs_small[dfs_small["Age"] > 17]

model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
sm.stats.anova_lm(model, typ=2)
# }}}

if run_type != "dry":
    dfs_grouped.to_csv(f"{path}{organ}_{shm_type}_values_{new_day}.csv", sep = ";", index = False)

# ## Checking isotypes

## Checking isotypes
dfi = df.loc[:, ["sample_id", "mu_count_seq_r", "mu_count_seq_s", "isotype"]]
dfi_A = dfi[dfi["isotype"] == "IGA"]
dfi_M = dfi[dfi["isotype"] == "IGM"]
dfi_G = dfi[dfi["isotype"].isin(["IGG12", "IGG3"])]
dfi_D = dfi[dfi["isotype"] == "IGD"]
dfi_E = dfi[dfi["isotype"] == "IGE"]


# {{{
iga = dfi_A.drop(columns = "isotype")
iga.columns = ["sample_id", "mu_count_seq_r_IGA", "mu_count_seq_s_IGA"]
iga["mu_freq"] = iga["mu_count_seq_r_IGA"] + iga["mu_count_seq_s_IGA"]
iga = iga.groupby(by = "sample_id").mean()
iga["mut_status_IGA"] = np.where(iga["mu_freq"] >= 0.01, 1, 0)
iga = iga.drop(columns = ["mu_freq"])
iga = iga.drop_duplicates()

igm = dfi_M.drop(columns = "isotype")
igm.columns = ["sample_id", "mu_count_seq_r_IGM", "mu_count_seq_s_IGM"]
igm["mu_freq"] = igm["mu_count_seq_r_IGM"] + igm["mu_count_seq_s_IGM"]
igm = igm.groupby(by = "sample_id").mean()
igm["mut_status_IGM"] = np.where(igm["mu_freq"] >= 0.01, 1, 0)
igm = igm.drop(columns = ["mu_freq"])
igm = igm.drop_duplicates()

igg = dfi_G.drop(columns = "isotype")
igg.columns = ["sample_id", "mu_count_seq_r_IGG", "mu_count_seq_s_IGG"]
igg["mu_freq"] = igg["mu_count_seq_r_IGG"] + igg["mu_count_seq_s_IGG"]
igg = igg.groupby(by = "sample_id").mean()
igg["mut_status_IGG"] = np.where(igg["mu_freq"] >= 0.01, 1, 0)
igg = igg.drop(columns = ["mu_freq"])
igg = igg.drop_duplicates()

igd = dfi_D.drop(columns = "isotype")
igd.columns = ["sample_id", "mu_count_seq_r_IGD", "mu_count_seq_s_IGD"]
igd["mu_freq"] = igd["mu_count_seq_r_IGD"] + igd["mu_count_seq_s_IGD"]
igd = igd.groupby(by = "sample_id").mean()
igd["mut_status_IGD"] = np.where(igd["mu_freq"] >= 0.01, 1, 0)
igd = igd.drop(columns = ["mu_freq"])
igd = igd.drop_duplicates()

ige = dfi_E.drop(columns = "isotype")
ige.columns = ["sample_id", "mu_count_seq_r_IGE", "mu_count_seq_s_IGE"]
ige = ige.groupby(by = "sample_id").mean()
ige["mu_freq"] = ige["mu_count_seq_r_IGE"] + ige["mu_count_seq_s_IGE"]
ige["mut_status_IGE"] = np.where(ige["mu_freq"] >= 0.01, 1, 0)
ige = ige.drop(columns = ["mu_freq"])
ige = ige.drop_duplicates()
# }}}

# {{{
igdf = pd.merge(iga, igm, on = "sample_id", how = "inner")
igdf = pd.merge(igdf, igg, on = "sample_id", how = "inner")
igdf = pd.merge(igdf, igd, on = "sample_id", how = "inner")
igdf = pd.merge(igdf, ige, on = "sample_id", how = "inner")
igdf = igdf.drop_duplicates()

igdf.to_csv("../../analysis/results_tables/SHM_{}_isotypes_sample_{}.csv".format(organ, new_day))
# }}}

# {{{
dfi_A = dfi_A.groupby(by = "sample_id").mean()
dfi_A = dfi_A.merge(metadata.loc[:, ["File", "Diet", "biogroup", "Age"]], left_index = True, right_on = "File")

dfi_M = dfi_M.groupby(by = "sample_id").mean()
dfi_M = dfi_M.merge(metadata.loc[:, ["File", "Diet", "biogroup", "Age"]], left_index = True, right_on = "File")

dfi_G = dfi_G.groupby(by = "sample_id").mean()
dfi_G = dfi_G.merge(metadata.loc[:, ["File", "Diet", "biogroup", "Age"]], left_index = True, right_on = "File")

dfi_D = dfi_D.groupby(by = "sample_id").mean()
dfi_D = dfi_D.merge(metadata.loc[:, ["File", "Diet", "biogroup", "Age"]], left_index = True, right_on = "File")

dfi_E = dfi_E.groupby(by = "sample_id").mean()
dfi_E = dfi_E.merge(metadata.loc[:, ["File", "Diet", "biogroup", "Age"]], left_index = True, right_on = "File")
# }}}

l_dfi = [(dfi_A, "IGA"), (dfi_M, "IGM"), (dfi_G, "IGG"), (dfi_D, "IGD"), (dfi_E, "IGE")]

for d,e in l_dfi:
    fig, ax = plt.subplots(figsize = (10, 5))
    sns.boxplot(data = d, x = "Age", y = "mu_count_seq_r", hue = "Diet", palette = palette2, showfliers = False, 
               hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)
    sns.swarmplot(data = d, x = "Age", y = "mu_count_seq_r", hue = "Diet", dodge = True, color = ".25", 
               hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)

    ax.tick_params(axis = "x", labelsize=14)
    ax.tick_params(axis = "y", labelsize=14)
    ax.set_xlabel("Age [Months]", fontsize = 18)
    if shm_type != "frequency":
        ax.set_ylabel("Mu count [Non-synonymous]", fontsize = 18)
    else:
        ax.set_ylabel("Mu freq [Non-synonymous]", fontsize = 18)

    ax.get_legend().remove()


    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    #plt.show()
    plt.savefig("{}{}_SHM_{}_NOTsynonym_{}_{}.pdf".format(path, organ, e, shm_type, new_day))

# {{{
method = "mu_count_seq_r"

for d,iso in l_dfi:
    print("\n\n" + iso)
    print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
    for e in [5, 20, 24]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Age"] == e]

        print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))
# }}}

method = "mu_count_seq_r"
for d,iso in l_dfi:
    print("\n\n" + iso)
    print("//K-W BETWEEN timepoints WITHIN treatments//")
    for e in ["AL_lifelong", "DR_lifelong"]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Diet"] == e]

        print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))

for d,iso in l_dfi:
    print("\n\n" + iso)
    print("NON-synonymous: AL")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
                 "mu_count_seq_r")
    print("\nNON-synonymous: DR")
    boxplots_cond(d, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
                 "mu_count_seq_r")
    print("\nNON-synonymous: AL_DR16M")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
                 "mu_count_seq_r")
    print("\nNON-synonymous: AL_DR20M")
    boxplots_cond(d, ["20AL_lifelong", "24AL_DR20M"], 
                 "mu_count_seq_r")

for d,iso in l_dfi:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

    model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in l_dfi:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in l_dfi:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('mu_count_seq_r ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

# ##### Synonymous

for d,e in l_dfi:
    #print(e)
    fig, ax = plt.subplots(figsize = (10, 5))
    sns.boxplot(data = d, x = "Age", y = "mu_count_seq_s", hue = "Diet", palette = palette2, showfliers = False, 
               hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)
    sns.swarmplot(data = d, x = "Age", y = "mu_count_seq_s", hue = "Diet", dodge = True, color = ".25", 
               hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], ax = ax)

    ax.tick_params(axis = "x", labelsize=14)
    ax.tick_params(axis = "y", labelsize=14)
    ax.set_xlabel("Age [Months]", fontsize = 18)
    if shm_type != "frequency":
        ax.set_ylabel("Mu count [Synonymous]", fontsize = 18)
    else:
        ax.set_ylabel("Mu freq [Synonymous]", fontsize = 18)

    ax.get_legend().remove()
    


    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    #plt.show()
    plt.savefig("{}{}_SHM_{}_synonym_{}_{}.pdf".format(path, organ, e, shm_type, new_day))

# {{{
method = "mu_count_seq_s"

for d,iso in l_dfi:
    print("\n\n" + iso)
    print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
    for e in [5, 20, 24]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Age"] == e]

        print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))
# }}}

method = "mu_count_seq_s"
for d,iso in l_dfi:
    print("\n\n" + iso)
    print("//K-W BETWEEN timepoints WITHIN treatments//")
    for e in ["AL_lifelong", "DR_lifelong"]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Diet"] == e]

        print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))

for d,iso in l_dfi:
    print("\n\n" + iso)
    print("synonymous: AL")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
                 "mu_count_seq_s")
    print("\nsynonymous: DR")
    boxplots_cond(d, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
                 "mu_count_seq_s")
    print("\nsynonymous: AL_DR16M")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
                 "mu_count_seq_s")
    print("\nsynonymous: AL_DR20M")
    boxplots_cond(d, ["20AL_lifelong", "24AL_DR20M"], 
                 "mu_count_seq_s")

for d,iso in l_dfi:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

    model = ols('mu_count_seq_s ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in l_dfi:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('mu_count_seq_s ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in l_dfi:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('mu_count_seq_s ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

if run_type != "dry":
    for dfi, iso in l_dfi:
        dfi.to_csv(f"{path}{organ}_{iso}_{shm_type}_values_{new_day}.csv", sep = ";", index = False)

# {{{
## Mutational status plots
# }}}






