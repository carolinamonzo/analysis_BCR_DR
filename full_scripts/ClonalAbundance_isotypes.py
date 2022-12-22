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
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
import matplotlib
import os
from scipy import stats
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm
from IPython.display import display
import scikit_posthocs as sp

new_day = datetime.datetime.now().strftime("%Y%m%d")
path_ori = "../../analysis/plots/"
path = "../../analysis/plots/clonal_abundance/"


#organ = "ILE"
organ = "SPL"

#run_type = "dry"
run_type = "wet"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}
palette4 = {"5AL_lifelong":"#D3E0F1", "16AL_lifelong":"#A2C2DD", "20AL_lifelong":"#6897C6", "24AL_lifelong":"#3D64A8", "5DR_lifelong":"#ECBAA1", 
            "20DR_lifelong":"#DD694D", "24DR_lifelong":"#AD1E22","20AL_DR16M":"#779D45", "24AL_DR16M":"#416F6F", "24AL_DR20M":"#EDD853"}
# }}}

# {{{
igA = pd.read_csv(path + 'clab_IGA_{}.tsv'.format(organ), sep = " ")
igM = pd.read_csv(path + 'clab_IGM_{}.tsv'.format(organ), sep = " ")
igD = pd.read_csv(path + 'clab_IGD_{}.tsv'.format(organ), sep = " ")
igE = pd.read_csv(path + 'clab_IGE_{}.tsv'.format(organ), sep = " ")
igG12 = pd.read_csv(path + 'clab_IGG12_{}.tsv'.format(organ), sep = " ")
igG3 = pd.read_csv(path + 'clab_IGG3_{}.tsv'.format(organ), sep = " ")

igG = pd.concat([igG12, igG3])
# }}}

# {{{
# Read metadata so we can get age and diet from there
nam = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8", "Illumina", "Illumina2", "Barcode"])
nam.columns = ["Mouse", "Diet", "Age", "sample_id"]

nam['Diet'] = nam['Diet'].str.replace(' ','_')

nam["biogroup"] = nam["Age"].astype(str) + nam["Diet"].astype(str)
# }}}

# {{{
migA = nam.merge(igA, on = "sample_id")
migM = nam.merge(igM, on = "sample_id")
migD = nam.merge(igD, on = "sample_id")
migE = nam.merge(igE, on = "sample_id")
migG = nam.merge(igG, on = "sample_id")

tup_dfs = [(migA, "IGA"), (migM, "IGM"), (migD, "IGD"), (migE, "IGE"), (migG, "IGG")]
# }}}

migA["Diet"].unique()

for d, iso in tup_dfs:

    fig, ax = plt.subplots(figsize = (13,6))
    sns.lineplot(data = d, x = d["rank"].astype(float), y = d["p"].astype(float), hue = migA["biogroup"], ax = ax, linewidth = 2.5, palette = palette4, hue_order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
              "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"])
    plt.xscale("log")
    plt.yscale("log")

    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Clonal rank in repertoire", fontsize = 20)
    ax.set_ylabel("Relative clonal abundance [Log]", fontsize = 20)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(2.5)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig(f"{path}/{organ}_{iso}_LogClonalAbundance_{new_day}.pdf")
    else:
        plt.show()

for d, iso in tup_dfs:

    fig, ax = plt.subplots(figsize = (13,6))
    sns.lineplot(data = d, x = d["rank"].astype(float), y = d["p"].astype(float), hue = migA["biogroup"], ax = ax, linewidth = 2.5, palette = palette4, hue_order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
              "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"])

    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Clonal rank in repertoire", fontsize = 20)
    ax.set_ylabel("Relative clonal abundance", fontsize = 20)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(2.5)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig(f"{path}/{organ}_{iso}_ClonalAbundance_{new_day}.pdf")
    else:
        plt.show()



## Clonal expansion is the sum of frequencies from all clones with rank below or equal to 20
for d, iso in tup_dfs:
    d["rank"] = d["rank"].astype(float)
    d["p"] = d["p"].astype(float)
    p20_df = d[d["rank"] <= 20]
    p20 = pd.DataFrame(p20_df.loc[:, ["biogroup", "sample_id", "p", "Age", "Diet"]].groupby(['biogroup', 'sample_id', "Age", "Diet"])["p"].sum())
    p20.reset_index(inplace = True)
    
    # Plot it all in one line
    fig, ax = plt.subplots(figsize = (9,7))
    ax = sns.boxplot(data = p20, x = "biogroup", y = "p", palette = palette4, order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
              "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"], showfliers = False, ax = ax)
    ax = sns.swarmplot(data = p20, x = "biogroup", y = "p", dodge = True, color=".25", ax = ax, order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
              "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"])
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    plt.xticks(rotation = 80)

    if organ == "ILE":
        plt.title(f"P20 Ileum [{iso}]\n", fontsize = 20)
    elif organ == "SPL":
        plt.title(f"P20 Spleen [{iso}]\n", fontsize = 20)

    ax.set_xlabel("", fontsize = 20)
    ax.set_ylabel("P20 frequency", fontsize = 20)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig(f"{path}/{organ}_{iso}_P20_biogroups_{new_day}.pdf")
    else:
        plt.show()
        
    # Plot it by hue
    fig, ax = plt.subplots(figsize = (13,5))
    ax = sns.boxplot(data = p20, x = "Age", y = "p", hue = "Diet", palette = palette2, showfliers = False, ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
    ax = sns.swarmplot(data = p20, x = "Age", y = "p", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("P20 frequency", fontsize = 20)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig("{}{}_{}_P20_AGE_{}.pdf".format(path, organ, iso, new_day))
    else:
        plt.show

# {{{
migA["rank"] = migA["rank"].astype(float)
migA["p"] = migA["p"].astype(float)
p20_igA = migA[migA["rank"] <= 20]
p20_igA = pd.DataFrame(p20_igA.loc[:, ["biogroup", "sample_id", "p", "Diet", "Age"]].groupby(['biogroup', 'sample_id', "Diet", "Age"])["p"].sum())
p20_igA.reset_index(inplace = True)

migM["rank"] = migM["rank"].astype(float)
migM["p"] = migM["p"].astype(float)
p20_igM = migM[migM["rank"] <= 20]
p20_igM = pd.DataFrame(p20_igM.loc[:, ["biogroup", "sample_id", "p", "Diet", "Age"]].groupby(['biogroup', 'sample_id', "Diet", "Age"])["p"].sum())
p20_igM.reset_index(inplace = True)

migD["rank"] = migD["rank"].astype(float)
migD["p"] = migD["p"].astype(float)
p20_igD = migD[migD["rank"] <= 20]
p20_igD = pd.DataFrame(p20_igD.loc[:, ["biogroup", "sample_id", "p", "Diet", "Age"]].groupby(['biogroup', 'sample_id', "Diet", "Age"])["p"].sum())
p20_igD.reset_index(inplace = True)

migE["rank"] = migE["rank"].astype(float)
migE["p"] = migE["p"].astype(float)
p20_igE = migE[migE["rank"] <= 20]
p20_igE = pd.DataFrame(p20_igE.loc[:, ["biogroup", "sample_id", "p", "Diet", "Age"]].groupby(['biogroup', 'sample_id', "Diet", "Age"])["p"].sum())
p20_igE.reset_index(inplace = True)

migG["rank"] = migG["rank"].astype(float)
migG["p"] = migG["p"].astype(float)
p20_igG = migG[migG["rank"] <= 20]
p20_igG = pd.DataFrame(p20_igG.loc[:, ["biogroup", "sample_id", "p", "Diet", "Age"]].groupby(['biogroup', 'sample_id', "Diet", "Age"])["p"].sum())
p20_igG.reset_index(inplace = True)

tup_dfs2 = [(p20_igA, "IGA"), (p20_igM, "IGM"), (p20_igD, "IGD"), (p20_igE, "IGE"), (p20_igG, "IGG")]
# }}}

def boxplots_cond(df_prov, order, method):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup"].isin(order)]
    
    #Linear regression
    x = df.Age.to_list()
    y = df[method].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")


# {{{
from scipy import stats
method = "p"

for d,iso in tup_dfs2:
    print("\n\n" + iso)
    print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
    for e in [5, 20, 24]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Age"] == e]

        print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))
# }}}

for d,iso in tup_dfs2:
    print("\n\n" + iso)
    print("//K-W BETWEEN timepoints WITHIN treatments//")
    for e in ["AL_lifelong", "DR_lifelong"]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Diet"] == e]

        print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))

for d,iso in tup_dfs2:
    print("\n\n" + iso)
    print("AL")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
                 "p")
    print("\nDR")
    boxplots_cond(d, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
                 "p")
    print("\nAL_DR16M")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
                 "p")
    print("\nAL_DR20M")
    boxplots_cond(d, ["20AL_lifelong", "24AL_DR20M"], 
                 "p")

for d,iso in tup_dfs2:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

    model = ols('p ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs2:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('p ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs2:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('p ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))




