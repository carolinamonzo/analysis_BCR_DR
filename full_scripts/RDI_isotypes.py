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
import re

new_day = datetime.datetime.now().strftime("%Y%m%d")
path_ori = "../../analysis/results_tables/"
path = "../../analysis/plots/RDI/"

run_type = "dry"
#run_type = "wet"

organ = "ILE"
#organ = "SPL"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"teal", "AL_DR20M":"gold", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}
palette4 = {"5AL_lifelong":"#D3E0F1", "16AL_lifelong":"#A2C2DD", "20AL_lifelong":"#6897C6", "24AL_lifelong":"#3D64A8", "5DR_lifelong":"#ECBAA1", 
            "20DR_lifelong":"#DD694D", "24DR_lifelong":"#AD1E22","20AL_DR16M":"#779D45", "24AL_DR16M":"#416F6F", "24AL_DR20M":"#EDD853"}
# }}}

igA = pd.read_csv(path_ori + "RDI_IGA_samp_{}_isosum.tsv".format(organ), sep = " ")
igM = pd.read_csv(path_ori + "RDI_IGM_samp_{}_isosum.tsv".format(organ), sep = " ")
igD = pd.read_csv(path_ori + "RDI_IGD_samp_{}_isosum.tsv".format(organ), sep = " ")
igE = pd.read_csv(path_ori + "RDI_IGE_samp_{}_isosum.tsv".format(organ), sep = " ")
igG = pd.read_csv(path_ori + "RDI_IGG_samp_{}_isosum.tsv".format(organ), sep = " ")

# {{{
# Drop duplicates
igA = igA[igA['row'] != igA['col']]
igM = igM[igM['row'] != igM['col']]
igD = igD[igD['row'] != igD['col']]
igE = igE[igE['row'] != igE['col']]
igG = igG[igG['row'] != igG['col']]

# Read metadata so we can get age and diet from there
nam = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8", "Illumina", "Illumina2", "Barcode"])
nam.columns = ["Mouse", "Diet", "Age", "sample_id"]

nam['Diet'] = nam['Diet'].str.replace(' ','_')

nam["biogroup"] = nam["Age"].astype(str) + nam["Diet"].astype(str)
# }}}

def format_with_metadata(ig, nam):
# Bring biogroups of the two in the comparison

    mig = pd.merge(nam.loc[:, ["sample_id", "biogroup"]], ig, left_on="sample_id", right_on="row")
    mig.columns = ["sample_id", "biogroup_row", "row", "col", "value"]
    mig = mig.drop(columns = ["sample_id"])
    mig = pd.merge(nam.loc[:, ["sample_id", "biogroup", "Diet", "Age"]], mig, left_on="sample_id", right_on="col")
    mig.columns = ['sample_id', 'biogroup_col', 'Diet', 'Age', 'biogroup_row', 'row', 'col', 'value']
    mig = mig.drop(columns = ["sample_id"])
    # Keep only what is between samples of the same biogroup
    mig = mig[mig["biogroup_col"] == mig["biogroup_row"]]
    
    return(mig)


# {{{
migA = format_with_metadata(igA, nam)
migM = format_with_metadata(igM, nam)
migD = format_with_metadata(igD, nam)
migE = format_with_metadata(igE, nam)
migG = format_with_metadata(igG, nam)

tup_dfs = [(migA, "IGA"), (migM, "IGM"), (migD, "IGD"), (migE, "IGE"), (migG, "IGG")]
# }}}

def boxplots_cond(df_prov, order, method):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup_col"].isin(order)]
    
    #Linear regression
    x = df.Age.to_list()
    y = df[method].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")


for d,iso in tup_dfs:
    print("\n\n" + iso)
    print("AL")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
                 "value")
    print("\nDR")
    boxplots_cond(d, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
                 "value")
    print("\nAL_DR16M")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
                 "value")
    print("\nAL_DR20M")
    boxplots_cond(d, ["20AL_lifelong", "24AL_DR20M"], 
                 "value")

for d,iso in tup_dfs:
    fig, ax = plt.subplots(figsize = (13,5))
    ax = sns.boxplot(data = d, x = "Age", y = "value", hue = "Diet", palette = palette2, showfliers = False, ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
    ax = sns.swarmplot(data = d, x = "Age", y = "value", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("RDI [{}]".format(iso), fontsize = 20)

    #str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
    #ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig("{}SPL_RDI_{}_Box_{}.pdf".format(path, iso, new_day))
    else:

        plt.show()

# {{{
from scipy import stats
method = "value"

for d,iso in tup_dfs:
    print("\n\n" + iso)
    print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
    for e in [5, 20, 24]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Age"] == e]

        print(stats.kruskal(test[method].to_list(), test["Diet"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Diet", p_adjust = 'fdr_bh'))
# }}}

for d,iso in tup_dfs:
    print("\n\n" + iso)
    print("//K-W BETWEEN timepoints WITHIN treatments//")
    for e in ["AL_lifelong", "DR_lifelong"]:
        print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
        test = d[d["Diet"] == e]

        print(stats.kruskal(test[method].to_list(), test["Age"].to_list()))
        print("Mann-Whitney comparing treatments at {} Months {}\n".format(e, method))
        print(sp.posthoc_mannwhitney(test, val_col = method, group_col = "Age", p_adjust = 'fdr_bh'))

for d,iso in tup_dfs:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

    model = ols('value ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('value ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('value ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))


