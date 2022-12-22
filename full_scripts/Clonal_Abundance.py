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

# ## CLONAL ABUNDANCE, values calculed on rstudio

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


tissue = "ILE"
#tissue = "SPL"

run_type = "dry"
#run_type = "wet"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}
palette4 = {"5AL_lifelong":"#D3E0F1", "16AL_lifelong":"#A2C2DD", "20AL_lifelong":"#6897C6", "24AL_lifelong":"#3D64A8", "5DR_lifelong":"#ECBAA1", 
            "20DR_lifelong":"#DD694D", "24DR_lifelong":"#AD1E22","20AL_DR16M":"#779D45", "24AL_DR16M":"#416F6F", "24AL_DR20M":"#EDD853"}
# }}}

if tissue == "ILE":
    df = pd.read_csv(path_ori + 'clonal_abundance_values_ILE.tsv', sep = " ")
elif tissue == "SPL":
    df = pd.read_csv(path_ori + "clonal_abundance_values_SPL.tsv", sep = " ")

if tissue == "ILE":
    files = []
    for i in os.listdir(path):
        if os.path.isfile(os.path.join(path,i)) and i.startswith("clab") and i.endswith("ILE.tsv"):
            files.append(i)
    print(files)
elif tissue == "SPL":
    files = []
    for i in os.listdir(path):
        if os.path.isfile(os.path.join(path,i)) and i.startswith("clab") and i.endswith("SPL.tsv"):
            files.append(i)
    print(files)

# {{{
if tissue == "ILE":
    AL5 = pd.read_csv(path + 'clab_5AL_lifelong_ILE.tsv', sep = " ")
    DR5 = pd.read_csv(path + 'clab_5DR_lifelong_ILE.tsv', sep = " ")
    AL16 = pd.read_csv(path + 'clab_16AL_lifelong_ILE.tsv', sep = " ")
    AL20 = pd.read_csv(path + 'clab_20AL_lifelong_ILE.tsv', sep = " ")
    DR20 = pd.read_csv(path + 'clab_20DR_lifelong_ILE.tsv', sep = " ")
    AL24 = pd.read_csv(path + 'clab_24AL_lifelong_ILE.tsv', sep = " ")
    DR24 = pd.read_csv(path + 'clab_24DR_lifelong_ILE.tsv', sep = " ")
    ALDR16M20 = pd.read_csv(path + 'clab_20AL_DR16M_ILE.tsv', sep = " ")
    ALDR16M24 = pd.read_csv(path + 'clab_24AL_DR16M_ILE.tsv', sep = " ")
    ALDR20M24 = pd.read_csv(path + 'clab_24AL_DR20M_ILE.tsv', sep = " ")   

elif tissue == "SPL":
    AL5 = pd.read_csv(path + 'clab_5AL_lifelong_SPL.tsv', sep = " ")
    DR5 = pd.read_csv(path + 'clab_5DR_lifelong_SPL.tsv', sep = " ")
    AL16 = pd.read_csv(path + 'clab_16AL_lifelong_SPL.tsv', sep = " ")
    AL20 = pd.read_csv(path + 'clab_20AL_lifelong_SPL.tsv', sep = " ")
    DR20 = pd.read_csv(path + 'clab_20DR_lifelong_SPL.tsv', sep = " ")
    AL24 = pd.read_csv(path + 'clab_24AL_lifelong_SPL.tsv', sep = " ")
    DR24 = pd.read_csv(path + 'clab_24DR_lifelong_SPL.tsv', sep = " ")
    ALDR16M20 = pd.read_csv(path + 'clab_20AL_DR16M_SPL.tsv', sep = " ")
    ALDR16M24 = pd.read_csv(path + 'clab_24AL_DR16M_SPL.tsv', sep = " ")
    ALDR20M24 = pd.read_csv(path + 'clab_24AL_DR20M_SPL.tsv', sep = " ")

gr = {"5AL_lifelong":AL5, "16AL_lifelong":AL16, "20AL_lifelong":AL20, "24AL_lifelong":AL24, 
          "5DR_lifelong":DR5, "20DR_lifelong":DR20, "24DR_lifelong":DR24, "20AL_DR16M":ALDR16M20, "24AL_DR16M":ALDR16M24, "24AL_DR20M":ALDR20M24}
# }}}



# {{{
da = pd.concat(gr, names=["biogroup"])
da = da.reset_index(level =0)

da = da[da["p"] != "p_sd"]

da["p"] = da["p"].astype(float)
da["rank"] = da["rank"].astype(float)
# }}}

def plot_abundance_log(df, text):
    fig, ax = plt.subplots(figsize = (13,6))
    sns.lineplot(data = df, x = df["rank"], y = df["p"], hue = df["biogroup"], ax = ax, linewidth = 2.5, palette = palette4, hue_order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
          "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"])
    plt.xscale("log")
    plt.yscale("log")

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
        plt.savefig("{}2Clonal_AbundanceLOG_{}_{}.pdf".format(path, text, new_day))
    else:
        plt.show()


def plot_abundance(df, text):
    fig, ax = plt.subplots(figsize = (13,6))
    sns.lineplot(data = df, x = df["rank"], y = df["p"], hue = df["biogroup"], ax = ax, linewidth = 2.5, palette = palette4, hue_order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
          "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"])
    plt.xscale("log")

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
        plt.savefig("{}2Clonal_Abundance_{}_{}.pdf".format(path, text, new_day))
    else:
        plt.show()


# Plot Ileum clonal abundance
plot_abundance(da, tissue)
plot_abundance_log(da, tissue)

## Clonal expansion is the sum of frequencies from all clones with rank below or equal to 20
p20_df = da[da["rank"] <= 20]

p20 = pd.DataFrame(p20_df.loc[:, ["biogroup", "sample_id", "p"]].groupby(['biogroup', 'sample_id'])["p"].sum())
p20.reset_index(inplace = True)


p20.to_csv("{}/{}_P20_values.csv".format(path, tissue), index = False)

# {{{
fig, ax = plt.subplots(figsize = (9,7))
ax = sns.boxplot(data = p20, x = "biogroup", y = "p", palette = palette4, order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
          "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"], showfliers = False, ax = ax)
ax = sns.swarmplot(data = p20, x = "biogroup", y = "p", dodge = True, color=".25", ax = ax, order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong", 
          "5DR_lifelong", "20DR_lifelong", "24DR_lifelong", "20AL_DR16M", "24AL_DR16M", "24AL_DR20M"])
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
plt.xticks(rotation = 80)

if tissue == "ILE":
    plt.title("P20 Ileum\n", fontsize = 20)
elif tissue == "SPL":
    plt.title("P20 Spleen\n", fontsize = 20)


ax.set_xlabel("", fontsize = 20)
ax.set_ylabel("P20 frequency", fontsize = 20)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.show()
# }}}

# {{{
fig, ax = plt.subplots(figsize = (13,5))
ax = sns.boxplot(data = p20, x = "Age", y = "p", hue = "Diet", palette = palette2, showfliers = False, ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
ax = sns.swarmplot(data = p20, x = "Age", y = "p", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("P20 frequency", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}ILE_P20_AGE_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

def boxplots_cond(df_prov, order, str_vals, text):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup"].isin(order)]
    
    fig, ax = plt.subplots(figsize = (6,5))
    ax = sns.boxplot(data = df, x = "biogroup", y = "p", palette = palette4, order = order, showfliers = False, ax = ax)
    ax = sns.swarmplot(data = df, x = "biogroup", y = "p", order = order, dodge = True, color=".25", ax = ax)
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("P20 frequency", fontsize = 20)

    ax.set_xticklabels(str_vals)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if text == "ILE":
        plt.ylim(0, 0.6)
    elif text == "SPL":
        plt.ylim(0, 0.9)
    
    nam = order[-1][2:]
    plt.title(nam, fontsize = 20)
    
    #Linear regression
    x = [float(w.split("_")[0][:-2]) for w in df["biogroup"].to_list()]
    y = df.p.to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    
    if run_type != "dry":
        plt.savefig("{}{}_P20_{}_Box_{}.pdf".format(path,text, nam, new_day))
    else:
        plt.show()


order = ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"]
str_vals = [5, 20, 24]
boxplots_cond(p20, order, str_vals, tissue)

# {{{
order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(p20, order, str_vals, tissue)

order = ["5DR_lifelong", "16AL_lifelong", "20DR_lifelong", "24DR_lifelong"]
boxplots_cond(p20, order, str_vals, tissue)

order = ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"]
boxplots_cond(p20, order, str_vals, tissue)

order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"]
boxplots_cond(p20, order, str_vals, tissue)
# }}}



# {{{
# Next do stats on the P20
# }}}
print(stats.kruskal(p20["p"].to_list(), p20["biogroup"].to_list()))
sp.posthoc_mannwhitney(p20, val_col = "p", group_col = "biogroup")


sp.posthoc_conover(p20, val_col = "p", group_col = "biogroup")

# {{{
# Getting the info we need into our dataframe
diets = []
for e in p20["biogroup"].to_list():
    if e[2:] == "R_lifelong":
        diets.append("DR_lifelong")
    elif e[2:] == "L_lifelong":
        diets.append("AL_lifelong")
    else:
        diets.append(e[2:])

ages = []
for e in p20["biogroup"].to_list():
    ages.append(e[:2])
ages = [w.replace('5D', "5") for w in ages]
ages = [w.replace('5A', "5") for w in ages]

ages = [ int(x) for x in ages ]

# Add to dataframe
p20["Age"] = ages
p20["Diet"] = diets
# }}}

# {{{
sh_24 = p20[p20["Age"] == 24]

print(stats.kruskal(sh_24["p"].to_list(), sh_24["Diet"].to_list()))

print("Mann whitney")

display(sp.posthoc_mannwhitney(sh_24, val_col = "p", group_col = "Diet"))

print("Conover")
display(sp.posthoc_conover(sh_24, val_col = "p", group_col = "Diet"))
# }}}

# {{{
sh_20 = p20[p20["Age"] == 20]

print(stats.kruskal(sh_20["p"].to_list(), sh_20["Diet"].to_list()))

print("Mann whitney")

display(sp.posthoc_mannwhitney(sh_20, val_col = "p", group_col = "Diet"))

print("Conover")
display(sp.posthoc_conover(sh_20, val_col = "p", group_col = "Diet"))
# }}}

# {{{
sh_5 = p20[p20["Age"] == 5]

print(stats.kruskal(sh_5["p"].to_list(), sh_5["Diet"].to_list()))

print("Mann whitney")

display(sp.posthoc_mannwhitney(sh_5, val_col = "p", group_col = "Diet"))

print("Conover")
display(sp.posthoc_conover(sh_5, val_col = "p", group_col = "Diet"))
# }}}

anova_lm(ols("p ~ Age * Diet + Age + Diet", data = p20[p20["Diet"].isin(["AL_DR16M", "DR_lifelong"])]).fit(), typ=2)

# {{{
sw16 = p20[p20["biogroup"].isin(["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"])]
sw16["Diet"] = "AL_DR16M"
AL_sw16 = pd.concat([sw16, p20[p20["Diet"] == "AL_lifelong"]])

anova_lm(ols("p ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2)
# }}}

# {{{
sw20 = p20[p20["biogroup"].isin(["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"])]
sw20["Diet"] = "AL_DR20M"
AL_sw20 = pd.concat([sw20, p20[p20["Diet"] == "AL_lifelong"]])

anova_lm(ols("p ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2)
# }}}

DR_sw16 = pd.concat([sw16, p20[p20["Diet"] == "DR_lifelong"]])
anova_lm(ols("p ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2)

DR_sw20 = pd.concat([sw20, p20[p20["Diet"] == "DR_lifelong"]])
anova_lm(ols("p ~ Age * Diet + Age + Diet", data = DR_sw20).fit(), typ=2)

sw16_sw20 = pd.concat([sw16, sw20])
anova_lm(ols("p ~ Age * Diet + Age + Diet", data = sw16_sw20).fit(), typ=2)

# {{{
## Kruskall wallis on the significance of the estimated age effect on alpha diversity
print("AL lifelong stats: \n")
test = p20[p20["Diet"] == "AL_lifelong"]

print(stats.kruskal(test["p"].to_list(), test["Age"].to_list()))

display(sp.posthoc_mannwhitney(test, val_col = "p", group_col = "Age"))

print("DR lifelong stats: \n")
test = p20[p20["Diet"] == "DR_lifelong"]

print(stats.kruskal(test["p"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "p", group_col = "Age"))

print("AL_DR16M stats: \n")
test = sw16[sw16["Diet"] == "AL_DR16M"]

print(stats.kruskal(test["p"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "p", group_col = "Age"))

print("AL_DR20M stats: \n")
test = sw20[sw20["Diet"] == "AL_DR20M"]

print(stats.kruskal(test["p"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "p", group_col = "Age"))
# }}}


