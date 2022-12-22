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
#     display_name: MPI_cmc
#     language: python
#     name: mpi_cmc
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

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}
palette4 = {"5AL_lifelong":"#D3E0F1", "16AL_lifelong":"#A2C2DD", "20AL_lifelong":"#6897C6", "24AL_lifelong":"#3D64A8", "5DR_lifelong":"#ECBAA1", 
            "20DR_lifelong":"#DD694D", "24DR_lifelong":"#AD1E22","20AL_DR16M":"#779D45", "24AL_DR16M":"#416F6F", "24AL_DR20M":"#EDD853"}
# }}}

spl = pd.read_csv(path_ori + "RDI_biogroups_SPL.tsv", sep = " ")
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
spl = spl[spl["value"] != 0.0]
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
spl = spl.drop_duplicates()
spl = spl[~spl['col'].isin(['value'])]

sns.boxplot(data = spl, x = "col", y = "value", palette = palette4, showfliers = False)



def boxplots_cond(df_prov, order, str_vals, text):
    
    # Subset the info we need
    df = df_prov[df_prov["col"].isin(order)]
    
    fig, ax = plt.subplots(figsize = (6,5))
    ax = sns.boxplot(data = df, x = "col", y = "value", palette = palette4, order = order, showfliers = False, ax = ax)
    ax = sns.swarmplot(data = df, x = "col", y = "value", order = order, dodge = True, color=".25", ax = ax)
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("RDI", fontsize = 20)

    ax.set_xticklabels(str_vals)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylim(0, 12)
    
    nam = order[-1][2:]
    plt.title(nam, fontsize = 20)
    
    #Linear regression
    x = [float(w.split("_")[0][:-2]) for w in df["col"].to_list()]
    y = df.value.to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    
    if run_type == "dry":
        plt.savefig("{}{}_RDI_{}_Box_{}.pdf".format(path,text, nam, new_day))
    else:
        plt.show()


order = ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"]
str_vals = [5, 20, 24]
boxplots_cond(spl, order, str_vals, "SPL")

# {{{
order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(spl, order, str_vals, "SPL")

order = ["5DR_lifelong", "16AL_lifelong", "20DR_lifelong", "24DR_lifelong"]
boxplots_cond(spl, order, str_vals, "SPL")

order = ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"]
boxplots_cond(spl, order, str_vals, "SPL")

order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"]
boxplots_cond(spl, order, str_vals, "SPL")
# }}}

print(stats.kruskal(spl["value"].to_list(), spl["col"].to_list()))
sp.posthoc_mannwhitney(spl, val_col = "value", group_col = "col")

# {{{
## For kruskall add metadata of age and diet
# }}}
# {{{
# Getting the info we need into our dataframe
diets = []
for e in spl["row"].to_list():
    if e[2:] == "R_lifelong":
        diets.append("DR_lifelong")
    elif e[2:] == "L_lifelong":
        diets.append("AL_lifelong")
    else:
        diets.append(e[2:])

ages = []
for e in spl["row"].to_list():
    ages.append(e[:2])
ages = [w.replace('5D', "5") for w in ages]
ages = [w.replace('5A', "5") for w in ages]

ages = [ int(x) for x in ages ]

# Add to dataframe
spl["Age"] = ages
spl["Diet"] = diets
# }}}

# {{{
fig, ax = plt.subplots(figsize = (13,5))
ax = sns.boxplot(data = spl, x = "Age", y = "", hue = "Diet", palette = palette2, showfliers = False, ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
ax = sns.swarmplot(data = spl, x = "Age", y = "value", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("RDI", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type == "dry":
    plt.savefig("{}SPL_RDI_Box_all_AGE_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

# {{{
splsub = spl[spl["Age"] != 16]
palette3 = {"DR_lifelong":"#CA4B1F", "AL_lifelong":"#1A5DA4", "AL_DR16M":"#2D915F", "AL_DR20M":"#EAE437"}


fig, ax = plt.subplots(figsize = (10,6))
sns.lineplot(data = splsub, x = "Age", y = "value", hue = "Diet", ax=ax, palette=palette3,
    linewidth = 3.5, hue_order = ["DR_lifelong", "AL_DR16M", "AL_DR20M", "AL_lifelong"], marker = "o")

plt.axvline(x = 16, alpha = 0.3, color = "teal")
plt.axvline(x = 20, alpha = 0.3, color = "gold")
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Repertoire dissimilarity index", fontsize = 20)
ax.set_ylabel(dic_test[e], fontsize = 20)
plt.ylim(0, None)
plt.xlim(0, 25)

leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for line in leg.get_lines():
    line.set_linewidth(3.5)

matplotlib.rcParams['pdf.fonttype'] = 42 

plt.tight_layout()

if run_type != "dry":

    plt.savefig("{}{}_RDI_line_{}.pdf".format(path, organ, new_day))

else:
    plt.show()
    
# }}}

# {{{
sh_24 = spl[spl["Age"] == 24]

print(stats.kruskal(sh_24["value"].to_list(), sh_24["Diet"].to_list()))

sp.posthoc_conover(sh_24, val_col = "value", group_col = "Diet")
# }}}

# {{{
sh_20 = spl[spl["Age"] == 20]

print(stats.kruskal(sh_20["value"].to_list(), sh_20["Diet"].to_list()))

sp.posthoc_conover(sh_20, val_col = "value", group_col = "Diet", p_adjust=None)
# }}}

# {{{
sh_5 = spl[spl["Age"] == 5]

print(stats.kruskal(sh_5["value"].to_list(), sh_5["Diet"].to_list()))

sp.posthoc_conover(sh_5, val_col = "value", group_col = "Diet", p_adjust="fdr_bh")
# }}}

anova_lm(ols("value ~ Age * Diet + Age + Diet", data = spl[spl["Diet"].isin(["DR_lifelong", "AL_DR20M"])]).fit(), typ=2)

# {{{
sw16 = spl[spl["row"].isin(["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"])]
sw16["Diet"] = "AL_DR16M"
AL_sw16 = pd.concat([sw16, spl[spl["Diet"] == "AL_lifelong"]])

anova_lm(ols("value ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2)
# }}}

# {{{
sw20 = spl[spl["row"].isin(["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"])]
sw20["Diet"] = "AL_DR20M"
AL_sw20 = pd.concat([sw20, spl[spl["Diet"] == "AL_lifelong"]])

anova_lm(ols("value ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2)
# }}}

DR_sw16 = pd.concat([sw16, spl[spl["Diet"] == "DR_lifelong"]])
anova_lm(ols("value ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2)

DR_sw20 = pd.concat([sw20, spl[spl["Diet"] == "DR_lifelong"]])
anova_lm(ols("value ~ Age * Diet + Age + Diet", data = DR_sw20).fit(), typ=2)

sw16_sw20 = pd.concat([sw16, sw20])
anova_lm(ols("value ~ Age * Diet + Age + Diet", data = sw16_sw20).fit(), typ=2)

# {{{
## Kruskall wallis on the significance of the estimated age effect on alpha diversity
print("AL lifelong stats: \n")
test = spl[spl["Diet"] == "AL_lifelong"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))

display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))

print("DR lifelong stats: \n")
test = spl[spl["Diet"] == "DR_lifelong"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))

print("AL_DR16M stats: \n")
test = sw16[sw16["Diet"] == "AL_DR16M"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))

print("AL_DR20M stats: \n")
test = sw20[sw20["Diet"] == "AL_DR20M"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))
# }}}


# ## Same but with ileum

spl = pd.DataFrame()

spl = pd.read_csv(path_ori + "RDI_biogroups_ILE.tsv", sep = " ")
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
spl = spl[spl["value"] != 0.0]
spl = spl.drop_duplicates()
spl = spl[~spl['col'].isin(['value'])]

spl.shape

order = ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"]
str_vals = [5, 20, 24]
boxplots_cond(spl, order, str_vals, "ILE")

# {{{
order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(spl, order, str_vals, "ILE")

order = ["5DR_lifelong", "16AL_lifelong", "20DR_lifelong", "24DR_lifelong"]
boxplots_cond(spl, order, str_vals, "ILE")

order = ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"]
boxplots_cond(spl, order, str_vals, "ILE")

order = ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"]
boxplots_cond(spl, order, str_vals, "ILE")
# }}}

print(stats.kruskal(spl["value"].to_list(), spl["col"].to_list()))
sp.posthoc_mannwhitney(spl, val_col = "value", group_col = "col")

# {{{
# Getting the info we need into our dataframe
diets = []
for e in spl["col"].to_list():
    if e[2:] == "R_lifelong":
        diets.append("DR_lifelong")
    elif e[2:] == "L_lifelong":
        diets.append("AL_lifelong")
    else:
        diets.append(e[2:])

ages = []
for e in spl["col"].to_list():
    ages.append(e[:2])
ages = [w.replace('5D', "5") for w in ages]
ages = [w.replace('5A', "5") for w in ages]

ages = [ int(x) for x in ages ]

# Add to dataframe
spl["Age"] = ages
spl["Diet"] = diets
# }}}

# {{{
fig, ax = plt.subplots(figsize = (13,5))
ax = sns.boxplot(data = spl, x = "Age", y = "value", hue = "Diet", palette = palette2, showfliers = False, ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
ax = sns.swarmplot(data = spl, x = "Age", y = "value", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"])
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("RDI", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type == "dry":
    plt.savefig("{}ILE_RDI_Box_all_AGE_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

# {{{
sh_24 = spl[spl["Age"] == 24]

print(stats.kruskal(sh_24["value"].to_list(), sh_24["Diet"].to_list()))

sp.posthoc_conover(sh_24, val_col = "value", group_col = "Diet")
# }}}

# {{{
sh_20 = spl[spl["Age"] == 20]

print(stats.kruskal(sh_20["value"].to_list(), sh_20["Diet"].to_list()))

sp.posthoc_conover(sh_20, val_col = "value", group_col = "Diet", p_adjust=None)
# }}}

# {{{
sh_5 = spl[spl["Age"] == 5]

print(stats.kruskal(sh_5["value"].to_list(), sh_5["Diet"].to_list()))

sp.posthoc_conover(sh_5, val_col = "value", group_col = "Diet", p_adjust="fdr_bh")
# }}}

anova_lm(ols("value ~ Age * Diet + Age + Diet", data = spl[spl["Diet"].isin(["AL_lifelong", "AL_DR20M"])]).fit(), typ=2)

# {{{
sw16 = spl[spl["row"].isin(["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR16M"])]
sw16["Diet"] = "AL_DR16M"
AL_sw16 = pd.concat([sw16, spl[spl["Diet"] == "AL_lifelong"]])

anova_lm(ols("value ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2)
# }}}

# {{{
sw20 = spl[spl["row"].isin(["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_DR20M"])]
sw20["Diet"] = "AL_DR20M"
AL_sw20 = pd.concat([sw20, spl[spl["Diet"] == "AL_lifelong"]])

anova_lm(ols("value ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2)
# }}}

DR_sw16 = pd.concat([sw16, spl[spl["Diet"] == "DR_lifelong"]])
anova_lm(ols("value ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2)

DR_sw20 = pd.concat([sw20, spl[spl["Diet"] == "DR_lifelong"]])
anova_lm(ols("value ~ Age * Diet + Age + Diet", data = DR_sw20).fit(), typ=2)

sw16_sw20 = pd.concat([sw16, sw20])
anova_lm(ols("value ~ Age * Diet + Age + Diet", data = sw16_sw20).fit(), typ=2)

# {{{
## Kruskall wallis on the significance of the estimated age effect on alpha diversity
print("AL lifelong stats: \n")
test = spl[spl["Diet"] == "AL_lifelong"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))

display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))

print("DR lifelong stats: \n")
test = spl[spl["Diet"] == "DR_lifelong"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))

print("AL_DR16M stats: \n")
test = sw16[sw16["Diet"] == "AL_DR16M"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))

print("AL_DR20M stats: \n")
test = sw20[sw20["Diet"] == "AL_DR20M"]

print(stats.kruskal(test["value"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "value", group_col = "Age"))
# }}}




