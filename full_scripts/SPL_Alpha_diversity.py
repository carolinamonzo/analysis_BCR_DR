# ---
# jupyter:
#   jupytext:
#     cell_markers: '{{{,}}}'
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
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

new_day = datetime.datetime.now().strftime("%Y%m%d")
path_ori = "../../analysis/plots/"
path = "../../analysis/plots/alpha_diversity/"

run_type = "dry"
#run_type = "wet"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}
# }}}

# {{{
files = []
for i in os.listdir(path_ori):
    if os.path.isfile(os.path.join(path_ori,i)) and i.startswith("alpha_values") and i.endswith("SPL.tsv"):
        files.append(i)

files = [x for x in files if "stats" not in x]
#files.remove('alpha_values_ILE_small.tsv')
files.remove('alpha_values_SPL.tsv')
print(files)
# }}}

# {{{
AL5 = pd.read_csv(path_ori + 'alpha_values_5AL_SPL.tsv', sep = " ")
DR5 = pd.read_csv(path_ori + 'alpha_values_5DR_SPL.tsv', sep = " ")
AL16 = pd.read_csv(path_ori + 'alpha_values_16AL_SPL.tsv', sep = " ")
AL20 = pd.read_csv(path_ori + 'alpha_values_20AL_SPL.tsv', sep = " ")
DR20 = pd.read_csv(path_ori + 'alpha_values_20DR_SPL.tsv', sep = " ")
AL24 = pd.read_csv(path_ori + 'alpha_values_24AL_SPL.tsv', sep = " ")
DR24 = pd.read_csv(path_ori + 'alpha_values_24DR_SPL.tsv', sep = " ")
ALDR16M20 = pd.read_csv(path_ori + 'alpha_values_20AL_DR16M_SPL.tsv', sep = " ")
ALDR16M24 = pd.read_csv(path_ori + 'alpha_values_24AL_DR16M_SPL.tsv', sep = " ")
ALDR20M24 = pd.read_csv(path_ori + 'alpha_values_24AL_DR20M_SPL.tsv', sep = " ")

gr = {"5M AL_lifelong":AL5, "16M AL_lifelong":AL16, "20M AL_lifelong":AL20, "24M AL_lifelong":AL24, 
      "5M DR_lifelong":DR5, "20M DR_lifelong":DR20, "24M DR_lifelong":DR24, "20M AL_DR16M":ALDR16M20, "24M AL_DR16M":ALDR16M24, "24M AL_DR20M":ALDR20M24}

df = pd.concat(gr, names=["biogroup"])
df = df.reset_index(level =0)

new = df["biogroup"].str.split("M ", n=1, expand=True)
df["Age"] = [float(i) for i in new[0]]
df["Diet"] = new[1]
# }}}



# {{{
fig, ax = plt.subplots(figsize = (11,6))
sns.lineplot(data = df, x = df["q"], y = df["d"], hue = df["biogroup"], palette = palette3, ax = ax, linewidth = 2.5)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Diversity order (q)", fontsize = 20)
ax.set_ylabel("Alpha-diversity (qD)", fontsize = 20)

leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for line in leg.get_lines():
    line.set_linewidth(2.5)
    
plt.axvline(x = 0, alpha = 0.5, color = "black", linestyle = "--")
plt.axvline(x = 1, alpha = 0.5, color = "black", linestyle = "--")
plt.axvline(x = 2, alpha = 0.5, color = "black", linestyle = "--")
plt.text(x = -0.12, y = 7500, s = "Richness", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
plt.text(x = 0.88, y = 7500, s = "Shannon", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
plt.text(x = 1.88, y = 7500, s = "Simpson", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)

str_vals = [0, 1, 2, 3]
ax.set_xticks([0.0, 1.0, 2.0, 3.0])
ax.set_xticklabels(str_vals)

str_yvals = [0, 5000, 10000, 15000]
ax.set_yticks([0, 5000, 10000, 15000])
ax.set_yticklabels(str_yvals)    
    
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig("{}Alpha_Hill_SPL_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

# {{{
shannon_df = df[df["q"] == 1]
richness_df = df[df["q"] == 0]

fig, ax = plt.subplots(figsize = (12,7))
order = ["5M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M", "20M AL_lifelong", 
      "24M AL_DR20M", "24M AL_lifelong"]
ax = sns.boxplot(data = shannon_df, x = "biogroup", y = "d", palette = palette3, order = order, showfliers = False, ax = ax)
ax = sns.swarmplot(data = shannon_df, x = "biogroup", y = "d", order = order, dodge = True, color=".25", ax = ax)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Shannon - Alpha diversity", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])
ax.set_xticklabels(order, rotation = 80)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}ShannonBox_all_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

# {{{
fig, ax = plt.subplots(figsize = (13,5))
ax = sns.boxplot(data = shannon_df, x = "Age", y = "d", hue = "Diet", palette = palette2, showfliers = False, ax = ax)
ax = sns.swarmplot(data = shannon_df, x = "Age", y = "d", hue = "Diet", dodge = True, color=".25", ax = ax)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Shannon - Alpha diversity", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}SPL_ShannonBox_all_AGE_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

# {{{
fig, ax = plt.subplots(figsize = (13,5))
ax = sns.boxplot(data = richness_df, x = "Age", y = "d", hue = "Diet", palette = palette2, showfliers = False, ax = ax)
ax = sns.swarmplot(data = richness_df, x = "Age", y = "d", hue = "Diet", dodge = True, color=".25", ax = ax)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Richness - Alpha diversity", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type == "dry":
    plt.savefig("{}SPL_RichnessBox_all_AGE_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

def boxplots_cond(df_prov, order, str_vals, text):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup"].isin(order)]
    
    fig, ax = plt.subplots(figsize = (6,5))
    ax = sns.boxplot(data = df, x = "biogroup", y = "d", palette = palette3, order = order, showfliers = False, ax = ax)
    ax = sns.swarmplot(data = df, x = "biogroup", y = "d", order = order, dodge = True, color=".25", ax = ax)
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("{} - Alpha diversity".format(text), fontsize = 20)

    ax.set_xticklabels(str_vals)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    if text == "Shannon":
        plt.ylim(0, 1020)
    elif text == "Simpson":
        plt.ylim(0, 200)
    else:
        exit(1)
    
    nam = order[-1].split(" ")[1]
    plt.title(nam, fontsize = 20)
    
    #Linear regression
    x = [float(w.split("M")[0]) for w in df["biogroup"].to_list()]
    y = df.d.to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    
    if run_type != "dry":
        plt.savefig("{}{}Box_{}_{}.pdf".format(path, text, nam, new_day))
    else:
        plt.show()


# {{{
# REAL DR VALUES

order = ["5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
str_vals = [5, 20, 24]
boxplots_cond(richness_df, order, str_vals, "Richness")
# }}}

# {{{
order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(richness_df, order, str_vals, "Richness")

order = ["5M DR_lifelong", "16M AL_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
boxplots_cond(richness_df, order, str_vals, "Richness")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"]
boxplots_cond(richness_df, order, str_vals, "Richness")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"]
boxplots_cond(richness_df, order, str_vals, "Richness")
# }}}





order = ["5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
str_vals = [5, 20, 24]
#df = boxplots_cond(shannon_df, order, str_vals, "Shannon")



# {{{
order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(shannon_df, order, str_vals, "Shannon")

order = ["5M DR_lifelong", "16M AL_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
boxplots_cond(shannon_df, order, str_vals, "Shannon")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"]
boxplots_cond(shannon_df, order, str_vals, "Shannon")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"]
boxplots_cond(shannon_df, order, str_vals, "Shannon")
# }}}

# {{{
simpson_df = df[df["q"] == 2]

fig, ax = plt.subplots(figsize = (12,7))
order = ["5M DR_lifelong", "5M AL_lifelong", "20M DR_lifelong", "16M AL_lifelong", "24M DR_lifelong", "20M AL_DR16M","20M AL_lifelong", "24M AL_DR16M", 
      "24M AL_DR20M", "24M AL_lifelong"]
ax = sns.boxplot(data = simpson_df, x = "biogroup", y = "d", palette = palette3, order = order, showfliers = False, ax = ax)
ax = sns.swarmplot(data = simpson_df, x = "biogroup", y = "d", order = order, dodge = True, color=".25", ax = ax)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Simpson - Alpha diversity", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])
ax.set_xticklabels(order, rotation = 80)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}SimpsonBox_all_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}
# {{{
fig, ax = plt.subplots(figsize = (13,5))
ax = sns.boxplot(data = simpson_df, x = "Age", y = "d", hue = "Diet", palette = palette2, showfliers = False, ax = ax)
ax = sns.swarmplot(data = simpson_df, x = "Age", y = "d", hue = "Diet", dodge = True, color=".25", ax = ax)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Simpson - Alpha diversity", fontsize = 20)

#str_vals = [5, 16, 20, 24, 5, 20, 24, 20, 24, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type == "dry":
    plt.savefig("{}SPL_SimpsonBox_all_AGE_{}.pdf".format(path, new_day))
else:

    plt.show()
# }}}

simpson_df[simpson_df["biogroup"] == "16M AL_lifelong"]["d"].median()


order = ["20M AL_DR16M", "24M AL_DR16M"]
str_vals = [5, 20, 24]
boxplots_cond(simpson_df, order, str_vals, "Shannon")

# {{{
order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(simpson_df, order, str_vals, "Simpson")

order = ["5M DR_lifelong", "16M AL_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
boxplots_cond(simpson_df, order, str_vals, "Simpson")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"]
boxplots_cond(simpson_df, order, str_vals, "Simpson")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"]
boxplots_cond(simpson_df, order, str_vals, "Simpson")
# }}}

# ### Trying to fit the linear regression on top of our plots

test_df = df[df["q"] == 2]
order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong"]
test_df = test_df[test_df["biogroup"].isin(order)]

test_df

# {{{
fig, ax = plt.subplots(figsize = (12,5))

#Linear regression
x = [float(w.split("M")[0]) for w in test_df["biogroup"].to_list()]
y = test_df.d.to_list()
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

ax = sns.lineplot(x = x, y = intercept + [slope *i for i in x])

str_vals = [5.0, 16.0, 20.0, 24.0]
str_col = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8"]
color_dict = dict(zip(str_vals, str_col))

ax = sns.boxplot(data = test_df, x = "Age", y = "d", palette = color_dict, ax = ax)
ax = sns.swarmplot(data = test_df, x = "Age", y = "d", dodge = True, color=".25", ax = ax)
ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Simpson - Alpha diversity", fontsize = 20)

#str_vals = [5, 16, 20, 24]
#ax.set_xticks(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong", "5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong", "20M AL_DR16M", "24M AL_DR16M", "24M AL_DR20M"])
#ax.set_xticklabels(str_vals)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)



matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.plot()
# }}}

sns.lineplot(x = x, y = intercept + [slope *i for i in x])

print(stats.kruskal(shannon_df["d"].to_list(), shannon_df["biogroup"].to_list()))
sp.posthoc_conover(shannon_df, val_col = "d", group_col = "biogroup")

# {{{
## Kruskall wallis on the significance of the estimated age effect on alpha diversity
print("AL lifelong stats: \n")
test = shannon_df[shannon_df["Diet"] == "AL_lifelong"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))

display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))

print("DR lifelong stats: \n")
test = shannon_df[shannon_df["Diet"] == "DR_lifelong"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))

print("AL_DR16M stats: \n")
test = sw16[sw16["Diet"] == "AL_DR16M"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))

print("AL_DR20M stats: \n")
test = sw20[sw20["Diet"] == "AL_DR20M"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))
# }}}

# {{{
sh_24 = shannon_df[shannon_df["Age"] == 24]

print(stats.kruskal(sh_24["d"].to_list(), sh_24["Diet"].to_list()))

sp.posthoc_conover(sh_24, val_col = "d", group_col = "Diet")
# }}}
# {{{
sh_20 = shannon_df[shannon_df["Age"] == 20]

print(stats.kruskal(sh_20["d"].to_list(), sh_20["Diet"].to_list()))

sp.posthoc_conover(sh_20, val_col = "d", group_col = "Diet", p_adjust=None)
# }}}

# {{{
sh_5 = shannon_df[shannon_df["Age"] == 5]

print(stats.kruskal(sh_5["d"].to_list(), sh_5["Diet"].to_list()))

sp.posthoc_conover(sh_5, val_col = "d", group_col = "Diet", p_adjust="fdr_bh")
# }}}

anova_lm(ols("d ~ Age * Diet + Age + Diet", data = shannon_df[shannon_df["Diet"].isin(["DR_lifelong", "AL_DR20M"])]).fit(), typ=2)

# {{{
sw16 = shannon_df[shannon_df["biogroup"].isin(["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"])]
sw16["Diet"] = "AL_DR16M"
AL_sw16 = pd.concat([sw16, shannon_df[shannon_df["Diet"] == "AL_lifelong"]])

anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2)
# }}}

# {{{
sw20 = shannon_df[shannon_df["biogroup"].isin(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"])]
sw20["Diet"] = "AL_DR20M"
AL_sw20 = pd.concat([sw20, shannon_df[shannon_df["Diet"] == "AL_lifelong"]])

anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2)
# }}}

DR_sw16 = pd.concat([sw16, shannon_df[shannon_df["Diet"] == "DR_lifelong"]])
anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2)

DR_sw20 = pd.concat([sw20, shannon_df[shannon_df["Diet"] == "DR_lifelong"]])
anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw20).fit(), typ=2)

sw16_sw20 = pd.concat([sw16, sw20])
anova_lm(ols("d ~ Age * Diet + Age + Diet", data = sw16_sw20).fit(), typ=2)

# ## And now with simpson

# {{{
## Kruskall wallis on the significance of the estimated age effect on alpha diversity
print("AL lifelong stats: \n")
test = simpson_df[simpson_df["Diet"] == "AL_lifelong"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))

print("DR lifelong stats: \n")
test = simpson_df[simpson_df["Diet"] == "DR_lifelong"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))

print("AL_DR16M stats: \n")
test = sw16[sw16["Diet"] == "AL_DR16M"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))

print("AL_DR20M stats: \n")
test = sw20[sw20["Diet"] == "AL_DR20M"]

print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
display(sp.posthoc_mannwhitney(test, val_col = "d", group_col = "Age"))
# }}}

anova_lm(ols("d ~ Age * Diet + Age + Diet", data = simpson_df[simpson_df["Diet"].isin(["AL_lifelong", "AL_DR20M"])]).fit(), typ=2)

# {{{
sw16 = simpson_df[simpson_df["biogroup"].isin(["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"])]
sw16["Diet"] = "AL_DR16M"
AL_sw16 = pd.concat([sw16, simpson_df[simpson_df["Diet"] == "AL_lifelong"]])

anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2)
# }}}

# {{{
sw20 = simpson_df[simpson_df["biogroup"].isin(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"])]
sw20["Diet"] = "AL_DR20M"
AL_sw20 = pd.concat([sw20, simpson_df[simpson_df["Diet"] == "AL_lifelong"]])

anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2)
# }}}

DR_sw16 = pd.concat([sw16, simpson_df[simpson_df["Diet"] == "DR_lifelong"]])
anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2)

DR_sw20 = pd.concat([sw20, simpson_df[simpson_df["Diet"] == "DR_lifelong"]])
anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw20).fit(), typ=2)

sw16_sw20 = pd.concat([sw16, sw20])
anova_lm(ols("d ~ Age * Diet + Age + Diet", data = sw16_sw20).fit(), typ=2)

print(stats.kruskal(simpson_df["d"].to_list(), simpson_df["biogroup"].to_list()))
sp.posthoc_mannwhitney(simpson_df, val_col = "d", group_col = "biogroup")

# {{{
sh_24 = simpson_df[simpson_df["Age"] == 24]

print(stats.kruskal(sh_24["d"].to_list(), sh_24["Diet"].to_list()))

sp.posthoc_conover(sh_24, val_col = "d", group_col = "Diet")
# }}}

# {{{
sh_20 = simpson_df[simpson_df["Age"] == 20]

print(stats.kruskal(sh_20["d"].to_list(), sh_20["Diet"].to_list()))

sp.posthoc_conover(sh_20, val_col = "d", group_col = "Diet", p_adjust=None)
# }}}

# {{{
sh_5 = simpson_df[simpson_df["Age"] == 5]

print(stats.kruskal(sh_5["d"].to_list(), sh_5["Diet"].to_list()))

sp.posthoc_conover(sh_5, val_col = "d", group_col = "Diet", p_adjust="fdr_bh")
# }}}
def do_stats_isotype(df):
    
    print(stats.kruskal(df["d"].to_list(), df["biogroup"].to_list()))
    display(sp.posthoc_conover(df, val_col = "d", group_col = "biogroup"))
    
    print("Anova AL/DR\n")
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = df[df["Diet"].isin(["AL_lifelong", "DR_lifelong"])]).fit(), typ=2))
    
    print("ANOVA with AL before switch\n")
    print("Anova AL/AL_DR16M\n")
    
    sw16 = df[df["biogroup"].isin(["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"])]
    sw16["Diet"] = "AL_DR16M"
    AL_sw16 = pd.concat([sw16, df[df["Diet"] == "AL_lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2))
    
    print("Anova AL/AL_DR20M\n")
    sw20 = shannon_df[shannon_df["biogroup"].isin(["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"])]
    sw20["Diet"] = "AL_DR20M"
    AL_sw20 = pd.concat([sw20, df[df["Diet"] == "AL_lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2))
    
    print("Anova AL_DR16M/DR\n")
    DR_sw16 = pd.concat([sw16, df[df["Diet"] == "DR_lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2))
    
    print("Anova AL_DR20M/DR\n")
    DR_sw20 = pd.concat([sw20, df[df["Diet"] == "DR_lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw20).fit(), typ=2))
    
    print("Anova AL_DR16M/AL_DR20M\n")
    sw16_sw20 = pd.concat([sw16, sw20])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = sw16_sw20).fit(), typ=2))
    
    print("\nKruskal/conover at 24M\n")
    
    sh_24 = df[df["Age"] == 24]

    print(stats.kruskal(sh_24["d"].to_list(), sh_24["Diet"].to_list()))

    display(sp.posthoc_conover(sh_24, val_col = "d", group_col = "Diet"))
    
    print("Kruskal/conover at 20M\n")
    
    sh_20 = df[df["Age"] == 20]

    print(stats.kruskal(sh_20["d"].to_list(), sh_20["Diet"].to_list()))

    display(sp.posthoc_conover(sh_20, val_col = "d", group_col = "Diet"))
    
    print("Kruskal/conover at 5M\n")
    
    sh_5 = df[df["Age"] == 5]

    print(stats.kruskal(sh_5["d"].to_list(), sh_5["Diet"].to_list()))

    display(sp.posthoc_conover(sh_5, val_col = "d", group_col = "Diet"))
    
    print("\nKruskal between ages within diet\n")
    ## Kruskall wallis on the significance of the estimated age effect on alpha diversity
    print("AL lifelong stats: \n")
    test = df[df["Diet"] == "AL_lifelong"]

    print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))

    display(sp.posthoc_conover(test, val_col = "d", group_col = "Age"))

    print("DR lifelong stats: \n")
    test = df[df["Diet"] == "DR_lifelong"]

    print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
    display(sp.posthoc_conover(test, val_col = "d", group_col = "Age"))

    print("AL_DR16M stats: \n")
    test = sw16[sw16["Diet"] == "AL_DR16M"]

    print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
    display(sp.posthoc_conover(test, val_col = "d", group_col = "Age"))

    print("AL_DR20M stats: \n")
    test = sw20[sw20["Diet"] == "AL_DR20M"]

    print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))
    display(sp.posthoc_conover(test, val_col = "d", group_col = "Age"))


do_stats_isotype(richness_df)


