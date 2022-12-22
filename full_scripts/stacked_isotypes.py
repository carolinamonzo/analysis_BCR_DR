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
path = "../../analysis/plots/isotypes_abundance/"

run_type = "dry"
#run_type = "wet"

#organ = "ILE"
organ = "SPL"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"teal", "AL_DR20M":"gold", "AL_DR24M":"grey"}
palette3 = {"5M AL lifelong":"#D3E0F1", "16M AL lifelong":"#A2C2DD", "20M AL lifelong":"#6897C6", "24M AL lifelong":"#3D64A8", "5M DR lifelong":"#ECBAA1", 
            "20M DR lifelong":"#DD694D", "24M DR lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}

iso_palette = {"IGA":"tab:blue", "IGD":"tab:orange", "IGE":"tab:green", "IGG":"tab:red", "IGM":"tab:brown"}
# }}}

# Read metadata so we can get age and diet from there
nam = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8"])
nam.columns = ["Mouse", "Diet", "Age", "Illumina", "Illumina2", "Barcode", "sample_id"]

# Read full file to get counts
iso = pd.read_csv("../../analysis/created_germlines/reannotated_db/merged_changeo_{}_isosum.tsv".format(organ), sep = "\t")

# {{{
cou = iso.loc[:, ["sample_id", "biogroup", "sum_iso"]]

# Getting the info we need into our dataframe
diets = []
for e in cou["biogroup"].to_list():
    if e[2:] == "R_lifelong":
        diets.append("DR_lifelong")
    elif e[2:] == "L_lifelong":
        diets.append("AL_lifelong")
    else:
        diets.append(e[2:])

ages = []
for e in cou["biogroup"].to_list():
    ages.append(e[:2])
ages = [w.replace('5D', "5") for w in ages]
ages = [w.replace('5A', "5") for w in ages]

ages = [ int(x) for x in ages ]

# Add to dataframe
cou["Age"] = ages
cou["Diet"] = diets
# }}}

# Get counts so we have an amount per sample
counts = cou.groupby(["biogroup", 'Age','Diet', "sample_id", "sum_iso"]).size()
counts = counts.reset_index()
counts.columns = ["biogroup", "Age", "Diet", "sample_id", "sum_iso", "counts"]
# Now normalize to mean
counts = counts.groupby(["biogroup", "Age", "Diet", "sum_iso"])["counts"].describe()[["mean"]]
counts = counts.reset_index()

# {{{
tostack = counts.pivot(columns="sum_iso", values='mean', index = "biogroup").reset_index()
tostack = tostack.set_index("biogroup")

# Getting the info we need into our dataframe
diets = []
for e in tostack.index.to_list():
    if e[2:] == "R_lifelong":
        diets.append("DR_lifelong")
    elif e[2:] == "L_lifelong":
        diets.append("AL_lifelong")
    else:
        diets.append(e[2:])

ages = []
for e in tostack.index.to_list():
    ages.append(e[:2])
ages = [w.replace('5D', "5") for w in ages]
ages = [w.replace('5A', "5") for w in ages]

ages = [ int(x) for x in ages ]

# Add to dataframe
tostack["Age"] = ages
tostack["Diet"] = diets
# }}}

# Get counts so we have an amount per sample
counts_samp = cou.groupby(["biogroup", 'Age','Diet', "sample_id", "sum_iso"]).size()
counts_samp = counts_samp.reset_index()
counts_samp.columns = ["biogroup", "Age", "Diet", "sample_id", "sum_iso", "counts"]

# {{{
disamp = {}

for e in counts_samp["sample_id"].unique():
    disamp[e] = counts_samp[counts_samp["sample_id"] == e]["counts"].sum()
    
# Map this to the dataframe
counts_samp["Total"] = counts_samp["sample_id"].map(disamp)
# Make this as percentages for plotting
counts_samp["perc"] = (counts_samp["counts"].astype(float) / counts_samp["Total"].astype(float)) * 100
# }}}

counts_samp[counts_samp["biogroup"] == "5AL_lifelong"].groupby(["sum_iso"])["perc"].mean()


def plot_perc_iso(isot):
    fig, ax = plt.subplots(figsize = (10,5))

    sns.boxplot(data = counts_samp[counts_samp["sum_iso"] == isot.upper()], x = "Age", y = "perc", hue = "Diet", palette = palette2, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], showfliers = False)
    sns.swarmplot(data = counts_samp[counts_samp["sum_iso"] == isot.upper()], x = "Age", y = "perc", hue = "Diet", hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], dodge = True, color = ".25")

    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("{} frequency [%]".format(isot), fontsize = 20)
    ax.get_legend().remove()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig("{}{}_{}_frequencyBox_{}.pdf".format(path, organ, isot, new_day))
    else:
        plt.show()


plot_perc_iso("IgA")
plot_perc_iso("IgG")
plot_perc_iso("IgM")
plot_perc_iso("IgE")
plot_perc_iso("IgD")





# {{{
def stacked_percent_df(df, diets):
    r = [0,1,2,3,4]

    # From raw value to percentage
    totals = [i+j+k+h+u for i,j,k,h,u in zip(df['IGM'], df['IGG'], df['IGA'], df["IGD"], df["IGE"])]
    Bars0 = [i / j * 100 for i,j in zip(df['IGM'], totals)]
    Bars1 = [i / j * 100 for i,j in zip(df['IGG'], totals)]
    Bars2 = [i / j * 100 for i,j in zip(df['IGA'], totals)]
    Bars3 = [i / j * 100 for i,j in zip(df['IGD'], totals)]
    Bars4 = [i / j * 100 for i,j in zip(df['IGE'], totals)]

    # Make dataframe with data, percentage of mice with x number of tumours per diet
    #st = pd.DataFrame([Bars0, Bars1, Bars2, Bars3, Bars4], columns = ["AL", "DR", "12M", "16M", "20M", "24M"])
    #st = pd.DataFrame([Bars0, Bars1, Bars2, Bars3, Bars4], columns = ["AL", "DR", "12M", "16M", "20M"])
    st = pd.DataFrame([Bars0, Bars1, Bars2, Bars3, Bars4], columns = diets)
    st.index = ["IGM", "IGG", "IGA", "IGD", "IGE"]
    return(st)

def make_staked_tumourload(df, text, figsize):
    # Stacked bars of tumour presence at time of death
    fig, ax = plt.subplots(figsize = figsize)

    lisa = ["tab:brown", "tab:red", "tab:blue", "tab:orange", "tab:green"]

    df.T.plot(kind='bar', stacked = True, ax = ax, color = lisa, width = 0.85, edgecolor = "black", 
             linewidth = 1.5)
    plt.xticks(rotation = 0)

    ax.tick_params(axis = "x", labelsize=20)
    ax.tick_params(axis = "y", labelsize=20)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
    plt.title("{}".format(text), fontsize = 18)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    m = text.split(" ")[0]
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig("{}{}_stacked_proportion_{}M_{}.pdf".format(path, organ, m, new_day))
    else:

        plt.show()
# }}}


groups = ["24DR_lifelong", "24AL_DR16M", "24AL_DR20M", "24AL_lifelong"]
toplot = tostack.loc[groups, :]
#diets = toplot["Diet"].tolist()
diets = ["DR", "16.sw", "20.sw", "AL"]
st_all = stacked_percent_df(toplot, diets)
make_staked_tumourload(st_all, "24 months of age", (7,7))

groups = ["20DR_lifelong", "20AL_DR16M", "20AL_lifelong"]
toplot = tostack.loc[groups, :]
#diets = toplot["Diet"].tolist()
diets = ["DR", "16.sw", "AL"]
st_all = stacked_percent_df(toplot, diets)
make_staked_tumourload(st_all, "20 months of age", (6,7))

groups = ["5DR_lifelong", "5AL_lifelong"]
toplot = tostack.loc[groups, :]
#diets = toplot["Diet"].tolist()
diets = ["DR", "AL"]
st_all = stacked_percent_df(toplot, diets)
make_staked_tumourload(st_all, "5 months of age", (5,7))

groups = ["16AL_lifelong"]
toplot = tostack.loc[groups, :]
#diets = toplot["Diet"].tolist()
diets = ["16"]
st_all = stacked_percent_df(toplot, diets)
make_staked_tumourload(st_all, "16 months of age", (3.9,7))

counts_samp = cou.groupby(["biogroup", 'Age','Diet', "sample_id", "sum_iso"]).size()
counts_samp = counts_samp.reset_index()
counts_samp.columns = ["biogroup", "Age", "Diet", "sample_id", "sum_iso", "counts"]

# {{{
samples = list(counts_samp["sample_id"].unique())
tot = []
a = []

for sam in samples:
    val = counts_samp[counts_samp["sample_id"]== sam]["counts"].sum()
    a.append(val.astype(str))
    tot = tot + (a * counts_samp[counts_samp["sample_id"]== sam].shape[0])
    a = []
    
counts_samp["Total"] = tot
# }}}

counts_samp["rel"] = (counts_samp["counts"] / counts_samp["Total"].astype(float))*100

# {{{
desc = {"IGA":"classwitched", "IGG":"classwitched", "IGE":"classwitched", "IGM":"naive", "IGD":"naive"}
def check_desc(x):
    for key in desc:
        if key.lower() in x.lower():
            return(desc[key])
    return ''

counts_samp["status"] = counts_samp["sum_iso"].map(lambda x: check_desc(x))
# }}}



counts_samp_status = counts_samp.groupby(["status", "sample_id", "biogroup", "Age", "Diet"]).sum().reset_index()
counts_samp_mer = counts_samp_status.pivot(index = "sample_id", columns = "status", values = "counts").merge(counts_samp_status.loc[:, ["sample_id", "biogroup", "Age", "Diet"]], left_index = True, right_on = "sample_id")
counts_samp_mer["rel_csn"] = counts_samp_mer["naive"].astype(float) / counts_samp_mer["classwitched"].astype(float)
counts_samp_mer.drop_duplicates(inplace = True)

# {{{
fig, ax = plt.subplots(figsize = (10,5))

sns.boxplot(data = counts_samp_mer, x = "Age", y = "rel_csn", hue = "Diet", palette = palette2, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], showfliers = False)
sns.swarmplot(data = counts_samp_mer, x = "Age", y = "rel_csn", hue = "Diet", hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], dodge = True, color = ".25")

ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Age [Months]", fontsize = 20)
ax.set_ylabel("Relative abundance [naive/class-switched]", fontsize = 14)
ax.get_legend().remove()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig("{}{}_naiveCSrelab_boxplot_{}.pdf".format(path, organ, new_day))
else:
    plt.show()
# }}}

print("//K-W BETWEEN TREATMENTS WITHIN TIMEPOINTS//")
for e in [5, 20, 24]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = counts_samp_mer[counts_samp_mer["Age"] == e]

    print(stats.kruskal(test["rel_csn"].to_list(), test["Diet"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months\n".format(e))
    print(sp.posthoc_mannwhitney(test, val_col = "rel_csn", group_col = "Diet", p_adjust = 'fdr_bh'))

print("//K-W BETWEEN timepoints WITHIN treatments//")
for e in ["AL_lifelong", "DR_lifelong", "AL_DR16M"]:
    print("Kruskall-Wallis comparing treatments at {} Months\n".format(e))
    test = counts_samp_mer[counts_samp_mer["Diet"] == e]

    print(stats.kruskal(test["rel_csn"].to_list(), test["Age"].to_list()))
    print("Mann-Whitney comparing treatments at {} Months\n".format(e))
    print(sp.posthoc_mannwhitney(test, val_col = "rel_csn", group_col = "Age", p_adjust = 'fdr_bh'))


def boxplots_cond(df_prov, order, method):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup"].isin(order)]
    
    #Linear regression
    x = df.Age.to_list()
    y = df[method].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")


print("\nAL")
boxplots_cond(counts_samp_mer, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
                 "rel_csn")
print("\nDR")
boxplots_cond(counts_samp_mer, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
                 "rel_csn")
print("\nAL_DR16M")
boxplots_cond(counts_samp_mer, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
                 "rel_csn")
print("\nAL_DR20M")
boxplots_cond(counts_samp_mer, ["20AL_lifelong", "24AL_DR20M"], 
                 "rel_csn")

# {{{
dfs_small = counts_samp_mer[counts_samp_mer["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

model = ols('rel_csn ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
display(sm.stats.anova_lm(model, typ=2))
# }}}

# {{{
dfs_small = counts_samp_mer[counts_samp_mer["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
dfs_small = dfs_small[dfs_small["Age"] > 17]

model = ols('rel_csn ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
display(sm.stats.anova_lm(model, typ=2))
# }}}

# {{{
dfs_small = counts_samp_mer[counts_samp_mer["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
dfs_small = dfs_small[dfs_small["Age"] > 17]

model = ols('rel_csn ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
display(sm.stats.anova_lm(model, typ=2))
# }}}

forbar = counts_samp_mer.groupby("biogroup").median().reset_index()

# {{{
AL_life16 = forbar[forbar["biogroup"] == "16AL_lifelong"]["rel_csn"]
AL_life20 = forbar[forbar["biogroup"] == "20AL_lifelong"]["rel_csn"]
AL16DR20 = forbar[forbar["biogroup"] == "20AL_DR16M"]["rel_csn"]

dicforbar = pd.DataFrame({"AL_lifelong": [(float(AL_life20) - float(AL_life16))], "AL_DR16M":[(float(AL16DR20)- float(AL_life16))]})
# }}}

# {{{
fig, ax = plt.subplots(figsize = (4,5))
sns.barplot(data = dicforbar, palette = palette2)

ax.tick_params(axis = "x", labelsize=14)
ax.tick_params(axis = "y", labelsize=16)
ax.set_xlabel("")
ax.set_ylabel("\u0394 (RelAbund 20M - RelAbund 16M)", fontsize = 17)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.title("[naive/class-switched]\n", fontsize = 18)

plt.ylim(-3, 3)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig("{}{}_delta16-20_naiveCS_{}.pdf".format(path, organ, new_day))
else:
    plt.show()
# }}}

# {{{
AL_life20 = forbar[forbar["biogroup"] == "20AL_lifelong"]["rel_csn"]
AL_life24 = forbar[forbar["biogroup"] == "24AL_lifelong"]["rel_csn"]
AL20DR24 = forbar[forbar["biogroup"] == "24AL_DR20M"]["rel_csn"]

dicforbar2 = pd.DataFrame({"AL_lifelong": [(float(AL_life24) - float(AL_life20))], "AL_DR20M":[(float(AL20DR24)- float(AL_life20))]})
# }}}

# {{{
fig, ax = plt.subplots(figsize = (4,5))
sns.barplot(data = dicforbar2, palette = palette2)

ax.tick_params(axis = "x", labelsize=14)
ax.tick_params(axis = "y", labelsize=16)
ax.set_xlabel("")
ax.set_ylabel("\u0394 (RelAbund 24M - RelAbund 20M)", fontsize = 17)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.title("[naive/class-switched]\n", fontsize = 18)

plt.ylim(-3, 3)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig("{}{}_delta20-24_naiveCS_{}.pdf".format(path, organ, new_day))
else:
    plt.show()
# }}}

# {{{
## Check stats

import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
stats = importr('stats')
m = np.array([[(float(dicforbar["AL_DR16M"]))+20, (float(dicforbar["AL_lifelong"]))+20], [(float(dicforbar2["AL_DR20M"]))+20, (float(dicforbar2["AL_lifelong"]))+20]])
res = stats.fisher_test(m)
print("p-value change in young switch vs old switch: {}".format(res[0][0]))
# }}}

m


def barplots_iso(d, t, iso):
    
    fig, ax = plt.subplots(figsize = (4,5))
    sns.barplot(data = d, palette = palette2)

    ax.tick_params(axis = "x", labelsize=14)
    ax.tick_params(axis = "y", labelsize=16)
    ax.set_xlabel("")
    ax.set_ylabel("\u0394 (RelAbund {}M - RelAbund {}M)".format(t[1], t[0]), fontsize = 17)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.title("{}\n".format(iso), fontsize = 18)

    if organ == "SPL":
        if iso in ["IGA"]:
            plt.ylim(-4.5, 1.5)
        elif iso in ["IGD"]:
            plt.ylim(-2.1, 1.1)
        elif iso in ["IGE"]:
            plt.ylim(-0.18, 0.14)
        elif iso in ["IGG"]:
            plt.ylim(-19, 21)
        else:
            plt.ylim(-24, 20)
    elif organ == "ILE":
        if iso in ["IGA"]:
            plt.ylim(-7, 9)
        elif iso in ["IGD"]:
            plt.ylim(-0.5, 0.02)
        elif iso in ["IGE"]:
            plt.ylim(-0.20, 0.14)
        elif iso in ["IGG"]:
            plt.ylim(-1.6, 1.6)
        else:
            plt.ylim(-3, 0)
    
    time = "-".join(t)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig("{}{}_{}delta{}_bar_{}.pdf".format(path, organ, iso, time, new_day))
    else:
        plt.show()


# ## Repeating this for all isotypes

# {{{
igA = counts_samp[counts_samp['sum_iso'] == "IGA"].groupby("biogroup").median().reset_index()
igD = counts_samp[counts_samp['sum_iso'] == "IGD"].groupby("biogroup").median().reset_index()
igE = counts_samp[counts_samp['sum_iso'] == "IGE"].groupby("biogroup").median().reset_index()
igG = counts_samp[counts_samp['sum_iso'] == "IGG"].groupby("biogroup").median().reset_index()
igM = counts_samp[counts_samp['sum_iso'] == "IGM"].groupby("biogroup").median().reset_index()

tup_dfs = [(igA, "IGA"), (igD, "IGD"), (igE, "IGE"), (igG, "IGG"), (igM, "IGM")]
# }}}

# {{{
t = ["16", "20"]
lateryoung = []
for d, iso in tup_dfs:
    AL_life16 = d[d["biogroup"] == "16AL_lifelong"]["rel"]
    AL_life20 = d[d["biogroup"] == "20AL_lifelong"]["rel"]
    AL16DR20 = d[d["biogroup"] == "20AL_DR16M"]["rel"]

    dicforbar = pd.DataFrame({"AL_lifelong": [(float(AL_life20) - float(AL_life16))], "AL_DR16M":[(float(AL16DR20)- float(AL_life16))]})
    
    lateryoung.append([iso] + dicforbar.values.tolist()[0])
    
    barplots_iso(dicforbar, t, iso)
    
lateryoung = pd.DataFrame(lateryoung)
lateryoung.columns = ["isotype", "AL_lifelong", "AL_DR16M"]
# }}}



t = ["20", "24"]
laterold = []
for d, iso in tup_dfs:
    AL_life20 = d[d["biogroup"] == "20AL_lifelong"]["rel"]
    AL_life24 = d[d["biogroup"] == "24AL_lifelong"]["rel"]
    AL20DR24 = d[d["biogroup"] == "24AL_DR20M"]["rel"]

    dicforbar2 = pd.DataFrame({"AL_lifelong": [(float(AL_life24) - float(AL_life20))], "AL_DR20M":[(float(AL20DR24)- float(AL_life20))]})
    
    barplots_iso(dicforbar2, t, iso)
    
    laterold.append([iso] + dicforbar2.values.tolist()[0])
laterold = pd.DataFrame(laterold)
laterold.columns = ["isotype", "AL_lifelong", "AL_DR20M"]

for e in ["IGA", "IGD", "IGE", "IGG", "IGM"]:
    
    baryoung = lateryoung[lateryoung["isotype"] == e]
    barold = laterold[laterold["isotype"] == e]
    
    m = np.array([[(float(baryoung["AL_DR16M"].iloc[0]))+23, (float(baryoung["AL_lifelong"].iloc[0]))+23], [(float(barold["AL_DR20M"].iloc[0]))+23, (float(barold["AL_lifelong"].iloc[0]))+23]])
    print(m)
    res = stats.fisher_test(m)
    print("{} p-value change in young switch vs old switch: {}".format(e, res[0][0]))

# ## Relative abundance through age and diet

# {{{
igA = counts_samp[counts_samp['sum_iso'] == "IGA"]
igD = counts_samp[counts_samp['sum_iso'] == "IGD"]
igE = counts_samp[counts_samp['sum_iso'] == "IGE"]
igG = counts_samp[counts_samp['sum_iso'] == "IGG"]
igM = counts_samp[counts_samp['sum_iso'] == "IGM"]

tup_dfs = [(igA, "IGA"), (igD, "IGD"), (igE, "IGE"), (igG, "IGG"), (igM, "IGM")]
# }}}

for d, iso in tup_dfs:
    fig, ax = plt.subplots(figsize = (10,5))

    sns.boxplot(data = d, x = "Age", y = "rel", hue = "Diet", palette = palette2, hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], showfliers = False)
    sns.swarmplot(data = d, x = "Age", y = "rel", hue = "Diet", hue_order = ["AL_lifelong", "DR_lifelong", "AL_DR16M", "AL_DR20M"], dodge = True, color = ".25")

    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("Relative abundance % [{}]".format(iso), fontsize = 18)
    ax.get_legend().remove()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig("{}{}_{}_relAbund_{}.pdf".format(path, organ, iso, new_day))
    else:
        plt.show()

# {{{
from scipy import stats
method = "rel"

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
    print("AL")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_lifelong", "24AL_lifelong"], 
                 "rel")
    print("\nDR")
    boxplots_cond(d, ["5DR_lifelong", "20DR_lifelong", "24DR_lifelong"], 
                 "rel")
    print("\nAL_DR16M")
    boxplots_cond(d, ["5AL_lifelong", "16AL_lifelong", "20AL_DR16M", "24AL_DR_16M"], 
                 "rel")
    print("\nAL_DR20M")
    boxplots_cond(d, ["20AL_lifelong", "24AL_DR20M"], 
                 "rel")

for d,iso in tup_dfs:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "DR_lifelong"])]

    model = ols('rel ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["AL_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('rel ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs:
    print("\n\n" + iso)
    dfs_small = d[d["Diet"].isin(["DR_lifelong", "AL_DR16M"])]
    dfs_small = dfs_small[dfs_small["Age"] > 17]

    model = ols('rel ~ C(Age) + C(Diet) + C(Age):C(Diet)', data=dfs_small).fit()
    display(sm.stats.anova_lm(model, typ=2))

for d,iso in tup_dfs:
    d.to_csv(f"{path}{organ}_{iso}_freqValues.csv", sep = ";", index = False)








