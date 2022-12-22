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

new_day = datetime.datetime.now().strftime("%Y%m%d")
path_ori = "../../analysis/plots/"
path = "../../analysis/plots/alpha_isotypes/"

#run_type = "dry"
run_type = "wet"

palette2 = {"DR lifelong":"red", "AL lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL lifelong":"#D3E0F1", "16M AL lifelong":"#A2C2DD", "20M AL lifelong":"#6897C6", "24M AL lifelong":"#3D64A8", "5M DR lifelong":"#ECBAA1", 
            "20M DR lifelong":"#DD694D", "24M DR lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}

iso_palette = {"IGA":"tab:blue", "IGD":"tab:orange", "IGE":"tab:green", "IGG":"tab:red", "IGM":"tab:brown"}
# }}}

files = []
for i in os.listdir(path):
    if os.path.isfile(os.path.join(path,i)) and i.startswith("alpha") and i.endswith("SPL_isosum.tsv") and "stats" not in i:
        files.append(i)
print(files)

df = pd.DataFrame()
for e in files:
    s = e[9:]
    dt = pd.read_csv("{}{}".format(path, e), sep = " ")
    dt["sample_id"] = s[:-15]
    df = pd.concat([df, dt])

# Read metadata so we can get age and diet from there
nam = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8"])
nam.columns = ["Mouse", "Diet", "Age", "Illumina", "Illumina2", "Barcode", "sample_id"]
iso = pd.merge(nam.loc[:, ["Diet", "Age", "sample_id"]], df, on = "sample_id")
iso["biogroup"] = iso["Age"].astype(str)+"M "+iso["Diet"]
iso.head()

# {{{
ls_biogroup = iso["biogroup"].unique().tolist()

for e in ls_biogroup:

    fig, ax = plt.subplots(figsize = (11,6))
    sns.lineplot(data = iso[iso["biogroup"]== e], x = "q", y = "d", hue = "sum_iso", ax = ax, linewidth = 2.5, palette = iso_palette)
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Diversity order (q)", fontsize = 20)
    ax.set_ylabel("Alpha-diversity (qD)", fontsize = 20)
    plt.title("{}\n".format(e), fontsize = 20)
    plt.ylim(0, 1400)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(2.5)

    plt.axvline(x = 0, alpha = 0.5, color = "black", linestyle = "--")
    plt.axvline(x = 1, alpha = 0.5, color = "black", linestyle = "--")
    plt.axvline(x = 2, alpha = 0.5, color = "black", linestyle = "--")
    plt.text(x = -0.12, y = 80, s = "Richness", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
    plt.text(x = 0.88, y = 80, s = "Shannon", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
    plt.text(x = 1.88, y = 80, s = "Simpson", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)

    str_vals = [0, 1, 2, 3]
    ax.set_xticks([0.0, 1.0, 2.0, 3.0])
    ax.set_xticklabels(str_vals)


    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig("{}SPL_Alpha_ISO-Hill_{}_{}.pdf".format(path, e, new_day))
    else:
        plt.show()
# }}}

ls_biogroup = iso["biogroup"].unique().tolist()
for a in ls_biogroup:
    ls_samp = iso[iso["biogroup"] == a]["sample_id"].unique().tolist()
    
    fig, ax = plt.subplots(figsize = (11,6))

    for e in ls_samp:
        sns.lineplot(data = iso[iso["sample_id"]== e], x = "q", y = "d", hue = "sum_iso", ax = ax, linewidth = 2.5, palette = iso_palette)
        
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Diversity order (q)", fontsize = 20)
    ax.set_ylabel("Alpha-diversity (qD)", fontsize = 20)
    plt.title("{}\n".format(a), fontsize = 20)
    plt.ylim(0, 1400)
    
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[:6], labels[:6], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(2.5)

    plt.axvline(x = 0, alpha = 0.5, color = "black", linestyle = "--")
    plt.axvline(x = 1, alpha = 0.5, color = "black", linestyle = "--")
    plt.axvline(x = 2, alpha = 0.5, color = "black", linestyle = "--")
    plt.text(x = -0.12, y = 80, s = "Richness", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
    plt.text(x = 0.88, y = 80, s = "Shannon", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
    plt.text(x = 1.88, y = 80, s = "Simpson", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)

    str_vals = [0, 1, 2, 3]
    ax.set_xticks([0.0, 1.0, 2.0, 3.0])
    ax.set_xticklabels(str_vals)


    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig("{}SPL_Alpha_ISO-Hill_Biorep-{}_{}.pdf".format(path, a, new_day))
    else:
        plt.show()

ls_isotype = iso["sum_iso"].unique().tolist()
ls_age = iso["Age"].unique().tolist()
for e in ls_age:
    for a in ls_isotype:
        sm = iso[iso["Age"] == e]
        fig, ax = plt.subplots(figsize = (11, 6))
        sns.lineplot(data = sm[sm["sum_iso"] == a], x = "q", y = "d", hue = "Diet", ax = ax, linewidth = 2.5, palette = palette2)
        ax.tick_params(axis = "x", labelsize=18)
        ax.tick_params(axis = "y", labelsize=18)
        ax.set_xlabel("Diversity order (q)", fontsize = 20)
        ax.set_ylabel("Alpha-diversity (qD)", fontsize = 20)
        plt.title("{}M - {}\n".format(e, a), fontsize = 20)

        leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        for line in leg.get_lines():
            line.set_linewidth(2.5)

        plt.axvline(x = 0, alpha = 0.5, color = "black", linestyle = "--")
        plt.axvline(x = 1, alpha = 0.5, color = "black", linestyle = "--")
        plt.axvline(x = 2, alpha = 0.5, color = "black", linestyle = "--")
        plt.text(x = -0.12, y = 80, s = "Richness", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
        plt.text(x = 0.88, y = 80, s = "Shannon", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)
        plt.text(x = 1.88, y = 80, s = "Simpson", rotation = 90, fontsize = 18, color = "black", alpha = 0.5)

        str_vals = [0, 1, 2, 3]
        ax.set_xticks([0.0, 1.0, 2.0, 3.0])
        ax.set_xticklabels(str_vals)


        matplotlib.rcParams['pdf.fonttype'] = 42
        plt.tight_layout()
        if run_type != "dry":
            plt.savefig("{}SPL_Alpha_ISO-Hill_{}M-{}_{}.pdf".format(path, e, a, new_day))
        else:
            plt.show()

iso.head()

# Make a smaller dataframe to work easily with the independent alpha values
shannon_df = iso[iso["q"] == 1.0]
simpson_df = iso[iso["q"] == 2.0]
richness_df = iso[iso["q"] == 0.0]


def box_alpha_iso(df, text):
    for a in ls_isotype:

        fig, ax = plt.subplots(figsize = (13,5))
        ax = sns.boxplot(data = df[df["sum_iso"] == a], x = "Age", y = "d", hue = "Diet", palette = palette2, showfliers = False, ax = ax, hue_order = ["AL lifelong", "DR lifelong", "AL_DR16M", "AL_DR20M"])
        ax = sns.swarmplot(data = df[df["sum_iso"] == a], x = "Age", y = "d", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL lifelong", "DR lifelong", "AL_DR16M", "AL_DR20M"])
        ax.tick_params(axis = "x", labelsize=18)
        ax.tick_params(axis = "y", labelsize=18)
        ax.set_xlabel("Age [Months]", fontsize = 20)
        ax.set_ylabel("{} - Alpha Diversity".format(text), fontsize = 20)
        plt.title("{}\n".format(a), fontsize = 20)


        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

        matplotlib.rcParams['pdf.fonttype'] = 42
        plt.tight_layout()

        if run_type != "dry":
            plt.savefig("{}SPL_{}_{}_{}.pdf".format(path, text, a, new_day))
        else:

            plt.show()

box_alpha_iso(shannon_df, "Shannon")
box_alpha_iso(simpson_df, "Simpson")
box_alpha_iso(richness_df, "Richness")


def do_stats_isotype(da, ig):
    
    df = da[da["sum_iso"] == ig]
    print("CHECKING {}\n".format(ig))
    print(stats.kruskal(df["d"].to_list(), df["biogroup"].to_list()))
    display(sp.posthoc_conover(df, val_col = "d", group_col = "biogroup"))
    
    print("Anova AL/DR\n")
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = df[df["Diet"].isin(["AL lifelong", "DR lifelong"])]).fit(), typ=2))
    
    print("ANOVA with AL before switch\n")
    print("Anova AL/AL_DR16M\n")
    
    sw16 = df[df["biogroup"].isin(["5M AL lifelong", "16M AL lifelong", "20M AL_DR16M", "24M AL_DR16M"])]
    sw16["Diet"] = "AL_DR16M"
    AL_sw16 = pd.concat([sw16, df[df["Diet"] == "AL lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw16).fit(), typ=2))
    
    print("Anova AL/AL_DR20M\n")
    sw20 = df[df["biogroup"].isin(["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL_DR20M"])]
    sw20["Diet"] = "AL_DR20M"
    AL_sw20 = pd.concat([sw20, df[df["Diet"] == "AL lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = AL_sw20).fit(), typ=2))
    
    print("Anova AL_DR16M/DR\n")
    DR_sw16 = pd.concat([sw16, df[df["Diet"] == "DR lifelong"]])
    display(anova_lm(ols("d ~ Age * Diet + Age + Diet", data = DR_sw16).fit(), typ=2))
    
    print("Anova AL_DR20M/DR\n")
    DR_sw20 = pd.concat([sw20, df[df["Diet"] == "DR lifelong"]])
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
    test = df[df["Diet"] == "AL lifelong"]

    print(stats.kruskal(test["d"].to_list(), test["Age"].to_list()))

    display(sp.posthoc_conover(test, val_col = "d", group_col = "Age"))

    print("DR lifelong stats: \n")
    test = df[df["Diet"] == "DR lifelong"]

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


print("\n\nIGA\n")
do_stats_isotype(shannon_df, "IGA")
print("\n\nIGD\n")
do_stats_isotype(shannon_df, "IGD")
print("\n\nIGE\n")
do_stats_isotype(shannon_df, "IGE")
print("\n\nIGG\n")
do_stats_isotype(shannon_df, "IGG")
print("\n\nIGM\n")
do_stats_isotype(shannon_df, "IGM")

print("\n\nIGA\n")
do_stats_isotype(simpson_df, "IGA")
print("\n\nIGD\n")
do_stats_isotype(simpson_df, "IGD")
print("\n\nIGE\n")
do_stats_isotype(simpson_df, "IGE")
print("\n\nIGG\n")
do_stats_isotype(simpson_df, "IGG")
print("\n\nIGM\n")
do_stats_isotype(simpson_df, "IGM")

print("\n\nIGA\n")
do_stats_isotype(richness_df, "IGA")
print("\n\nIGD\n")
do_stats_isotype(richness_df, "IGD")
print("\n\nIGE\n")
do_stats_isotype(richness_df, "IGE")
print("\n\nIGG\n")
do_stats_isotype(richness_df, "IGG")
print("\n\nIGM\n")
do_stats_isotype(richness_df, "IGM")

palette_new = {"5M AL_lifelong":"dodgerblue", "16M AL_lifelong":"dodgerblue", "20M AL_lifelong":"dodgerblue", "24M AL_lifelong":"dodgerblue", "5M DR_lifelong":"red", 
            "20M DR_lifelong":"red", "24M DR_lifelong":"red","20M AL_DR16M":"tab:green", "24M AL_DR16M":"tab:green", "24M AL_DR20M":"yellow"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}


def boxplots_cond(df_prov, order, str_vals, text, a):
    
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
        plt.ylim(0, 250)
    else:
        exit(1)
    
    nam = order[-1].split(" ")[1]
    plt.title(nam + " - " + a, fontsize = 20)
    
    #Linear regression
    x = [float(w.split("M")[0]) for w in df["biogroup"].to_list()]
    y = df.d.to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    
    if run_type == "dry":
        plt.savefig("{}SPL_{}Box_{}_{}_{}.pdf".format(path, text, a, nam, new_day))
        #plt.show()
    else:
        plt.show()


# ### Shannon linear regression

for a in ls_isotype:
    order = ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(shannon_df[shannon_df["sum_iso"] == a], order, str_vals, "Shannon", a)

# {{{
for a in ls_isotype:
    order = ["5M DR lifelong", "16M AL lifelong", "20M DR lifelong", "24M DR lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(shannon_df[shannon_df["sum_iso"] == a], order, str_vals, "Shannon", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(shannon_df[shannon_df["sum_iso"] == a], order, str_vals, "Shannon", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL_DR16M", "24M AL_DR16M"]
    str_vals = [5, 20, 24]
    boxplots_cond(shannon_df[shannon_df["sum_iso"] == a], order, str_vals, "Shannon", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL_DR20M"]
    str_vals = [5, 20, 24]
    boxplots_cond(shannon_df[shannon_df["sum_iso"] == a], order, str_vals, "Shannon", a)
# }}}

# ### Simpson linear regression

for a in ls_isotype:
    order = ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(simpson_df[simpson_df["sum_iso"] == a], order, str_vals, "Simpson", a)

# {{{
for a in ls_isotype:
    order = ["5M DR lifelong", "16M AL lifelong", "20M DR lifelong", "24M DR lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(simpson_df[simpson_df["sum_iso"] == a], order, str_vals, "Simpson", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(simpson_df[simpson_df["sum_iso"] == a], order, str_vals, "Simpson", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL_DR16M", "24M AL_DR16M"]
    str_vals = [5, 20, 24]
    boxplots_cond(simpson_df[simpson_df["sum_iso"] == a], order, str_vals, "Simpson", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL_DR20M"]
    str_vals = [5, 20, 24]
    boxplots_cond(simpson_df[simpson_df["sum_iso"] == a], order, str_vals, "Simpson", a)
# }}}

# ### Richness linear regression

for a in ls_isotype:
    order = ["5M DR lifelong", "20M DR lifelong", "24M DR lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(richness_df[richness_df["sum_iso"] == a], order, str_vals, "Richness", a)

# {{{
for a in ls_isotype:
    order = ["5M DR lifelong", "16M AL lifelong", "20M DR lifelong", "24M DR lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(richness_df[richness_df["sum_iso"] == a], order, str_vals, "Richness", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL lifelong"]
    str_vals = [5, 20, 24]
    boxplots_cond(richness_df[richness_df["sum_iso"] == a], order, str_vals, "Richness", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL_DR16M", "24M AL_DR16M"]
    str_vals = [5, 20, 24]
    boxplots_cond(richness_df[richness_df["sum_iso"] == a], order, str_vals, "Richness", a)
    
for a in ls_isotype:
    order = ["5M AL lifelong", "16M AL lifelong", "20M AL lifelong", "24M AL_DR20M"]
    str_vals = [5, 20, 24]
    boxplots_cond(richness_df[richness_df["sum_iso"] == a], order, str_vals, "Richness", a)
# }}}


