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
from natsort import natsorted
import math

new_day = datetime.datetime.now().strftime("%Y%m%d")
path = "../../analysis/plots/cdr3/"

run_type = "dry"
#run_type = "wet"

#organ = "ILE"
organ = "SPL"

palette2 = {"DR lifelong":"red", "AL lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"teal", "AL_DR20M":"gold", "AL_DR24M":"grey"}
# }}}

# Reading metadata file
eD = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";", usecols = ["Mouse", "Diet", 
                "Age", "File"])

# {{{
# Read Spleen germline annotated
ile = pd.read_csv("../../analysis/created_germlines/clone-germline-{}.tsv".format(organ), sep = "\t",
                  usecols = ["junction_length", 
                             "v_call_genotyped", "germline_d_call", "germline_j_call File", 
                             "isotype", "clone_id"])

# Read Ileum germline annotated
#ile = pd.read_csv("../created_germlines/clone-germline-ILE.tsv", sep = "\t",
#                  usecols = ["junction_length", 
#                             "v_call_genotyped", "germline_d_call", "germline_j_call File", 
#                             "isotype", "clone_id"])
                  
ile[['germline_j_call','File']] = ile['germline_j_call File'].str.split(' ',expand=True)
ile = ile.drop(columns = ['germline_j_call File'])
# }}}

ile = pd.merge(eD, ile, on = ["File"], how = "inner")
ile["CDR3aa_length"] = ile["junction_length"]/3
ile["biorep"] = ile["Age"].astype(str) + ile["Diet"]

ile.head()


def plot_CDR3_biorep(df, iso):
    
    fig, ax = plt.subplots(figsize = (8,5))
    gr = ['5AL lifelong', '16AL lifelong', '20AL lifelong', '24AL lifelong', '5DR lifelong', '20DR lifelong', '24DR lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
    co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
    for i in range(0, len(gr)):

        sns.distplot(df[df["biorep"] == gr[i]].loc[:, ["CDR3aa_length"]], hist = False, kde = False, fit = norm, 
                     fit_kws={"color":co[i], "linewidth":3}, label = gr[i])
        (mu, sigma) = norm.fit(df[df["biorep"] == gr[i]].loc[:, ["CDR3aa_length"]])
        print("{}: mu={}, sigma={}".format(gr[i], mu, sigma))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tick_params(axis = "x", labelsize = 18)
    plt.tick_params(axis = "y", labelsize = 18)
    ax.set_xlabel("CDR3 length (AA)", fontsize = 16)
    ax.set_ylabel("Frequency", fontsize = 16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 16)
    plt.xlim(0, 30)
    matplotlib.rcParams['pdf.fonttype'] = 42                 
    
    plt.tight_layout()
    
    if run_type != "dry":
        plt.savefig("{}{}_cdr3_{}_biogroup_{}.pdf".format(path, organ, iso, new_day))
    else:
        plt.show()


# {{{
dic_rep = {"IGG12":"IGG", "IGG3":"IGG"}
ile.replace({"isotype":dic_rep}, inplace = True)

igA = ile[ile["isotype"] == "IGA"]
igE = ile[ile["isotype"] == "IGE"]
igD = ile[ile["isotype"] == "IGD"]
igG = ile[ile["isotype"] == "IGG"]
igM = ile[ile["isotype"] == "IGM"]

dfs = [(igA, "IGA"), (igE, "IGE"), (igD, "IGD"), (igG, "IGG"), (igM, "IGM")]
# }}}

for d, iso in dfs:
    print(iso)
    plot_CDR3_biorep(d, iso)
# The MU is basically telling us the center of the gaussian normal distribution we are using
# They are basically all the same except AL at 24 months

# {{{
def cdr3_samp(ile, eD):
    #Fit CDR3 length distribution on a normal distribution, calculate mu (center) and sigma (standard deviation)
    lmouse = []
    lmu = []
    lsigma = []

    for e in ile["Mouse"].unique():
        (mu, sigma) = norm.fit(ile[ile["Mouse"] == e].loc[:, ["CDR3aa_length"]])
        lmouse.append(e)
        lmu.append(mu)
        lsigma.append(sigma)

    dens_samp = pd.DataFrame({"Mouse":lmouse, "Mu":lmu, "Sigma":lsigma})
    df_dens_mouse = pd.merge(eD.loc[:, ["Mouse", "Diet", "Age"]].drop_duplicates(), dens_samp, on = ["Mouse"], how = "inner")
    df_dens_mouse["biorep"] = df_dens_mouse["Age"].astype(str) + df_dens_mouse["Diet"]
    return(df_dens_mouse)

igAsamp = cdr3_samp(igA, eD)
igEsamp = cdr3_samp(igE, eD)
igDsamp = cdr3_samp(igD, eD)
igGsamp = cdr3_samp(igG, eD)
igMsamp = cdr3_samp(igM, eD)

dfs_samp = [(igAsamp, "IGA"), (igEsamp, "IGE"), (igDsamp, "IGD"), (igGsamp, "IGG"), (igMsamp, "IGM")]
# }}}
# {{{
iga = igAsamp.loc[:, ["Mouse", "Mu", "Sigma"]]
ige = igEsamp.loc[:, ["Mouse", "Mu", "Sigma"]]
igd = igDsamp.loc[:, ["Mouse", "Mu", "Sigma"]]
igg = igGsamp.loc[:, ["Mouse", "Mu", "Sigma"]]
igm = igMsamp.loc[:, ["Mouse", "Mu", "Sigma"]]
iga.columns = ["Animal_ID", "cdr3_mu_IGA", "cdr3_sigma_IGA"]
ige.columns = ["Animal_ID", "cdr3_mu_IGE", "cdr3_sigma_IGE"]
igd.columns = ["Animal_ID", "cdr3_mu_IGD", "cdr3_sigma_IGD"]
igg.columns = ["Animal_ID", "cdr3_mu_IGG", "cdr3_sigma_IGG"]
igm.columns = ["Animal_ID", "cdr3_mu_IGM", "cdr3_sigma_IGM"]

igdf = pd.merge(iga, ige, on = "Animal_ID", how = "inner")
igdf = pd.merge(igdf, igd, on = "Animal_ID", how = "inner")
igdf = pd.merge(igdf, igg, on = "Animal_ID", how = "inner")
igdf = pd.merge(igdf, igm, on = "Animal_ID", how = "inner")

igdf.to_csv("../../analysis/results_tables/cdr3_isotypes_{}_{}.csv".format(organ,new_day))
# }}}












# {{{
# Mean of the mu and sigma for each biogroup
igAbiogroup = igAsamp.groupby("biorep").agg({"Mu":np.mean, "Sigma":np.mean}).reset_index()
igEbiogroup = igEsamp.groupby("biorep").agg({"Mu":np.mean, "Sigma":np.mean}).reset_index()
igDbiogroup = igDsamp.groupby("biorep").agg({"Mu":np.mean, "Sigma":np.mean}).reset_index()
igGbiogroup = igGsamp.groupby("biorep").agg({"Mu":np.mean, "Sigma":np.mean}).reset_index()
igMbiogroup = igMsamp.groupby("biorep").agg({"Mu":np.mean, "Sigma":np.mean}).reset_index()

dfs_biogroup = [(igAbiogroup, "IGA"), (igEbiogroup, "IGE"), (igDbiogroup, "IGD"), (igGbiogroup, "IGG"), (igMbiogroup, "IGM")]
# }}}

for df_dens_biogroup, iso in dfs_biogroup:
    fig, ax = plt.subplots(figsize = (5,5))
    gr = ['5AL lifelong', '16AL lifelong', '20AL lifelong', '24AL lifelong', '5DR lifelong', '20DR lifelong', '24DR lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
    co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
    for i in range(0, len(gr)):
        mu = df_dens_biogroup[df_dens_biogroup["biorep"] == gr[i]]["Mu"]
        sigma = df_dens_biogroup[df_dens_biogroup["biorep"] == gr[i]]["Sigma"]

        x = np.linspace(mu - 3*sigma, mu+3*sigma, 100)
        plt.plot(x, stats.norm.pdf(x, mu, sigma), color = co[i], linewidth=3)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tick_params(axis = "x", labelsize = 18)
    plt.tick_params(axis = "y", labelsize = 18)
    ax.set_xlabel("CDR3 length (AA)", fontsize = 16)
    ax.set_ylabel("Frequency", fontsize = 16)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 16)
    plt.xlim(0, 30)
    matplotlib.rcParams['pdf.fonttype'] = 42                 

    plt.tight_layout()

    if run_type != "dry":
        plt.savefig("{}{}_cdr3_{}_biogroup_{}.pdf".format(path, organ, iso, new_day))
    else:
        plt.show()

# Stats to compare these
def compare_curvesKW(df_dens_biogroup, a, b):
    print(a)
    print(b)
    mu1 = df_dens_biogroup[df_dens_biogroup["biorep"] == a]["Mu"]
    sigma1 = df_dens_biogroup[df_dens_biogroup["biorep"] == a]["Sigma"]
    l1 = [item for sublist in np.linspace(mu1 - 3*sigma1, mu1+3*sigma1, 100).tolist() for item in sublist]

    mu2 = df_dens_biogroup[df_dens_biogroup["biorep"] == b]["Mu"]
    sigma2 = df_dens_biogroup[df_dens_biogroup["biorep"] == b]["Sigma"]
    l2 = [item for sublist in np.linspace(mu2 - 3*sigma2, mu2+3*sigma2, 100).tolist() for item in sublist]

    print(stats.ks_2samp(l1, l2))
    print("\\ \n")


for df_dens_biogroup, iso in dfs_biogroup:
    print("//{}//".format(iso))


    biorepl = ['5DR lifelong', '5AL lifelong', '16AL lifelong', '20DR lifelong', 
                '20AL_DR16M', '20AL lifelong', '24DR lifelong', '24AL_DR16M', '24AL_DR20M', '24AL lifelong']

    for i in range(0, len(biorepl)):
        for j in range(1, len(biorepl)):
            if j < len(biorepl):
                compare_curvesKW(df_dens_biogroup, biorepl[i], biorepl[j])

for ile, iso in dfs:

    ile = ile[ile["biorep"] != '16AL lifelong']
    g = sns.FacetGrid(ile, col = "biorep", col_order = ['5DR lifelong', '5AL lifelong', '20DR lifelong', 
                '20AL_DR16M', '20AL lifelong', '24DR lifelong', '24AL_DR16M', '24AL_DR20M', '24AL lifelong'], col_wrap = 5, hue = "File", palette = "husl")
    g = g.map(sns.distplot, "CDR3aa_length", hist = False, kde = False, fit = norm, fit_kws={"linewidth":1.5})
    for ax in g.axes.flat:
        if ax.get_title():
            ax.set_title(ax.get_title().split("=")[1], fontsize = 16)
        ax.tick_params(axis = "x", labelsize = "14")
        ax.tick_params(axis = "y", labelsize = "14")
        ax.set_xlabel(ax.get_xlabel(), fontsize=16)
        ax.set_ylabel(ax.get_ylabel(), fontsize=16)
    plt.xlim(0, 30)
    matplotlib.rcParams['pdf.fonttype'] = 42                 

    plt.tight_layout()

    if run_type != "dry":
        plt.savefig("{}{}_cdr3_{}_allgrid_{}.pdf".format(path, organ, iso, new_day))

    plt.show()

def compare_distr(df, a, b):
    mu1 = df[df["Mouse"] == a]["Mu"]
    sigma1 = df[df["Mouse"] == a]["Sigma"]
    l1 = [item for sublist in np.linspace(mu1 - 3*sigma1, mu1+3*sigma1, 100).tolist() for item in sublist]

    mu2 = df[df["Mouse"] == b]["Mu"]
    sigma2 = df[df["Mouse"] == b]["Sigma"]
    l2 = [item for sublist in np.linspace(mu2 - 3*sigma2, mu2+3*sigma2, 100).tolist() for item in sublist]

    return(stats.ks_2samp(l1, l2)[1])


# {{{
def run_compare_distr(df_dens_mouse):

    biorep = ['20AL_DR16M', '5DR lifelong', '20AL lifelong', '5AL lifelong', '20DR lifelong', '24AL lifelong', '24DR lifelong',
           '24AL_DR16M', '24AL_DR20M']
    # Store number of significant comparisons within each group
    num_sig = 0

    dic_ratios = {"Biogroup":[], "Num_comparisons":[], "Num_significant":[]}

    for j in range(0, len(biorep)):
        df = df_dens_mouse[df_dens_mouse["biorep"] == biorep[j]]
        print(f'\\\\{biorep[j]}\\\\')
        mice = df["Mouse"].to_list()

        for e in [compare_distr(df, mice[0], mice[1]), compare_distr(df, mice[0], mice[2]), compare_distr(df, mice[0], mice[3]), compare_distr(df, mice[0], mice[4]), 
                  compare_distr(df, mice[1], mice[2]), compare_distr(df, mice[1], mice[3]), compare_distr(df, mice[1], mice[4]), compare_distr(df, mice[2], mice[3]), 
                  compare_distr(df, mice[2], mice[4]), compare_distr(df, mice[3], mice[4])]:
            if e < 0.05:
                num_sig = num_sig + 1

        dic_ratios["Biogroup"].append(biorep[j])
        dic_ratios["Num_comparisons"].append(len(mice)*2)
        dic_ratios["Num_significant"].append(num_sig)
        print("Total possible comparisons: {}".format(len(mice)*2))
        print("Significant comparisons:{}".format(num_sig))
        num_sig = 0

    ratios_df = pd.DataFrame(dic_ratios)
    biogroups = ['20AL_DR16M', '5DR lifelong', '20AL lifelong', '5AL lifelong', '16AL lifelong', '20DR lifelong', '24AL lifelong', '24DR lifelong', '24AL_DR16M', '24AL_DR20M']
    Ages = [20, 5, 20, 5, 16, 20, 24, 24, 24, 24]
    Diets = ["AL_DR16M", "DR lifelong", "AL lifelong", "AL lifelong", "AL lifelong", "DR lifelong", "AL lifelong", "DR lifelong", "AL_DR16M", "AL_DR20M"]

    dic_toAge = dict(zip(biogroups, Ages))
    dic_toDiet = dict(zip(biogroups, Diets))

    ratios_df["Age"] = ratios_df["Biogroup"].map(dic_toAge)
    ratios_df["Diet"] = ratios_df["Biogroup"].map(dic_toDiet)
    return(ratios_df)
    
print("\n //IGA//")
igAratio = run_compare_distr(igAsamp)
print("\n //IGE//")
igEratio = run_compare_distr(igEsamp)
print("\n //IGD//")
igDratio = run_compare_distr(igDsamp)
print("\n //IGG//")
igGratio = run_compare_distr(igGsamp)
print("\n //IGM//")
igMratio = run_compare_distr(igMsamp)

dfs_ratio = [(igAratio, "IGA"), (igEratio, "IGE"), (igDratio, "IGD"), (igGratio, "IGG"), (igMratio, "IGM")]
# }}}



for ratios_df, iso in dfs_ratio:
    print("\n //{}//".format(iso))
# For AL and DR check if the linear regression of the slopes 5, 20, 24 is significant
    order = ["5DR lifelong", "20DR lifelong", "24DR lifelong"]
    str_vals = [5, 20, 24]

    x = str_vals
    y = ratios_df[ratios_df["Biogroup"].isin(["5DR lifelong", "20DR lifelong", "24DR lifelong"])].sort_values(by = ["Age"], ascending = True)["Num_significant"].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("DR lifelong")
    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

for ratios_df, iso in dfs_ratio:
    print("\n //{}//".format(iso))
    # For AL and DR check if the linear regression of the slopes 5, 20, 24 is significant
    order = ["5DR lifelong", "20DR lifelong", "24DR lifelong"]
    str_vals = [5, 20, 24]

    x = str_vals
    y = ratios_df[ratios_df["Biogroup"].isin(["5AL lifelong", "20AL lifelong", "24AL lifelong"])].sort_values(by = ["Age"], ascending = True)["Num_significant"].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("AL lifelong")
    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

from statsmodels.stats.anova import anova_lm
for ratios_df, iso in dfs_ratio:
    print("\n //{}//".format(iso))
    # Comparison of the slopes between AL and dR through age
    display(anova_lm(ols("Num_significant ~ Age * Diet + Age + Diet", data = ratios_df[ratios_df["Diet"].isin(["AL lifelong", "DR lifelong"])]).fit(), typ=2))

# {{{
# Fisher to test if within timepoints there are differences in the ratio
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
statsr = importr('stats')

for ratios_df, iso in dfs_ratio:
    print("\n //{}//".format(iso))

    m = ratios_df[ratios_df["Biogroup"].isin(["5AL lifelong", "5DR lifelong"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('5M: p-value AL/DR: {}'.format(res[0][0]))

    m = ratios_df[ratios_df["Biogroup"].isin(["20AL lifelong", "20DR lifelong"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('20M: p-value AL/DR: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["20AL lifelong", "20AL_DR16M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('20M: p-value AL/ALDR16M: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["20DR lifelong", "20AL_DR16M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('20M: p-value DR/ALDR16M: {}'.format(res[0][0]))

    m = ratios_df[ratios_df["Biogroup"].isin(["24AL lifelong", "24DR lifelong"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('24M: p-value AL/DR: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["24AL lifelong", "24AL_DR16M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('24M: p-value AL/ALDR16M: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["24DR lifelong", "24AL_DR16M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('24M: p-value DR/ALDR16M: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["24AL lifelong", "24AL_DR20M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('24M: p-value AL/ALDR20M: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["24DR lifelong", "24AL_DR20M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('24M: p-value DR/ALDR20M: {}'.format(res[0][0]))
    m = ratios_df[ratios_df["Biogroup"].isin(["24AL_DR16M", "24AL_DR20M"])].loc[:, ["Num_significant", "Num_comparisons"]].to_numpy()
    res = statsr.fisher_test(m)
    print('24M: p-value ALDR16M/ALDR20M: {}'.format(res[0][0]))
# }}}
for df_dens_mouse, iso in dfs_samp:

    fig, ax = plt.subplots(figsize = (13,5))
    ax = sns.boxplot(data = df_dens_mouse, x = df_dens_mouse["Age"], y = df_dens_mouse["Mu"], hue = df_dens_mouse["Diet"], hue_order = ["AL lifelong", "DR lifelong", "AL_DR16M", "AL_DR20M"], 
                    palette = {"AL lifelong":"dodgerblue", "DR lifelong":"red", "AL_DR16M":"teal", "AL_DR20M":"gold"})
    ax = sns.swarmplot(data = df_dens_mouse, x = "Age", y = "Mu", hue = "Diet", dodge = True, color=".25", ax = ax, hue_order = ["AL lifelong", "DR lifelong", "AL_DR16M", "AL_DR20M"])
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("Mu of CDR3 length distribution", fontsize = 20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()

    if run_type != "dry":
        plt.savefig("{}{}_MuCDR3_{}_boxplot_{}.pdf".format(path, organ, iso, new_day))
    else:

        plt.show()
for df_dens_mouse, iso in dfs_samp:
    # For AL and DR check if the linear regression of the slopes 5, 20, 24 is significant
    order = ["5DR lifelong", "20DR lifelong", "24DR lifelong"]
    str_vals = [5, 20, 24]

    x = df_dens_mouse[df_dens_mouse["biorep"].isin(["5DR lifelong", "20DR lifelong", "24DR lifelong"])].sort_values(by = ["Age"], ascending = True)["Age"].to_list()
    y = df_dens_mouse[df_dens_mouse["biorep"].isin(["5DR lifelong", "20DR lifelong", "24DR lifelong"])].sort_values(by = ["Age"], ascending = True)["Mu"].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("DR lifelong")
    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")

for df_dens_mouse, iso in dfs_samp:
    # For AL and DR check if the linear regression of the slopes 5, 20, 24 is significant
    order = ["5AL lifelong", "20AL lifelong", "24AL lifelong"]
    str_vals = [5, 20, 24]

    x = df_dens_mouse[df_dens_mouse["biorep"].isin(["5DR lifelong", "20DR lifelong", "24DR lifelong"])].sort_values(by = ["Age"], ascending = True)["Age"].to_list()
    y = df_dens_mouse[df_dens_mouse["biorep"].isin(["5DR lifelong", "20DR lifelong", "24DR lifelong"])].sort_values(by = ["Age"], ascending = True)["Mu"].to_list()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("AL lifelong")
    print(f"Slope: {slope}\nIntercept: {intercept}\np-value: {p_value}\nstd_err: {std_err}")
























