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
path = "../../analysis/plots/lineplots_BCRmetrics/"

#run_type = "dry"
run_type = "wet"

#organ = "ILE"
organ = "SPL"

palette2 = {"DR_lifelong":"#D55E00", "AL_lifelong":"#0072B2", "AL_DR16M":"#009E73", "AL_DR20M":"#F0E442"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

# {{{
if organ == "SPL":
    df = pd.read_csv("../../analysis/results_tables/allBCRmetrics_tumpathoscore_20220415.csv", sep = ";")
else:
    df = pd.read_csv("../../analysis/results_tables/allBCRmetrics_tumpathoscore_ILE_20220415.csv", sep = ";")

df["Diet"] = df["Diet"].replace({0:"DR_lifelong", 3:"AL_lifelong", 1:"AL_DR16M", 2:"AL_DR20M"})
df = df[df["Age"] != 16]
# }}}



# {{{
lis_test = ['p20', 'mu_freq_seq_r','mu_freq_seq_s','RDI_uniqueness','Richness', 'Shannon', 'Simpson', "cdr3_mu", "cdr3_sigma"]

dic_test = {'p20':"P20 frequency", 'mu_freq_seq_r':"Mu freq [Non-synonymous]",'mu_freq_seq_s':"Mu freq [Synonymous]",
            'RDI_uniqueness':"RDI uniqueness",'Richness':"Richness (q = 0)", 'Shannon':"Shannon (q = 1)", 'Simpson': "Simpson (q = 2)",
           "cdr3_mu": "CDR3 mean mu", "cdr3_sigma": "CDR3 mean sigma"}

for e in lis_test:

    fig, ax = plt.subplots(figsize = (4,2.5))
    sns.lineplot(data = df, x = "Age", y = e, hue = "Diet", ax=ax, palette=palette2,
                linewidth = 1.5, hue_order = ["DR_lifelong", "AL_DR16M", "AL_DR20M", "AL_lifelong"], marker = "o")


    plt.axvline(x = 16, alpha = 0.3, color = "teal")
    plt.axvline(x = 20, alpha = 0.3, color = "gold")
    ax.tick_params(axis = "x", labelsize=8)
    ax.tick_params(axis = "y", labelsize=8)
    ax.set_xlabel("Age (month)", fontsize = 9)
    ax.set_ylabel(dic_test[e], fontsize = 9)
    plt.ylim(0, None)
    plt.xlim(0, 25)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 9)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(1.5)

    matplotlib.rcParams['pdf.fonttype'] = 42 

    plt.tight_layout()

    if run_type != "dry":

        plt.savefig("{}{}_{}_line_{}.pdf".format(path, organ, e, new_day))

    else:
        plt.show()
    
# }}}

spl = pd.read_csv("../../analysis/results_tables/RDI_biogroups_{}.tsv".format(organ), sep = " ")
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
spl = spl[spl["value"] != 0.0]
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
spl = spl.drop_duplicates()
spl = spl[~spl['col'].isin(['value'])]

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
splsub = spl[spl["Age"] != 16]


fig, ax = plt.subplots(figsize = (4,2.5))
sns.lineplot(data = splsub, x = "Age", y = "value", hue = "Diet", ax=ax, palette=palette2,
    linewidth = 1.5, hue_order = ["DR_lifelong", "AL_DR16M", "AL_DR20M", "AL_lifelong"], marker = "o")

plt.axvline(x = 16, alpha = 0.3, color = "teal")
plt.axvline(x = 20, alpha = 0.3, color = "gold")
ax.tick_params(axis = "x", labelsize=8)
ax.tick_params(axis = "y", labelsize=8)
ax.set_xlabel("Age (month)", fontsize = 9)
ax.set_ylabel("Repertoire dissimilarity index", fontsize = 9)
plt.ylim(0, None)
plt.xlim(0, 25)

leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for line in leg.get_lines():
    line.set_linewidth(1.5)

matplotlib.rcParams['pdf.fonttype'] = 42 

plt.tight_layout()

if run_type != "dry":

    plt.savefig("{}{}_RDI_line_{}.pdf".format(path, organ, new_day))

else:
    plt.show()
# }}}
df.columns = ['Unnamed: 0', 'Animal_ID', 'biogroup', 'Age', 'Diet', 'p20',
       'mu_freq_seq_r', 'mu_freq_seq_s', 'cdr3_mu', 'cdr3_sigma',
       'Post_antigenic', 'Clone_antigen_exposed', 'Naive', 'RDI_uniqueness',
       'IGHV1', 'IGHV10', 'IGHV11', 'IGHV12', 'IGHV13', 'IGHV14', 'IGHV15',
       'IGHV16', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7', 'IGHV8',
       'IGHV9', 'IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'relAb_IGA', 'relAb_IGD',
       'relAb_IGE', 'relAb_IGG', 'relAb_IGM', 'RDI_uniqueness_IGA',
       'RDI_uniqueness_IGD', 'RDI_uniqueness_IGE', 'RDI_uniqueness_IGM',
       'RDI_uniqueness_IGG', 'p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE',
       'p20_IGG', 'mu_freq_seq_r_IGA', 'mu_freq_seq_s_IGA',
       'mu_freq_seq_r_IGM', 'mu_freq_seq_s_IGM', 'mu_freq_seq_r_IGG',
       'mu_freq_seq_s_IGG', 'mu_freq_seq_r_IGD', 'mu_freq_seq_s_IGD',
       'mu_freq_seq_r_IGE', 'mu_freq_seq_s_IGE', 'cdr3_mu_IGA',
       'cdr3_sigma_IGA', 'cdr3_mu_IGE', 'cdr3_sigma_IGE', 'cdr3_mu_IGD',
       'cdr3_sigma_IGD', 'cdr3_mu_IGG', 'cdr3_sigma_IGG', 'cdr3_mu_IGM',
       'cdr3_sigma_IGM', 'presAbsPatho', 'presAbsTumor', 'Morbidity',
       'Richness', 'Shannon', 'Simpson', 'Richness_IGA', 'Richness_IGE',
       'Richness_IGD', 'Richness_IGM', 'Richness_IGG', 'Shannon_IGA',
       'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM', 'Shannon_IGG',
       'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD', 'Simpson_IGM',
       'Simpson_IGG']

# # ISOTYPES


df.columns

# {{{
if organ == "SPL":
    lis_iso = ['p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE','p20_IGG', 'Richness_IGA', 'Richness_IGE',
       'Richness_IGD', 'Richness_IGM', 'Richness_IGG', 'Shannon_IGA',
       'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM', 'Shannon_IGG',
       'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD', 'Simpson_IGM',
       'Simpson_IGG',  'mu_count_seq_r_IGA', 'mu_count_seq_s_IGA',
       'mu_count_seq_r_IGM', 'mu_count_seq_s_IGM', 'mu_count_seq_r_IGG',
       'mu_count_seq_s_IGG', 'mu_count_seq_r_IGD', 'mu_count_seq_s_IGD',
       'mu_count_seq_r_IGE', 'mu_count_seq_s_IGE']
    dic_iso = {'p20_IGA':"[IgA] P20 frequency", 'p20_IGM':"[IgM] P20 frequency", 'p20_IGD':"[IgD] P20 frequency", 
           'p20_IGE':"[IgE] P20 frequency",'p20_IGG':"[IgG] P20 frequency", 'Richness_IGA':"[IgA] Richness (q = 0)", 
           'Richness_IGE':"[IgE] Richness (q = 0)", 'Richness_IGD':"[IgD] Richness (q = 0)", 'Richness_IGM':"[IgM] Richness (q = 0)", 
           'Richness_IGG':"[IgG] Richness (q = 0)", 'Shannon_IGA':"[IgA] Shannon (q = 1)",
           'Shannon_IGE':"[IgE] Shannon (q = 1)", 'Shannon_IGD':"[IgD] Shannon (q = 1)", 'Shannon_IGM':"[IgM] Shannon (q = 1)", 
           'Shannon_IGG':"[IgG] Shannon (q = 1)",'Simpson_IGA':"[IgA] Simpson (q = 1)", 'Simpson_IGE':"[IgE] Simpson (q = 1)", 
           'Simpson_IGD':"[IgD] Simpson (q = 1)", 'Simpson_IGM':"[IgM] Simpson (q = 1)", 'Simpson_IGG':"[IgG] Simpson (q = 1)", 
          'mu_count_seq_r_IGA':"[IgA] Mu freq [Non-synonymous]", 'mu_count_seq_s_IGA':"[IgA] Mu freq [Synonymous]",
       'mu_count_seq_r_IGM':"[IgM] Mu freq [Non-synonymous]", 'mu_count_seq_s_IGM':"[IgM] Mu freq [Synonymous]", 'mu_count_seq_r_IGG':"[IgG] Mu freq [Non-synonymous]",
       'mu_count_seq_s_IGG':"[IgG] Mu freq [Synonymous]", 'mu_count_seq_r_IGD':"[IgD] Mu freq [Non-synonymous]", 'mu_count_seq_s_IGD':"[IgD] Mu freq [Synonymous]",
       'mu_count_seq_r_IGE':"[IgE] Mu freq [Non-synonymous]", 'mu_count_seq_s_IGE':"[IgE] Mu freq [Synonymous]"}
else:

    lis_iso = ['p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE','p20_IGG', 'Richness_IGA', 'Richness_IGE',
       'Richness_IGD', 'Richness_IGM', 'Richness_IGG', 'Shannon_IGA',
       'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM', 'Shannon_IGG',
       'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD', 'Simpson_IGM',
       'Simpson_IGG',  'mu_freq_seq_r_IGA', 'mu_freq_seq_s_IGA',
       'mu_freq_seq_r_IGM', 'mu_freq_seq_s_IGM', 'mu_freq_seq_r_IGG',
       'mu_freq_seq_s_IGG', 'mu_freq_seq_r_IGD', 'mu_freq_seq_s_IGD',
       'mu_freq_seq_r_IGE', 'mu_freq_seq_s_IGE']
    
    dic_iso = {'p20_IGA':"[IgA] P20 frequency", 'p20_IGM':"[IgM] P20 frequency", 'p20_IGD':"[IgD] P20 frequency", 
           'p20_IGE':"[IgE] P20 frequency",'p20_IGG':"[IgG] P20 frequency", 'Richness_IGA':"[IgA] Richness (q = 0)", 
           'Richness_IGE':"[IgE] Richness (q = 0)", 'Richness_IGD':"[IgD] Richness (q = 0)", 'Richness_IGM':"[IgM] Richness (q = 0)", 
           'Richness_IGG':"[IgG] Richness (q = 0)", 'Shannon_IGA':"[IgA] Shannon (q = 1)",
           'Shannon_IGE':"[IgE] Shannon (q = 1)", 'Shannon_IGD':"[IgD] Shannon (q = 1)", 'Shannon_IGM':"[IgM] Shannon (q = 1)", 
           'Shannon_IGG':"[IgG] Shannon (q = 1)",'Simpson_IGA':"[IgA] Simpson (q = 1)", 'Simpson_IGE':"[IgE] Simpson (q = 1)", 
           'Simpson_IGD':"[IgD] Simpson (q = 1)", 'Simpson_IGM':"[IgM] Simpson (q = 1)", 'Simpson_IGG':"[IgG] Simpson (q = 1)", 
          'mu_freq_seq_r_IGA':"[IgA] Mu freq [Non-synonymous]", 'mu_freq_seq_s_IGA':"[IgA] Mu freq [Synonymous]",
       'mu_freq_seq_r_IGM':"[IgM] Mu freq [Non-synonymous]", 'mu_freq_seq_s_IGM':"[IgM] Mu freq [Synonymous]", 'mu_freq_seq_r_IGG':"[IgG] Mu freq [Non-synonymous]",
       'mu_freq_seq_s_IGG':"[IgG] Mu freq [Synonymous]", 'mu_freq_seq_r_IGD':"[IgD] Mu freq [Non-synonymous]", 'mu_freq_seq_s_IGD':"[IgD] Mu freq [Synonymous]",
       'mu_freq_seq_r_IGE':"[IgE] Mu freq [Non-synonymous]", 'mu_freq_seq_s_IGE':"[IgE] Mu freq [Synonymous]"}



for e in lis_iso:
    
    df = df[df["Age"] != 16]

    fig, ax = plt.subplots(figsize = (4,2.5))
    sns.lineplot(data = df, x = "Age", y = e, hue = "Diet", ax=ax, palette=palette2,
                linewidth = 1.5, hue_order = ["DR_lifelong", "AL_DR16M", "AL_DR20M", "AL_lifelong"], marker = "o")


    plt.axvline(x = 16, alpha = 0.3, color = "teal")
    plt.axvline(x = 20, alpha = 0.3, color = "gold")
    ax.tick_params(axis = "x", labelsize=8)
    ax.tick_params(axis = "y", labelsize=8)
    ax.set_xlabel("Age (month)", fontsize = 9)
    ax.set_ylabel(dic_iso[e], fontsize = 9)
    plt.ylim(0, None)
    plt.xlim(0, 25)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 8)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(1.5)

    matplotlib.rcParams['pdf.fonttype'] = 42 

    plt.tight_layout()

    if run_type != "dry":

        plt.savefig("{}{}_{}_line_{}.pdf".format(path, organ, e, new_day))

    else:
        plt.show()
# }}}

# {{{
igA = pd.read_csv("../../analysis/results_tables/RDI_IGA_biogroup_{}_isosum.tsv".format(organ), sep = " ")
igD = pd.read_csv("../../analysis/results_tables/RDI_IGD_biogroup_{}_isosum.tsv".format(organ), sep = " ")
igE = pd.read_csv("../../analysis/results_tables/RDI_IGE_biogroup_{}_isosum.tsv".format(organ), sep = " ")
igG = pd.read_csv("../../analysis/results_tables/RDI_IGG_biogroup_{}_isosum.tsv".format(organ), sep = " ")
igM = pd.read_csv("../../analysis/results_tables/RDI_IGM_biogroup_{}_isosum.tsv".format(organ), sep = " ")

listup = [(igA, "IgA"), (igD, "IgD"), (igE, "IgE"), (igG, "IgG"), (igM, "IgM")]
# }}}

listup2 = []
for d, i in listup:
    d = pd.DataFrame(d)
    # We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
    d = d[d["value"] != 0.0]
    # We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
    d = d.drop_duplicates()
    d = d[~d['col'].isin(['value'])]
    listup2.append((d, i))

listup3 = []
for d, i in listup2:
    # Getting the info we need into our dataframe
    d = pd.DataFrame(d)
    diets = []
    for e in d["row"].to_list():
        if e[2:] == "R_lifelong":
            diets.append("DR_lifelong")
        elif e[2:] == "L_lifelong":
            diets.append("AL_lifelong")
        else:
            diets.append(e[2:])

    ages = []
    for e in d["row"].to_list():
        ages.append(e[:2])
    ages = [w.replace('5D', "5") for w in ages]
    ages = [w.replace('5A', "5") for w in ages]

    ages = [ int(x) for x in ages ]

    # Add to dataframe
    d["Age"] = ages
    d["Diet"] = diets
    listup3.append((d, i))

for d, i in listup3:
    spl = pd.DataFrame(d)

    splsub = spl[spl["Age"] != 16]


    fig, ax = plt.subplots(figsize = (4,2.5))
    sns.lineplot(data = splsub, x = "Age", y = "value", hue = "Diet", ax=ax, palette=palette2,
        linewidth = 1.5, hue_order = ["DR_lifelong", "AL_DR16M", "AL_DR20M", "AL_lifelong"], marker = "o")

    plt.axvline(x = 16, alpha = 0.3, color = "teal")
    plt.axvline(x = 20, alpha = 0.3, color = "gold")
    ax.tick_params(axis = "x", labelsize=8)
    ax.tick_params(axis = "y", labelsize=8)
    ax.set_xlabel("Age (month)", fontsize = 9)
    ax.set_ylabel("[{}] Repertoire dissimilarity index".format(i), fontsize = 9)
    plt.ylim(0, None)
    plt.xlim(0, 25)

    leg = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 8)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in leg.get_lines():
        line.set_linewidth(1.5)

    matplotlib.rcParams['pdf.fonttype'] = 42 

    plt.tight_layout()

    if run_type != "dry":

        plt.savefig("{}{}_RDI_{}_line_{}.pdf".format(path, organ, i, new_day))

    else:
        plt.show()




