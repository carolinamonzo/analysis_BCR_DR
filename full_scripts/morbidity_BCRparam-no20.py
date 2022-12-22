# ---
# jupyter:
#   jupytext:
#     cell_markers: '{{{,}}}'
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: MPI_cmc
#     language: python
#     name: mpi_cmc
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
path = "../../analysis/plots/path_BCR/"

#run_type = "dry"
run_type = "wet"

#organ = "ILE"
organ = "SPL"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", "AL_DR16M":"teal", "AL_DR20M":"gold"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

import scipy
print(scipy.__version__)

# {{{
df_patho = pd.read_csv("../../analysis/results_tables/pathology_allBCRmetrics_20220404.csv", sep = ";")
df_tum = pd.read_csv("../../analysis/results_tables/tumor_allBCRmetrics_20220404.csv", sep = ";")

old_patho = pd.read_csv("../../analysis/results_tables/{}_PathologyMergeAlphaP20SHMisoabundCS.csv".format(organ), sep = ";")
old_tum = pd.read_csv("../../analysis/results_tables/{}_tumourMergeAlphaP20SHMisoabundCS.csv".format(organ), sep = ";")
# }}}

# {{{
patho_tot = pd.DataFrame(old_patho.loc[:, ['Animal_ID',
       'Meager,_hardly_any_WAT', 'Crooked_back/_bad_habitus',
       'Enlarged_spleen', 'Small_spleen', 'Granulated_spleen',
       'WAT_brown', 'Pancreas_knotty',
       'Pancreas_pale/_discoloured',
       'Gastrointestinal_bleedings', 'Bloated_stomach', 'Fluid_filled_stomach',
       'Bloated_SI_and_colon', 'Hepatomegalia',
       'Liver_pale',
       'Liver_appears_granulated', 'Liver_yellow_/_Jaundice',
       'Uterine_cyst_left_horn', 'Uterine_cyst_right_horn',
       'Bleedings_into_uterus', 'Cyst_left_ovary',
       'Cyst_right_ovary', 'Enlarged_heart', 'Enlarged_kidneys',
       'Fully-filled_gall_bladder', 'Internal_bleeding_into_abdominal_cavity',
       'Internal_bleeding_into_pleural_cavity', 'WAT_cyst',
       'Kidney_spotty_or_granulated', 'Hind_leg_paralysis',
       'Agglutinated_vagina', 'Gasping_for_air/_agonal_respiration',
       'Cold_to_the_touch', 'On_the_verge_of_dying', 'Critical_weight_loss',
       'enlarged_adrenal_gland', 'Bloody_poo', 'Black_gall_bladder',
       'Indigestion', 'Other_cysts',
       'WAT_in_pleural_cavity', 'lymph_vessels_close_to_spine_thickened',
       'Ill-kept,_rumpled_furr', 'Uterus_thickened_and_hard',
       'Discoloured_kidneys', 'Kidneys_pale', 'Free_fluid_in_abdominal_cavity',
       'Enlarged_bladder', 'Open_wounds_on_tail', 'Enlarged_lungs',
       'Lungs_pale', 'Skull_bone_dissolved_due_to_eye_tumour',
       'Organs_dislocated/shifted_due_to_tumour']].sum(axis = 0))

tum_tot = pd.DataFrame(df_tum.loc[:, ['Animal_ID', 'Eye', 'Connective_tissue', 'Liver', 'Kidney', 'Ovary',
       'Uterine', 'Lymph_node', 'Lung', 'Colorectal', 'Brain', 'Heart']].sum(axis = 0))
# }}}

tum_tot[tum_tot[0] != 0].index

patho_tot[patho_tot[0] != 0].index

# {{{
affec_patho = old_patho.loc[:, ['Animal_ID', 'Crooked_back/_bad_habitus', 'Enlarged_spleen',
       'WAT_brown', 'Pancreas_knotty', 'Uterine_cyst_left_horn',
       'Uterine_cyst_right_horn', 'Enlarged_kidneys',
       'enlarged_adrenal_gland']].drop_duplicates()

affec_tum = df_tum.loc[:, ['Animal_ID', 'Connective_tissue', 'Liver', 'Lung']].drop_duplicates()
# }}}

affec_patho["total"] = affec_patho.sum(axis = 1)
affec_patho["Score_patho"] = affec_patho["total"]

affec_tum["total"] = affec_tum.sum(axis = 1)
affec_tum["Score_tum"] = affec_tum["total"].replace(3, 2)

totest = pd.merge(affec_patho.loc[:, ["Animal_ID", "Score_patho"]], affec_tum.loc[:, ["Animal_ID", "Score_tum"]], on = "Animal_ID", how = "inner")
totest["Morbidity"] = totest["Score_patho"] + totest["Score_tum"]

df = df_patho.loc[:, ["Animal_ID", 'biogroup',
       'q', 'd', 'Age', 'Diet', 'p20', 'mu_freq_seq_r', 'mu_freq_seq_s',
        'cdr3_mu', 'cdr3_sigma', 'IgM-IgD-', 'IgM+IgD+SHM+',
       'IgM+IgD+SHM-', 'RDI_uniqueness', 'IGHV1', 'IGHV10', 'IGHV11', 'IGHV12',
       'IGHV13', 'IGHV14', 'IGHV15', 'IGHV16', 'IGHV2', 'IGHV3', 'IGHV4',
       'IGHV5', 'IGHV6', 'IGHV7', 'IGHV8', 'IGHV9', 'IGHJ1', 'IGHJ2', 'IGHJ3',
       'IGHJ4', 'relAb_IGA', 'relAb_IGD', 'relAb_IGE', 'relAb_IGG',
       'relAb_IGM', 'alpha_IGA', 'alpha_IGD', 'alpha_IGE', 'alpha_IGM',
       'alpha_IGG', 'RDI_uniqueness_IGA', 'RDI_uniqueness_IGD',
       'RDI_uniqueness_IGE', 'RDI_uniqueness_IGM', 'RDI_uniqueness_IGG',
       'p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE', 'p20_IGG',
       'mu_count_seq_r_IGA', 'mu_count_seq_s_IGA', 
       'mu_count_seq_r_IGM', 'mu_count_seq_s_IGM', 
       'mu_count_seq_r_IGG', 'mu_count_seq_s_IGG', 
       'mu_count_seq_r_IGD', 'mu_count_seq_s_IGD', 
       'mu_count_seq_r_IGE', 'mu_count_seq_s_IGE', 
       'cdr3_mu_IGA', 'cdr3_sigma_IGA', 'cdr3_mu_IGE', 'cdr3_sigma_IGE',
       'cdr3_mu_IGD', 'cdr3_sigma_IGD', 'cdr3_mu_IGG', 'cdr3_sigma_IGG',
       'cdr3_mu_IGM', 'cdr3_sigma_IGM', 'presAbsPatho']]
df = pd.merge(df, df_tum.loc[:, ["Animal_ID", "presAbsTumor"]], on = "Animal_ID", how = "inner")
df = pd.merge(df, totest.loc[:, ["Animal_ID", "Morbidity"]], on = "Animal_ID", how = "inner")

# {{{
alpha = pd.DataFrame()
alpha["Animal_ID"] = df["Animal_ID"].drop_duplicates()
df["q"] = pd.to_numeric(df["q"])

alpha["Richness"] = list(df[df["q"] == 0].drop_duplicates()["d"].dropna())
alpha["Shannon"] = list(df[df["q"] == 1].drop_duplicates()["d"].dropna())
alpha["Simpson"] = list(df[df["q"] == 2].drop_duplicates()["d"].dropna())

alpha["Richness_IGA"] = list(df[df["q"] == 0].drop_duplicates()["alpha_IGA"].dropna())
alpha["Richness_IGE"] = list(df[df["q"] == 0].drop_duplicates()["alpha_IGE"].dropna())
alpha["Richness_IGD"] = list(df[df["q"] == 0].drop_duplicates()["alpha_IGD"].dropna())
alpha["Richness_IGM"] = list(df[df["q"] == 0].drop_duplicates()["alpha_IGM"].dropna())
alpha["Richness_IGG"] = list(df[df["q"] == 0].drop_duplicates()["alpha_IGG"].dropna())

alpha["Shannon_IGA"] = list(df[df["q"] == 1].drop_duplicates()["alpha_IGA"].dropna())
alpha["Shannon_IGE"] = list(df[df["q"] == 1].drop_duplicates()["alpha_IGE"].dropna())
alpha["Shannon_IGD"] = list(df[df["q"] == 1].drop_duplicates()["alpha_IGD"].dropna())
alpha["Shannon_IGM"] = list(df[df["q"] == 1].drop_duplicates()["alpha_IGM"].dropna())
alpha["Shannon_IGG"] = list(df[df["q"] == 1].drop_duplicates()["alpha_IGG"].dropna())

alpha["Simpson_IGA"] = list(df[df["q"] == 2].drop_duplicates()["alpha_IGA"].dropna())
alpha["Simpson_IGE"] = list(df[df["q"] == 2].drop_duplicates()["alpha_IGE"].dropna())
alpha["Simpson_IGD"] = list(df[df["q"] == 2].drop_duplicates()["alpha_IGD"].dropna())
alpha["Simpson_IGM"] = list(df[df["q"] == 2].drop_duplicates()["alpha_IGM"].dropna())
alpha["Simpson_IGG"] = list(df[df["q"] == 2].drop_duplicates()["alpha_IGG"].dropna())
# }}}

df = pd.merge(df.drop(columns = ["q", "d", "alpha_IGA", "alpha_IGE", "alpha_IGM", "alpha_IGD", "alpha_IGG"]).drop_duplicates(), alpha.drop_duplicates(), on = "Animal_ID", how = "inner")
df["Diet"] = df["Diet"].replace({"DR lifelong": 0, "AL lifelong": 3, "AL_DR16M":1, "AL_DR20M":2})

# {{{
df.columns = ['Animal_ID', 'biogroup', 'Age', 'Diet', 'p20', 'mu_freq_seq_r',
       'mu_freq_seq_s', 'cdr3_mu', 'cdr3_sigma', 'Post_antigenic',
       'Clone_antigen_exposed', 'Naive', 'RDI_uniqueness', 'IGHV1', 'IGHV10',
       'IGHV11', 'IGHV12', 'IGHV13', 'IGHV14', 'IGHV15', 'IGHV16', 'IGHV2',
       'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7', 'IGHV8', 'IGHV9', 'IGHJ1',
       'IGHJ2', 'IGHJ3', 'IGHJ4', 'relAb_IGA', 'relAb_IGD', 'relAb_IGE',
       'relAb_IGG', 'relAb_IGM', 'RDI_uniqueness_IGA', 'RDI_uniqueness_IGD',
       'RDI_uniqueness_IGE', 'RDI_uniqueness_IGM', 'RDI_uniqueness_IGG',
       'p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE', 'p20_IGG',
       'mu_count_seq_r_IGA', 'mu_count_seq_s_IGA', 
       'mu_count_seq_r_IGM', 'mu_count_seq_s_IGM', 
       'mu_count_seq_r_IGG', 'mu_count_seq_s_IGG', 
       'mu_count_seq_r_IGD', 'mu_count_seq_s_IGD', 
       'mu_count_seq_r_IGE', 'mu_count_seq_s_IGE', 
       'cdr3_mu_IGA', 'cdr3_sigma_IGA', 'cdr3_mu_IGE', 'cdr3_sigma_IGE',
       'cdr3_mu_IGD', 'cdr3_sigma_IGD', 'cdr3_mu_IGG', 'cdr3_sigma_IGG',
       'cdr3_mu_IGM', 'cdr3_sigma_IGM', 'presAbsPatho', 'presAbsTumor',
       'Morbidity', 'Richness', 'Shannon', 'Simpson', 'Richness_IGA',
       'Richness_IGE', 'Richness_IGD', 'Richness_IGM', 'Richness_IGG',
       'Shannon_IGA', 'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM',
       'Shannon_IGG', 'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD',
       'Simpson_IGM', 'Simpson_IGG']

# Since we cant have symbols, had to change these class switched states"
#'IgM-IgD-','IgM+IgD+SHM+', 'IgM+IgD+SHM-'
#'Post_antigenic','Clone_antigen_exposed', 'Naive'
# }}}

df.to_csv("../../analysis/results_tables/allBCRmetrics_tumpathoscore_{}.csv".format(new_day), sep = ";")

# {{{
p_value = []
rsquared = []
measure = []
#for e in ['p20', 'mu_freq', 'mut_status',
#       'relAb_IGA', 'relAb_IGD', 'relAb_IGE', 'relAb_IGG', 'relAb_IGM',
#       'Naive_perc', 'CS_perc', 'Richness', 'Shannon', 'Simpson']:

#for e in ["Age", "Diet",'p20', 'mu_freq_seq_r',
#       'mu_freq_seq_s', 'mut_status', 'cdr3_mu', 'cdr3_sigma', 'Post_antigenic',
#       'Clone_antigen_exposed', 'Naive', 'RDI_uniqueness', 'IGHV1', 'IGHV10',
#       'IGHV11', 'IGHV12', 'IGHV13', 'IGHV14', 'IGHV15', 'IGHV16', 'IGHV2',
#       'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7', 'IGHV8', 'IGHV9', 'IGHJ1',
#       'IGHJ2', 'IGHJ3', 'IGHJ4', 'relAb_IGA', 'relAb_IGD', 'relAb_IGE',
#       'relAb_IGG', 'relAb_IGM', 'RDI_uniqueness_IGA', 'RDI_uniqueness_IGD',
#       'RDI_uniqueness_IGE', 'RDI_uniqueness_IGM', 'RDI_uniqueness_IGG',
#       'p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE', 'p20_IGG',
#       'mu_count_seq_r_IGA', 'mu_count_seq_s_IGA', 'mut_status_IGA',
#       'mu_count_seq_r_IGM', 'mu_count_seq_s_IGM', 'mut_status_IGM',
#       'mu_count_seq_r_IGG', 'mu_count_seq_s_IGG', 'mut_status_IGG',
#       'mu_count_seq_r_IGD', 'mu_count_seq_s_IGD', 'mut_status_IGD',
#       'mu_count_seq_r_IGE', 'mu_count_seq_s_IGE', 'mut_status_IGE',
#       'cdr3_mu_IGA', 'cdr3_sigma_IGA', 'cdr3_mu_IGE', 'cdr3_sigma_IGE',
#       'cdr3_mu_IGD', 'cdr3_sigma_IGD', 'cdr3_mu_IGG', 'cdr3_sigma_IGG',
#       'cdr3_mu_IGM', 'cdr3_sigma_IGM',
#       'Richness', 'Shannon', 'Simpson', 'Richness_IGA',
#       'Richness_IGE', 'Richness_IGD', 'Richness_IGM', 'Richness_IGG',
#       'Shannon_IGA', 'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM',
#       'Shannon_IGG', 'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD',
#       'Simpson_IGM', 'Simpson_IGG']:

ls_test = ["Age", "Diet",'p20', 'mu_freq_seq_r',
       'mu_freq_seq_s','cdr3_mu', 'cdr3_sigma', 'Post_antigenic',
       'Clone_antigen_exposed', 'Naive', 'RDI_uniqueness',
       'Richness', 'Shannon', 'Simpson']

for e in ls_test:


    results = ols('Morbidity ~ {}'.format(e), data = df).fit()
    p_value.append(results.pvalues.tolist()[1])
    rsquared.append(results.rsquared)
    measure.append(e)
# }}}


toplot = pd.DataFrame({"p_value_linreg":p_value, "rsquared":rsquared, "Measure":measure})
toplot = toplot.sort_values(by = "rsquared", ascending = False)
#toplot.set_index("measure", inplace = True)

# {{{
p_value = []
corre = []
measure = []
for e in ls_test:

    rho, p_val = stats.spearmanr(df["Morbidity"].to_list(), df[e].to_list())
    p_value.append(p_val)
    corre.append(rho)
    measure.append(e)
    
df_spearman = pd.DataFrame({"Measure":measure, "Correlation":corre, "p_value_spearman":p_value})
df_spearman.sort_values(by = "p_value_spearman", ascending = True, inplace = True)
# }}}

dfmer = pd.merge(df_spearman, toplot, on = "Measure", how = "inner").dropna().sort_values(by = "rsquared", ascending = False)

# Enable R
# %load_ext rpy2.ipython

# {{{ magic_args="-i dfmer -w 5 -h 5 --units in -r 200 -o dfmer" language="R"
# dfmer$padjustSpearman <- p.adjust(dfmer$p_value_spearman, "bonferroni")
# dfmer$padjustLinreg <- p.adjust(dfmer$p_value_linreg, "bonferroni")
# }}}

dfmer.sort_values("p_value_spearman")

# {{{
# Min 0,05 from spearman
#dfmer = dfmer[dfmer["p_value_spearman"] < 0.05]
# Make positive and negative rsquared according to direction of spearman correlation
l = []
for e in dfmer["Correlation"]:
    if e > 0:
        l.append(1)
    elif e < 0:
        l.append(-1)
    else:
        print("Error")
        
dfmer["Direction"] = l
dfmer["rsquared"] = dfmer["rsquared"] * dfmer["Direction"]

l2 = []
for e in dfmer["p_value_spearman"]:
    if e < 0.05:
        l2.append("magenta")
    elif e >= 0.05:
        l2.append("lightgrey")
    else:
        print("error")
        
dfmer["Color"] = l2
# }}}

# {{{
# Blue are significant negative spearman correlation
# Orange are significant positive spearman correlation
fig, ax = plt.subplots(figsize = (6,6))
sns.barplot(data = dfmer, y = "rsquared", x = "Measure", palette = dfmer["Color"].to_list())

ax.tick_params(axis = "x", labelsize = 16, rotation = 90)
ax.tick_params(axis = "y", labelsize = 16)
ax.set_ylabel("Variance explained (r^2)",  fontsize = 18)
ax.set_xlabel(" ")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig(f"{path}/SPL_BCR_VarExplained_{new_day}.pdf")
else:
    plt.show()
# }}}

# ## Now for the isotypes

# {{{
p_value = []
rsquared = []
measure = []
#for e in ['p20', 'mu_freq', 'mut_status',
#       'relAb_IGA', 'relAb_IGD', 'relAb_IGE', 'relAb_IGG', 'relAb_IGM',
#       'Naive_perc', 'CS_perc', 'Richness', 'Shannon', 'Simpson']:

#for e in ["Age", "Diet",'p20', 'mu_freq_seq_r',
#       'mu_freq_seq_s', 'mut_status', 'cdr3_mu', 'cdr3_sigma', 'Post_antigenic',
#       'Clone_antigen_exposed', 'Naive', 'RDI_uniqueness', 'IGHV1', 'IGHV10',
#       'IGHV11', 'IGHV12', 'IGHV13', 'IGHV14', 'IGHV15', 'IGHV16', 'IGHV2',
#       'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7', 'IGHV8', 'IGHV9', 'IGHJ1',
#       'IGHJ2', 'IGHJ3', 'IGHJ4', 'relAb_IGA', 'relAb_IGD', 'relAb_IGE',
#       'relAb_IGG', 'relAb_IGM', 'RDI_uniqueness_IGA', 'RDI_uniqueness_IGD',
#       'RDI_uniqueness_IGE', 'RDI_uniqueness_IGM', 'RDI_uniqueness_IGG',
#       'p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE', 'p20_IGG',
#       'mu_count_seq_r_IGA', 'mu_count_seq_s_IGA', 'mut_status_IGA',
#       'mu_count_seq_r_IGM', 'mu_count_seq_s_IGM', 'mut_status_IGM',
#       'mu_count_seq_r_IGG', 'mu_count_seq_s_IGG', 'mut_status_IGG',
#       'mu_count_seq_r_IGD', 'mu_count_seq_s_IGD', 'mut_status_IGD',
#       'mu_count_seq_r_IGE', 'mu_count_seq_s_IGE', 'mut_status_IGE',
#       'cdr3_mu_IGA', 'cdr3_sigma_IGA', 'cdr3_mu_IGE', 'cdr3_sigma_IGE',
#       'cdr3_mu_IGD', 'cdr3_sigma_IGD', 'cdr3_mu_IGG', 'cdr3_sigma_IGG',
#       'cdr3_mu_IGM', 'cdr3_sigma_IGM',
#       'Richness', 'Shannon', 'Simpson', 'Richness_IGA',
#       'Richness_IGE', 'Richness_IGD', 'Richness_IGM', 'Richness_IGG',
#       'Shannon_IGA', 'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM',
#       'Shannon_IGG', 'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD',
#       'Simpson_IGM', 'Simpson_IGG']:

ls_test = ['relAb_IGA', 'relAb_IGD', 'relAb_IGE',
       'relAb_IGG', 'relAb_IGM', 'RDI_uniqueness_IGA', 'RDI_uniqueness_IGD',
       'RDI_uniqueness_IGE', 'RDI_uniqueness_IGM', 'RDI_uniqueness_IGG',
       'p20_IGA', 'p20_IGM', 'p20_IGD', 'p20_IGE', 'p20_IGG',
       'mu_count_seq_r_IGA', 'mu_count_seq_s_IGA', 
       'mu_count_seq_r_IGM', 'mu_count_seq_s_IGM', 
       'mu_count_seq_r_IGG', 'mu_count_seq_s_IGG', 
       'mu_count_seq_r_IGD', 'mu_count_seq_s_IGD', 
       'mu_count_seq_r_IGE', 'mu_count_seq_s_IGE', 
       'cdr3_mu_IGA', 'cdr3_sigma_IGA', 'cdr3_mu_IGE', 'cdr3_sigma_IGE',
       'cdr3_mu_IGD', 'cdr3_sigma_IGD', 'cdr3_mu_IGG', 'cdr3_sigma_IGG',
       'cdr3_mu_IGM', 'cdr3_sigma_IGM','Richness_IGA',
       'Richness_IGE', 'Richness_IGD', 'Richness_IGM', 'Richness_IGG',
       'Shannon_IGA', 'Shannon_IGE', 'Shannon_IGD', 'Shannon_IGM',
       'Shannon_IGG', 'Simpson_IGA', 'Simpson_IGE', 'Simpson_IGD',
       'Simpson_IGM', 'Simpson_IGG']

for e in ls_test:
    results = ols('Morbidity ~ {}'.format(e), data = df).fit()
    p_value.append(results.pvalues.tolist()[1])
    rsquared.append(results.rsquared)
    measure.append(e)
# }}}

# {{{
toplot = pd.DataFrame({"p_value_linreg":p_value, "rsquared":rsquared, "Measure":measure})
toplot = toplot.sort_values(by = "rsquared", ascending = False)

p_value = []
corre = []
measure = []
for e in ls_test:

    rho, p_val = stats.spearmanr(df["Morbidity"].to_list(), df[e].to_list())
    p_value.append(p_val)
    corre.append(rho)
    measure.append(e)
    
df_spearman = pd.DataFrame({"Measure":measure, "Correlation":corre, "p_value_spearman":p_value})
df_spearman.sort_values(by = "p_value_spearman", ascending = True, inplace = True)

dfmer = pd.merge(df_spearman, toplot, on = "Measure", how = "inner").dropna().sort_values(by = "rsquared", ascending = False)
# }}}

# {{{
# Make positive and negative rsquared according to direction of spearman correlation
l = []
for e in dfmer["Correlation"]:
    if e > 0:
        l.append(1)
    elif e < 0:
        l.append(-1)
    else:
        print("Error")
        
dfmer["Direction"] = l
dfmer["rsquared"] = dfmer["rsquared"] * dfmer["Direction"]

l2 = []
for e in dfmer["p_value_spearman"]:
    if e < 0.05:
        l2.append("magenta")
    elif e >= 0.05:
        l2.append("lightgrey")
    else:
        print("error")
        
dfmer["Color"] = l2
# }}}

# Enable R
# %load_ext rpy2.ipython

# {{{ magic_args="-i dfmer -w 5 -h 5 --units in -r 200 -o dfmer" language="R"
# dfmer$padjustSpearman <- p.adjust(dfmer$p_value_spearman, "bonferroni")
# dfmer$padjustLinreg <- p.adjust(dfmer$p_value_linreg, "bonferroni")
# }}}

dfmer.sort_values(by = "p_value_spearman")









# {{{
# Blue are significant negative spearman correlation
# Orange are significant positive spearman correlation
fig, ax = plt.subplots(figsize = (14,6))
sns.barplot(data = dfmer, y = "rsquared", x = "Measure", palette = dfmer["Color"].to_list())

ax.tick_params(axis = "x", labelsize = 16, rotation = 90)
ax.tick_params(axis = "y", labelsize = 16)
ax.set_ylabel("Variance explained (r^2)",  fontsize = 18)
ax.set_xlabel(" ")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig(f"{path}/SPL_isotypesBCR_VarExplained_{new_day}.pdf")
else:
    plt.show()
# }}}
# ## Repeat using only AL new


# {{{
p_value = []
rsquared = []
measure = []

df = df[df["Diet"] == 3]

ls_test = ["Age", "Diet",'p20', 'mu_freq_seq_r',
       'mu_freq_seq_s', 'cdr3_mu', 'cdr3_sigma', 'Post_antigenic',
       'Clone_antigen_exposed', 'Naive', 'RDI_uniqueness',
       'Richness', 'Shannon', 'Simpson']

for e in ls_test:


    results = ols('Morbidity ~ {}'.format(e), data = df).fit()
    p_value.append(results.pvalues.tolist()[1])
    rsquared.append(results.rsquared)
    measure.append(e)
# }}}

# {{{
toplot = pd.DataFrame({"p_value_linreg":p_value, "rsquared":rsquared, "Measure":measure})
toplot = toplot.sort_values(by = "rsquared", ascending = False)
#toplot.set_index("measure", inplace = True)

p_value = []
corre = []
measure = []
for e in ls_test:

    rho, p_val = stats.spearmanr(df["Morbidity"].to_list(), df[e].to_list())
    p_value.append(p_val)
    corre.append(rho)
    measure.append(e)
    
df_spearman = pd.DataFrame({"Measure":measure, "Correlation":corre, "p_value_spearman":p_value})
df_spearman.sort_values(by = "p_value_spearman", ascending = True, inplace = True)

dfmer = pd.merge(df_spearman, toplot, on = "Measure", how = "inner").dropna().sort_values(by = "rsquared", ascending = False)
# }}}

# {{{
# Min 0,05 from spearman
#dfmer = dfmer[dfmer["p_value_spearman"] < 0.05]
# Make positive and negative rsquared according to direction of spearman correlation
l = []
for e in dfmer["Correlation"]:
    if e > 0:
        l.append(1)
    elif e < 0:
        l.append(-1)
    else:
        print("Error")
        
dfmer["Direction"] = l
dfmer["rsquared"] = dfmer["rsquared"] * dfmer["Direction"]

l2 = []
for e in dfmer["p_value_spearman"]:
    if e < 0.05:
        l2.append("magenta")
    elif e >= 0.05:
        l2.append("lightgrey")
    else:
        print("error")
        
dfmer["Color"] = l2
# }}}

# {{{
# Blue are significant negative spearman correlation
# Orange are significant positive spearman correlation
fig, ax = plt.subplots(figsize = (6,6))
sns.barplot(data = dfmer, y = "rsquared", x = "Measure", palette = dfmer["Color"].to_list())

ax.tick_params(axis = "x", labelsize = 16, rotation = 90)
ax.tick_params(axis = "y", labelsize = 16)
ax.set_ylabel("Variance explained (r^2)",  fontsize = 18)
ax.set_xlabel(" ")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type != "dry":
    plt.savefig(f"{path}/SPL_ALpathoBCR_VarExplained_{new_day}.pdf")
else:
    plt.show()
# }}}



## REPEAT USING ONLY AL
p_value = []
rsquared = []
measure = []
for e in ['p20', 'mu_freq', 
       'relAb_IGA', 'relAb_IGD', 'relAb_IGE', 'relAb_IGG', 'relAb_IGM',
       'Naive_perc', 'CS_perc', 'Richness', 'Shannon', 'Simpson']:

    results = ols('Morbidity ~ {}'.format(e), data = totest[totest["Diet"] == "AL lifelong"]).fit()
    p_value.append(results.pvalues.tolist()[1])
    rsquared.append(results.rsquared)
    measure.append(e)

toplot = pd.DataFrame({"p_value":p_value, "rsquared":rsquared, "measure":measure})
toplot = toplot.sort_values(by = "rsquared", ascending = False)
toplot.set_index("measure", inplace = True)

# {{{
p_value = []
corre = []
measure = []
for e in ['p20', 'mu_freq', 
       'relAb_IGA', 'relAb_IGD', 'relAb_IGE', 'relAb_IGG', 'relAb_IGM',
       'Naive_perc', 'CS_perc', 'Richness', 'Shannon', 'Simpson']:

    rho, p_val = stats.spearmanr(totest[totest["Diet"] == "AL lifelong"]["Morbidity"].to_list(), totest[totest["Diet"] == "AL lifelong"][e].to_list())
    p_value.append(p_val)
    corre.append(rho)
    measure.append(e)
    
df_spearman = pd.DataFrame({"Measure":measure, "Correlation":corre, "p_value":p_value})
df_spearman.sort_values(by = "p_value", ascending = True, inplace = True)
# }}}

df_spearman

# {{{
# Blue are significant negative spearman correlation
# Orange are significant positive spearman correlation
fig, ax = plt.subplots(figsize = (8,6))
sns.barplot(data = toplot, y = "rsquared", x = toplot.index, palette = ["deepskyblue", 
        "deepskyblue", "lightgrey", "deepskyblue", "deepskyblue", "lightgrey", "lightgrey", 
        "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey"])

ax.tick_params(axis = "x", labelsize = 16, rotation = 90)
ax.tick_params(axis = "y", labelsize = 16)
ax.set_ylabel("Variance explained (r^2)",  fontsize = 18)
ax.set_xlabel(" ")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
if run_type == "dry":
    plt.savefig(f"{path}/SPL_ALpathoBCR_VarExplained_{new_day}.pdf")
else:
    plt.show()
# }}}

toplot


