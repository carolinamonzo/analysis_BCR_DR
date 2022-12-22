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
import scipy.stats as stats
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm
from IPython.display import display
import scikit_posthocs as sp

new_day = datetime.datetime.now().strftime("%Y%m%d")
path_ori = "../../analysis/plots/"
path = "../../analysis/plots/alpha_phenotype/"

run_type = "dry"
#run_type = "wet"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}

iso_palette = {"IGA":"tab:blue", "IGD":"tab:orange", "IGE":"tab:green", "IGG":"tab:red", "IGM":"tab:brown"}
# }}}

# ### Getting and formatting the F2 tumour loads

## Since this should be systemic comparisons, we are doing spleen vs el_diablo_F2
metadata = pd.read_csv("../../metadata/F2-cross-sectional-elDiablo.csv", sep = ";", index_col = "Animal_ID", nrows=288).fillna(0)
metadata = metadata.replace(1000, 10)

# {{{
# Now lets modify the string in the column
tomod = metadata["Cause_of_death"].tolist()
mod = []
for e in tomod:
    if "5M" in e:
        mod.append(5)
    elif "12M" in e:
        mod.append(12)
    elif "16M" in e:
        mod.append(16)
    elif "20M" in e:
        mod.append(20)
    elif "24M" in e:
        mod.append(24)
    elif "28M" in e:
        mod.append(28)
    else:
        print("Error, no months in list")

metadata["Cause_of_death"] = mod
metadata["Cause_of_death"] = metadata["Cause_of_death"].astype(float)
# }}}

dic_diet = {"chronic_AL":"AL lifelong", "chronic_DR":"DR lifelong", "AL_DR52W":"AL_DR12M", 
            "AL_DR69W":"AL_DR16M", "AL_DR87W":"AL_DR20M", "AL_DR104W":"AL_DR24M"}
metadata["Treatment"] = metadata["Cohort"].map(dic_diet)

# {{{
# Total tumour counts
df_col = list(metadata.columns)

tumour_col = []

for e in df_col:
    if "tumour" in e:
        tumour_col.append(e)
    else:
        continue
        

tumour_col = tumour_col[:-2]
tum_df = metadata.loc[:, tumour_col + ["Treatment", "Date_of_pathology", "Date_of_Birth", "Cause_of_death"]]
individuals, _ = tum_df.shape
tum_df = tum_df.replace(100, 6)
# }}}

# {{{
# Count number of tumours per diet group
tum_df.drop(columns = "Cause_of_death")
AL_tum = tum_df.loc[tum_df['Treatment'] == "AL lifelong"].replace(100, 10)
DR_tum = tum_df.loc[tum_df["Treatment"] == "DR lifelong"].replace(100, 10)
sw52_tum = tum_df.loc[tum_df["Treatment"] == "AL_DR12M"].replace(100, 10)
sw69_tum = tum_df.loc[tum_df["Treatment"] == "AL_DR16M"].replace(100, 10)
sw87_tum = tum_df.loc[tum_df["Treatment"] == "AL_DR20M"].replace(100, 10)
sw104_tum = tum_df.loc[tum_df["Treatment"] == "AL_DR24M"].replace(100, 10)


# Get sum of columns per mouse, cohort and date of death in a new dataframe
tumsum_df = tum_df.loc[:, ["Treatment", "Cause_of_death"]]
tumsum_df["tumour_count"] = tum_df.drop(columns = "Cause_of_death").sum(axis = 1, numeric_only = True)
# }}}

def per_tumour_df(_tum, d):
    pert_ = pd.DataFrame()
    pert_["Age"] = _tum.loc[:, "Cause_of_death"]
    pert_["Eye"] = _tum.loc[:, ["Harderian_gland_(eye)_tumour_left", 
                                    "Harderian_gland_(eye)_tumour_right"]].sum(axis = 1)
    pert_["Connective_tissue"] = _tum.loc[:, ['Connective_tissue_tumour_in_left_WAT',
           'Connective_tissue_tumour_in_right_WAT',
           'Connective_tissue_tumour_in_Pancreas',
           'Connective_tissue_tumour_in_mesenteric_WAT_along_SI/Colon',
           'Connective_tissue_tumour_associated_with_liver/stomach',
           'Connective_tissue_tumour_associated_with_kidneys',
           'Connective_tissue_tumour_associated_with_left_uterine_horn/_ovary',
           'Connective_tissue_tumour_associated_with_right_uterine_horn_/ovary',
           'Connective_tissue_tumour_on_heart',
           'Connective_tissue_tumour_in_pleural_cavity/_on_sternum']].sum(axis = 1)
    pert_["Liver"] = _tum.loc[:, ['Liver_tumour_big_lobe', 'Liver_tumour_small_lobe']].sum(axis = 1)
    pert_["Kidney"] = _tum.loc[:, ['Left_kidney_tumour', 'Right_kidney_tumour']].sum(axis = 1)
    pert_["Ovary"] = _tum.loc[:, ['Right_ovary_tumour','Left_ovary_tumour']].sum(axis = 1)
    pert_["Uterine"] = _tum.loc[:, ['Left_uterine_tumour', 'Right_uterine_tumour']].sum(axis = 1)
    pert_["Lymph_node"] = _tum.loc[:, ['Ulcerated_lymph_node_tumour?', 'Lymph_node_tumour_on_back',
           'Lymph_node_tumour_left_hind_leg', 'Lymph_node_tumour_right_hind_leg',
           'Lymph_node_tumour_left_front_leg', 'Lymph_node_tumour_right_front_leg',
           'Lymph_node_tumour_left_side', 'Lymph_node_tumour_right_side',
           'Lymph_node_tumour_in_proximity_to_kidneys', 'Lymph_node_tumour_neck',
           'Lymph_node_tumour_groin_area', 'Lymph_node_tumour_in_abdominal_cavity',
           'Lymph_node_tumour_on_stomach/_sternum',
           'Lymph_node_tumour_in_pleural_cavity', 'Lymph_node_tumour_anus',
           'Lymph_node_tumour_on_face/_skull']].sum(axis = 1)
    pert_["Lung"] = _tum.loc[:, ['Left_lung_tumour','Right_lung_tumour']].sum(axis = 1)
    pert_["Colorectal"] = _tum.loc[:, ['Colorectal_tumour']].sum(axis = 1)
    pert_["Brain"] = _tum.loc[:, ['Brain_tumour']].sum(axis = 1)
    pert_["Heart"] = _tum.loc[:, ['Heart_tumour']].sum(axis = 1)
    pert_["diet"] = d
    pert_ = pert_.reset_index()
    return(pert_)


# {{{
pert_AL = per_tumour_df(AL_tum, "AL_lifelong")
pert_DR = per_tumour_df(DR_tum, "DR_lifelong")
pert_sw52 = per_tumour_df(sw52_tum, "AL_DR12M")
pert_sw69 = per_tumour_df(sw69_tum, "AL_DR16M")
pert_sw87 = per_tumour_df(sw87_tum, "AL_DR20M")
pert_sw104 = per_tumour_df(sw104_tum, "AL_DR24M")

all_pert = pd.concat([pert_AL, pert_DR, pert_sw52, pert_sw69, pert_sw87, pert_sw104])
# }}}

def make_counts(pert_, d):
    counts = pd.DataFrame(pert_["Age"])
    counts["counts"] = pert_.iloc[:, 2:].sum(axis = 1)
    counts["diet"] = d
    counts["Animal_ID"] = pert_["Animal_ID"]
    return(counts)


# {{{
countsAL = make_counts(pert_AL, "AL_lifelong")
countsDR = make_counts(pert_DR, "DR_lifelong")
counts52 = make_counts(pert_sw52, "AL_DR12M")
counts69 = make_counts(pert_sw69, "AL_DR16M")
counts87 = make_counts(pert_sw87, "AL_DR20M")
counts104 = make_counts(pert_sw104, "AL_DR24M")

all_counts = pd.concat([countsAL, countsDR, counts52, counts69, counts87, counts104])
# }}}








# ### Getting and formatting SPL alpha diversity

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

# Keep only q values 0, 1, 2
df = df[df["q"].isin([0, 1, 2])]

# Read metadata so we can get mouse name
nam = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8"])
nam.columns = ["Mouse", "Diet", "Age", "Illumina", "Illumina2", "Barcode", "sample_id"]

df2 = pd.merge(nam.loc[:, ["Mouse", "sample_id"]], df, on = "sample_id", how = "inner")
# Rename columns
df2.columns = ['Animal_ID', 'sample_id', 'biogroup', 'q', 'd', 'd_sd', 'd_lower',
       'd_upper', 'e', 'e_lower', 'e_upper', 'Age', 'Diet']
# }}}

all_counts

# ### Mix tumour counts and tumour per tissue with alpha diversity values

# {{{
tum_counts_alpha = pd.merge(all_counts.loc[:, ["Animal_ID", "counts"]], df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
shannon_df = tum_counts_alpha[tum_counts_alpha["q"] == 1]

pa_path_all = all_counts.loc[:, ["Animal_ID", "counts"]]
pa_path_all[pa_path_all["counts"] > 1] = 1
# Generate same dataframe for pert dataframes
pert_path_alpha = pd.merge(pa_path_all, df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pert_path_alpha.to_csv("{}/df_tumour_countsPA_alpha.csv".format(path), index = False)
# }}}

# {{{
sns.lmplot(data = shannon_df[shannon_df["Age"] == 24], x = "counts", y = "d", aspect = 1.2, )

plt.xlabel("Number of tumours", fontsize = 20)
plt.ylabel("Shannon Alpha diversity", fontsize = 20)
plt.tick_params(axis = "x", labelsize = 18)
plt.tick_params(axis = "y", labelsize = 18)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.show()
# }}}

# Calulating pearson correlation for this:
stats.pearsonr(shannon_df[shannon_df["Age"] == 24]["counts"], 
               shannon_df[shannon_df["Age"] == 24]["d"])

# Save tumour counts with alpha from spleen for analysis in R (tumour_alpha_stats.Rmd)
tum_counts_alpha.to_csv("{}/df_tumour_counts_alpha.csv".format(path), index = False)

pert_path_alpha.head()

# {{{
# Generate same dataframe for pert dataframes
pert_path_alpha = pd.merge(all_pert.drop(columns = ["Age", "diet"]), df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pert_path_alpha.to_csv("{}/df_tumtis_alpha.csv".format(path), index = False)

pa_path_all = all_pert.drop(columns = ["Age", "diet"])
pa_path_all["Eye"].values[pa_path_all["Eye"] > 1] = 1
pa_path_all["Connective_tissue"].values[pa_path_all["Connective_tissue"] > 1] = 1
pa_path_all["Liver"].values[pa_path_all["Liver"] > 1] = 1
pa_path_all["Kidney"].values[pa_path_all["Kidney"] > 1] = 1
pa_path_all["Ovary"].values[pa_path_all["Ovary"] > 1] = 1
pa_path_all["Uterine"].values[pa_path_all["Uterine"] > 1] = 1
pa_path_all["Lymph_node"].values[pa_path_all["Lymph_node"] > 1] = 1
pa_path_all["Lung"].values[pa_path_all["Lung"] > 1] = 1
pa_path_all["Colorectal"].values[pa_path_all["Colorectal"] > 1] = 1
pa_path_all["Brain"].values[pa_path_all["Brain"] > 1] = 1
pa_path_all["Heart"].values[pa_path_all["Heart"] > 1] = 1

# Generate same dataframe for pert dataframes
pert_path_alpha = pd.merge(pa_path_all, df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pert_path_alpha.to_csv("{}/df_tumtisPA_alpha.csv".format(path), index = False)
# }}}

# {{{
# Check linear regression within each diet, age associated with tumour
# }}}

def boxplots_cond(df_prov, order, str_vals, text):
    
    # Subset the info we need
    df = df_prov[df_prov["biogroup"].isin(order)]
    
    fig, ax = plt.subplots(figsize = (6,5))
    ax = sns.boxplot(data = df, x = "biogroup", y = "counts", palette = palette3, order = order, showfliers = False, ax = ax)
    ax = sns.swarmplot(data = df, x = "biogroup", y = "counts", order = order, dodge = True, color=".25", ax = ax)
    ax.tick_params(axis = "x", labelsize=18)
    ax.tick_params(axis = "y", labelsize=18)
    ax.set_xlabel("Age [Months]", fontsize = 20)
    ax.set_ylabel("Total tumour counts", fontsize = 20)

    ax.set_xticklabels(str_vals)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
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
        plt.savefig("{}tumour_{}Box_{}_{}.pdf".format(path, text, nam, new_day))
    else:
        plt.show()


order = ["5M DR_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
str_vals = [5, 20, 24]
boxplots_cond(shannon_df, order, str_vals, "Richness")

# {{{
order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_lifelong"]
str_vals = [5, 16, 20, 24]
boxplots_cond(shannon_df, order, str_vals, "Richness")

order = ["5M DR_lifelong", "16M AL_lifelong", "20M DR_lifelong", "24M DR_lifelong"]
boxplots_cond(shannon_df, order, str_vals, "Richness")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_DR16M", "24M AL_DR16M"]
boxplots_cond(shannon_df, order, str_vals, "Richness")

order = ["5M AL_lifelong", "16M AL_lifelong", "20M AL_lifelong", "24M AL_DR20M"]
boxplots_cond(shannon_df, order, str_vals, "Richness")
# }}}

# ## Getting and formatting dataframes for F2 pathology

# {{{
# Select columns that are not tumours
f1_col = list(metadata.columns)

tumour_col = []

for e in f1_col:
    if "tumour" in e:
        continue
    else:
        tumour_col.append(e)

tumour_col.append("Skull_bone_dissolved_due_to_eye_tumour")
tumour_col.append("Organs_dislocated/shifted_due_to_tumour")
tumour_col.remove("Date_of_pathology")
tumour_col.remove('Cohort')
tumour_col.remove('Date_of_Birth')
tumour_col.remove('Natural_death/killed')
tumour_col.remove('Samples_taken')
tumour_col.remove('Corpse_frozen')
tumour_col.remove('Comments')
tumour_col.remove("Post-mortem_cannibalism")
tumour_col.remove("Sarcoma")
# }}}

tum_df = metadata.loc[:, tumour_col]
tum_df = tum_df.replace(2, 1)


# {{{
def describe_tumours(pert_):
    _desc = pd.DataFrame()
    pert_ = pert_.drop(columns = ["Treatment"])
    _desc["Total"] = pert_.T.sum(axis = 1)
    _desc["Mean"] = pert_.T.mean(axis = 1)
    _desc["Max"] = pert_.T.max(axis = 1)


    m1 = pert_.values > 0.0
    df = pd.DataFrame(np.select([m1], [1], default=0), index=pert_.index, columns=pert_.columns)
    row_num, _ = df.shape
    _desc["Percentage_affected_mice"] = (df.T.sum(axis = 1)/row_num)*100
    _desc["Num_affected"] = df.T.sum(axis = 1)
    _desc["Num_not-affected"] = row_num - (df.T.sum(axis = 1))
    
    return(_desc)

def plot_tumours(pert_, diet, new_day):
    sns.boxplot(data=pert_.drop(columns = ["Cause_of_death", "Treatment"]))
    plt.title("{}\n".format(diet))
    plt.xticks(rotation = 90)
    plt.tight_layout()


    plt.show()
# }}}

d_al = describe_tumours(tum_df[tum_df["Treatment"] == "AL lifelong"])
d_dr = describe_tumours(tum_df[tum_df["Treatment"] == "DR lifelong"])
d_sw52 = describe_tumours(tum_df[tum_df["Treatment"] == "AL_DR12M"])
d_sw69 = describe_tumours(tum_df[tum_df["Treatment"] == "AL_DR16M"])
d_sw87 = describe_tumours(tum_df[tum_df["Treatment"] == "AL_DR20M"])
d_sw104 = describe_tumours(tum_df[tum_df["Treatment"] == "AL_DR24M"])

plot_tumours(tum_df[tum_df["Treatment"] == "AL lifelong"], "AL", new_day)
plot_tumours(tum_df[tum_df["Treatment"] == "DR lifelong"], "DR", new_day)
plot_tumours(tum_df[tum_df["Treatment"] == "AL_DR12M"], "AL_DR12M", new_day)
plot_tumours(tum_df[tum_df["Treatment"] == "AL_DR16M"], "AL_DR16M", new_day)
plot_tumours(tum_df[tum_df["Treatment"] == "AL_DR20M"], "AL_DR20M", new_day)
plot_tumours(tum_df[tum_df["Treatment"] == "AL_DR24M"], "AL_DR24M", new_day)

path_AL = tum_df[tum_df["Treatment"] == "AL lifelong"]
path_DR = tum_df[tum_df["Treatment"] == "DR lifelong"]
path_sw52 = tum_df[tum_df["Treatment"] == "AL_DR12M"]
path_sw69 = tum_df[tum_df["Treatment"] == "AL_DR16M"]
path_sw87 = tum_df[tum_df["Treatment"] == "AL_DR20M"]
path_sw104 = tum_df[tum_df["Treatment"] == "AL_DR24M"]

all_path = pd.concat([path_AL, path_DR, path_sw52, path_sw69, path_sw87, path_sw104])

all_path.columns

# {{{
# Generate same dataframe for pert dataframes
pert_path_alpha = pd.merge(all_path.drop(columns = ["Cause_of_death", "Treatment"]), df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pert_path_alpha.to_csv("{}/df_pathologyall_alpha.csv".format(path), sep = ";", index = False)

## Making also presence/absence dataframes
pa_path_all = all_path.drop(columns = ["Cause_of_death", "Treatment"])
pa_path_all[pa_path_all > 1] = 1
pa_all_path_alpha = pd.merge(pa_path_all, df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pa_all_path_alpha.to_csv("{}/df_pathologyallPA_alpha.csv".format(path), sep = ";", index = False)
# }}}

pa_all_path_alpha


def per_tumour_df(_tum):
    pert_ = pd.DataFrame()
    pert_["All_organs"] = _tum.loc[:, ['Autolysis_of_organs', 'All_organs_discoloured',
                                       'Gasping_for_air/_agonal_respiration','Cold_to_the_touch',
                                       'On_the_verge_of_dying', 'Other_cysts', 
                                       'Organs_dislocated/shifted_due_to_tumour']].sum(axis = 1)
    pert_["Spleen"] = _tum.loc[:, ['Enlarged_spleen', 'Small_spleen', 'Granulated_spleen', 
                                   'Spleen_discoloured_black']].sum(axis = 1)
    pert_["WAT"] = _tum.loc[:, ['Greyish/greenish/_blackish_WAT', 'WAT_brown', 'WAT_cyst', 'WAT_in_pleural_cavity', 
                               'Meager,_hardly_any_WAT']].sum(axis = 1)
    pert_["Pancreas"] = _tum.loc[:, ['Pancreas_knotty','Pancreas_glassy','Pancreas_pale/_discoloured',
                                     'Pancreas_black']].sum(axis = 1)
    pert_["Internal_bleeding"] = _tum.loc[:, ["Gastrointestinal_bleedings", "Fluid_filled_stomach", 
                                             "Bleedings_into_uterus", 'Internal_bleeding_into_abdominal_cavity', 
                                             'Internal_bleeding_into_pleural_cavity', 'Free_fluid_in_abdominal_cavity']].sum(axis = 1)
    pert_["Intestinal_tract"] = _tum.loc[:, ['Bloated_stomach', 'Bloated_SI_and_colon', 'Discoloured_SI', 
                                            'Bloody_poo', 'Discoloured_cecum', 'Indigestion']].sum(axis = 1)
    pert_["Liver"] = _tum.loc[:, ['Hepatomegalia','Liver_lobes_discoloured_black','Liver_pale',
                                  'Liver_appears_granulated','Liver_yellow_/_Jaundice']].sum(axis = 1)
    pert_["Kidneys"] = _tum.loc[:, ['Enlarged_kidneys', 'Kidney_spotty_or_granulated', 'Discoloured_kidneys', 
                                   'Kidneys_pale']].sum(axis = 1)
    pert_["Reproductive_tract"] = _tum.loc[:, ['Uterine_cyst_left_horn','Uterine_cyst_right_horn',
                                               'Uterus_discoloured_black', 'Cyst_left_ovary','Cyst_right_ovary', 
                                              'Agglutinated_vagina', 'Uterus_thickened_and_hard']].sum(axis = 1)
    pert_["gall_bladder"] = _tum.loc[:, ['Fully-filled_gall_bladder', 'Black_gall_bladder', 'Enlarged_bladder']].sum(axis = 1)
    pert_["Lungs"] = _tum.loc[:, ['Enlarged_lungs','Lungs_pale']].sum(axis = 1)
    pert_["Heart"] = _tum.loc[:, ["Enlarged_heart"]].sum(axis = 1)
    pert_["Physical_appearance"] = _tum.loc[:, ['Crooked_back/_bad_habitus', 'Hind_leg_paralysis', 'Critical_weight_loss', 
                                               'Ill-kept,_rumpled_furr', 'Open_wounds_on_tail']].sum(axis = 1)
    pert_["Lymph_vessels"] = _tum.loc[:, ['lymph_vessels_close_to_spine_thickened']].sum(axis = 1)
    pert_["Adrenal_gland"] = _tum.loc[:, ['enlarged_adrenal_gland']].sum(axis = 1)
    return(pert_)


# {{{
sum_al = per_tumour_df(tum_df[tum_df["Treatment"] == "AL lifelong"])
sum_dr = per_tumour_df(tum_df[tum_df["Treatment"] == "DR lifelong"])
sum_sw52 = per_tumour_df(tum_df[tum_df["Treatment"] == "AL_DR12M"])
sum_sw69 = per_tumour_df(tum_df[tum_df["Treatment"] == "AL_DR16M"])
sum_sw87 = per_tumour_df(tum_df[tum_df["Treatment"] == "AL_DR20M"])
sum_sw104 = per_tumour_df(tum_df[tum_df["Treatment"] == "AL_DR24M"])

sum_all = pd.concat([sum_al, sum_dr, sum_sw52, sum_sw69, sum_sw87, sum_sw104])
# }}}

# {{{
# Generate same dataframe for pert dataframes
sum_path_alpha = pd.merge(sum_all, df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
sum_path_alpha.to_csv("{}/df_pathologysum_alpha.csv".format(path), sep = ";", index = False)

## Making also presence/absence dataframes
pa_path_all = sum_all.copy()
pa_path_all[pa_path_all > 1] = 1
pa_path_alpha = pd.merge(pa_path_all, df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pa_path_alpha.to_csv("{}/df_pathologysumPA_alpha.csv".format(path), sep = ";", index = False)
# }}}




pa_path_alpha


# {{{
fig, ax = plt.subplots(figsize = (6,5))

sns.boxplot(data = pa_path_alpha[pa_path_alpha["q"] == 2], y = "d", x = "Spleen", showfliers = False)
sns.swarmplot(data = pa_path_alpha[pa_path_alpha["q"] == 2], y = "d", x = "Spleen", color = ".25")
plt.ylim(0, None)

ax.tick_params(axis = "x", labelsize=18)
ax.tick_params(axis = "y", labelsize=18)
ax.set_xlabel("Enlarged Spleen", fontsize = 20)
ax.set_ylabel("Simpson - Alpha diversity", fontsize = 20)


ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
    
if run_type != "dry":
    plt.savefig("{}Path_ESpleen_Simpson_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

pa_path_alpha[pa_path_alpha["q"] == 2]["Reproductive_tract"].min()


