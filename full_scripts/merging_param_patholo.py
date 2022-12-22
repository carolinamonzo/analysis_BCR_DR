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
import os

new_day = datetime.datetime.now().strftime("%Y%m%d")
path = "../../analysis/plots/alpha_phenotype/"

#run_type = "dry"
run_type = "wet"

#organ = "ILE"
organ = "SPL"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", "AL_DR16M":"teal", "AL_DR20M":"gold"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette1 = dict(zip(gr, co))

palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}
# }}}

# ## Getting and formatting F2 tumour load

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

# {{{
def make_counts(pert_, d):
    counts = pd.DataFrame(pert_["Age"])
    counts["counts"] = pert_.iloc[:, 2:].sum(axis = 1)
    counts["diet"] = d
    counts["Animal_ID"] = pert_["Animal_ID"]
    return(counts)

countsAL = make_counts(pert_AL, "AL_lifelong")
countsDR = make_counts(pert_DR, "DR_lifelong")
counts52 = make_counts(pert_sw52, "AL_DR12M")
counts69 = make_counts(pert_sw69, "AL_DR16M")
counts87 = make_counts(pert_sw87, "AL_DR20M")
counts104 = make_counts(pert_sw104, "AL_DR24M")

all_counts = pd.concat([countsAL, countsDR, counts52, counts69, counts87, counts104])
# }}}



# ## Getting and formatting SPL alpha diversity

# {{{
path_ori = "../../analysis/plots/"
files = []
for i in os.listdir(path_ori):
    if os.path.isfile(os.path.join(path_ori,i)) and i.startswith("alpha_values") and i.endswith("{}.tsv".format(organ)):
        files.append(i)

files = [x for x in files if "stats" not in x]
#files.remove('alpha_values_ILE_small.tsv')
files.remove('alpha_values_{}.tsv'.format(organ))
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
nam.columns = ["Animal_ID", "Diet", "Age", "Illumina", "Illumina2", "Barcode", "sample_id"]

df2 = pd.merge(nam.loc[:, ["Animal_ID", "sample_id", "Diet", "Age"]], df.loc[:, ["sample_id", "q", "d", "biogroup"]], on = "sample_id", how = "inner")
# }}}

# ## Mix tumour count and per tissue with alpha diversity

# {{{
tum_counts_alpha = pd.merge(all_counts.loc[:, ["Animal_ID", "counts"]], df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
shannon_df = tum_counts_alpha[tum_counts_alpha["q"] == 1]
richness_df = tum_counts_alpha[tum_counts_alpha["q"] == 0]
simpson_df = tum_counts_alpha[tum_counts_alpha["q"] == 2]

pa_path_all = all_counts.loc[:, ["Animal_ID", "counts"]]
pa_path_all[pa_path_all["counts"] > 1] = 1
# Generate same dataframe for pert dataframes
pert_path_alpha = pd.merge(pa_path_all, df2.loc[:, ["Animal_ID", "biogroup", "q", "d", "Age", "Diet"]], on = "Animal_ID", how = "inner")
# Save pathology with alpha from spleen for analysis in R (pathology_alpha_stats.Rmd)
pert_path_alpha.to_csv("{}/df_tumour_countsPA_alpha.csv".format(path), index = False)

PAshannon = pert_path_alpha[pert_path_alpha["q"] == 1]
PArichness = pert_path_alpha[pert_path_alpha["q"] == 0]
PAsimpson = pert_path_alpha[pert_path_alpha["q"] == 2]
# }}}

def correlation_numtum_alpha(dfa, alpha):

    sns.lmplot(data = dfa, x = "d", y = "counts", aspect = 1.2, )

    plt.xlabel("[{}] Alpha diversity".format(alpha), fontsize = 20)
    plt.ylabel("Number of tumours", fontsize = 20)
    plt.tick_params(axis = "x", labelsize = 18)
    plt.tick_params(axis = "y", labelsize = 18)
    plt.ylim(0, None)
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    if run_type != "dry":
        plt.savefig("{}corr_{}_numtum_{}.pdf".format(path, alpha, new_day))
    else:
        plt.show()
    # Calulating pearson correlation for this:
    print(stats.pearsonr(dfa["counts"], dfa["d"]))


correlation_numtum_alpha(shannon_df, "Shannon")
correlation_numtum_alpha(richness_df, "Richness")
correlation_numtum_alpha(simpson_df, "Simpson")

# Save tumour counts with alpha from spleen for analysis in R (tumour_alpha_stats.Rmd)
if run_type != "dry":
    tum_counts_alpha.to_csv("{}/df_tumour_counts_alpha.csv".format(path), index = False)

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
if run_type != "dry":
    pert_path_alpha.to_csv("{}/df_tumtisPA_alpha.csv".format(path), index = False)
# }}}

# ## Read clonal abundance

p20 = pd.read_csv("../../analysis/plots/clonal_abundance/{}_P20_values.csv".format(organ), sep = ",")
p20 = pd.merge(p20.loc[:, ["sample_id", "p"]], nam.loc[:, ["sample_id", "Animal_ID"]], on = "sample_id", how = "inner")
p20.columns = ["sample_id", "p20", "Animal_ID"]
p20.drop(columns = ["sample_id"], inplace = True)

# Merge this data into the tumour file
tum_alpha_p20 = pd.merge(pert_path_alpha, p20, on = "Animal_ID", how = "inner")

# ## Read mean SHM

# if its ileum, the date is the 6th
shm = pd.read_csv(f"../../analysis/results_tables/{organ}_frequency_values_20211107.csv", sep = ";")
shm.columns = ["Age", "mu_freq_seq_r", "mu_freq_seq_s", "sample_id", "Diet", "biogroup"]
shm = pd.merge(nam.loc[:, ["Animal_ID", "sample_id"]], shm, on = "sample_id", how = "inner")
shm["mu_freq"] = shm["mu_freq_seq_r"] + shm["mu_freq_seq_s"]
# Adding mutational status column
# 1 is mutated, 0 not mutated
shm["mut_status"] = np.where(shm["mu_freq"] >= 0.01, 1, 0)
shm.drop(columns = ["sample_id", "Diet", "Age", "biogroup", "mu_freq"], inplace = True)

tum_alpha_p20_shm = pd.merge(tum_alpha_p20, shm, on = "Animal_ID", how = "inner")

# ## Read CDR3

cdr3 = pd.read_csv("../../analysis/results_tables/CDR3_SPL_20220404.csv", sep = ",")
cdr3.columns = ["Animal_ID", "cdr3_mu", "cdr3_sigma"]
tum_alpha_p20_shm_cdr3 = pd.merge(tum_alpha_p20_shm, cdr3, on = "Animal_ID", how = "inner")

# ## Read ClassSwitch

cs = pd.read_csv("../../analysis/plots/naiveClassSwitch/SPL_NaiveClassw_20211109.csv", sep = ";")
cs = cs.drop(columns = "biogroup")
cs.columns = ['IgM-IgD-', 'IgM+IgD+SHM+', 'IgM+IgD+SHM-', 'Animal_ID']

tum_alpha_p20_shm_cdr3_cs = pd.merge(tum_alpha_p20_shm_cdr3, cs, on = "Animal_ID", how = "inner")

# ## Read RDI  
# RDI has one value per comparison, so we will use uniqueness as a measure of beta diversity

beta = pd.read_csv("../../analysis/results_tables/RDI_samp_SPL.tsv", sep = " ")
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
beta = beta[beta["value"] != 0.0]
# We remove the ones that the RDI is duplicated because is the same comparison
beta = beta.drop_duplicates()
beta = beta[~beta['col'].isin(['value'])]
beta = beta.drop(columns = "col")
beta = pd.DataFrame(beta.groupby(["row"])["value"].agg("min")).reset_index()
beta.columns = ["sample_id", "RDI_uniqueness"]

metadata = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8"])
metadata.columns = ["Mouse", "Diet", "Age", "Illumina", "Illumina2", "Barcode", "sample_id"]
beta = pd.merge(beta, metadata, on = ["sample_id"], how = "inner")
beta = beta.drop(columns = ["sample_id", "Diet", "Age", "Illumina", "Illumina2", "Barcode"])
beta.columns = ["RDI_uniqueness", "Animal_ID"]

tum_alpha_p20_shm_cdr3_cs_beta = pd.merge(tum_alpha_p20_shm_cdr3_cs, beta, on = "Animal_ID", how = "inner")

# ## Read V and J usage

# {{{
# Read file with frequencies
v = pd.read_csv(f"../../analysis/results_tables/Vgene_usage_SPL.tsv", sep = " ")
# Unifing subtypes of gene into original one

## DONT DO FOR ALLELES
v["gene"] = [x.split("S")[0] for x in v["gene"]]
v = v.groupby(["sample_id", "gene"]).aggregate({"seq_count":"sum", "seq_freq":"sum"}).reset_index()

df_matrix = v.loc[:, ["sample_id", "gene", "seq_freq"]].fillna(0).pivot(index = "sample_id", columns = "gene")
df_matrix.columns = df_matrix.columns.get_level_values(1)
df_matrix.fillna(0, inplace = True)
df_matrix = df_matrix.reset_index()
v = pd.merge(df_matrix, metadata, on = ["sample_id"], how = "inner")
v = v.drop(columns = ["sample_id", "Diet", "Age", "Illumina", "Illumina2", "Barcode"])
v.columns = ['IGHV1', 'IGHV10', 'IGHV11', 'IGHV12', 'IGHV13', 'IGHV14', 'IGHV15',
       'IGHV16', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7', 'IGHV8',
       'IGHV9', 'Animal_ID']

tum_alpha_p20_shm_cdr3_cs_beta_v = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta, v, on = "Animal_ID", how = "inner")
# }}}

# {{{
# Read file with frequencies
j = pd.read_csv(f"../../analysis/results_tables/Jgene_usage_SPL.tsv", sep = " ")
# Unifing subtypes of gene into original one

## DONT DO FOR ALLELES
j["gene"] = [x.split("S")[0] for x in j["gene"]]
j = j.groupby(["sample_id", "gene"]).aggregate({"seq_count":"sum", "seq_freq":"sum"}).reset_index()

df_matrix = j.loc[:, ["sample_id", "gene", "seq_freq"]].fillna(0).pivot(index = "sample_id", columns = "gene")
df_matrix.columns = df_matrix.columns.get_level_values(1)
df_matrix.fillna(0, inplace = True)
df_matrix = df_matrix.reset_index()
j = pd.merge(df_matrix, metadata, on = ["sample_id"], how = "inner")
j = j.drop(columns = ["sample_id", "Diet", "Age", "Illumina", "Illumina2", "Barcode"])
j.columns = ['IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'Animal_ID']

tum_alpha_p20_shm_cdr3_cs_beta_v_j = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v, j, on = "Animal_ID", how = "inner")
# }}}

# ## Isotypes relative frequency

organ = "SPL"
freqA = pd.read_csv(f"../../analysis/plots/isotypes_abundance/{organ}_IGA_freqValues.csv", sep = ";")
freqA = pd.merge(freqA.loc[:, ["sample_id", "rel"]], nam.loc[:, ["Animal_ID", "sample_id"]], on = "sample_id", how = "inner")
freqA.columns = ["sample_id", "relAb_IGA", "Animal_ID"]
freqA.drop(columns = ["sample_id"], inplace = True)
freqD = pd.read_csv(f"../../analysis/plots/isotypes_abundance/{organ}_IGD_freqValues.csv", sep = ";")
freqD = pd.merge(freqD.loc[:, ["sample_id", "rel"]], nam.loc[:, ["Animal_ID", "sample_id"]], on = "sample_id", how = "inner")
freqD.columns = ["sample_id", "relAb_IGD", "Animal_ID"]
freqD.drop(columns = ["sample_id"], inplace = True)
freqE = pd.read_csv(f"../../analysis/plots/isotypes_abundance/{organ}_IGE_freqValues.csv", sep = ";")
freqE = pd.merge(freqE.loc[:, ["sample_id", "rel"]], nam.loc[:, ["Animal_ID", "sample_id"]], on = "sample_id", how = "inner")
freqE.columns = ["sample_id", "relAb_IGE", "Animal_ID"]
freqE.drop(columns = ["sample_id"], inplace = True)
freqG = pd.read_csv(f"../../analysis/plots/isotypes_abundance/{organ}_IGG_freqValues.csv", sep = ";")
freqG = pd.merge(freqG.loc[:, ["sample_id", "rel"]], nam.loc[:, ["Animal_ID", "sample_id"]], on = "sample_id", how = "inner")
freqG.columns = ["sample_id", "relAb_IGG", "Animal_ID"]
freqG.drop(columns = ["sample_id"], inplace = True)
freqM = pd.read_csv(f"../../analysis/plots/isotypes_abundance/{organ}_IGM_freqValues.csv", sep = ";")
freqM = pd.merge(freqM.loc[:, ["sample_id", "rel"]], nam.loc[:, ["Animal_ID", "sample_id"]], on = "sample_id", how = "inner")
freqM.columns = ["sample_id", "relAb_IGM", "Animal_ID"]
freqM.drop(columns = ["sample_id"], inplace = True)

tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j, freqA, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso, freqD, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso, freqE, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso, freqG, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso, freqM, on = "Animal_ID", how = "inner")

# ## Alpha isotypes

# {{{
files = []
for i in os.listdir("../../analysis/plots/alpha_isotypes/"):
    if os.path.isfile(os.path.join("../../analysis/plots/alpha_isotypes/",i)) and i.startswith("alpha") and i.endswith("SPL_isosum.tsv") and "stats" not in i:
        files.append(i)
df = pd.DataFrame()
for e in files:
    s = e[9:]
    dt = pd.read_csv("{}{}".format("../../analysis/plots/alpha_isotypes/", e), sep = " ")
    dt["sample_id"] = s[:-15]
    df = pd.concat([df, dt])
    
df = pd.merge(df, metadata, on = ["sample_id"], how = "inner")
df = df.drop(columns = ["d_sd", "d_lower", "d_upper", "e", "e_lower", "e_upper", "sample_id", "Diet", "Age", "Illumina", "Illumina2", "Barcode"])
iga = df[df["sum_iso"] == "IGA"].drop(columns = ["sum_iso"])
igd = df[df["sum_iso"] == "IGD"].drop(columns = ["sum_iso"])
ige = df[df["sum_iso"] == "IGE"].drop(columns = ["sum_iso"])
igm = df[df["sum_iso"] == "IGM"].drop(columns = ["sum_iso"])
igg = df[df["sum_iso"] == "IGG"].drop(columns = ["sum_iso"])

iga.columns = ["q", "alpha_IGA", "Animal_ID"]
igd.columns = ["q", "alpha_IGD", "Animal_ID"]
ige.columns = ["q", "alpha_IGE", "Animal_ID"]
igm.columns = ["q", "alpha_IGM", "Animal_ID"]
igg.columns = ["q", "alpha_IGG", "Animal_ID"]
# }}}

tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso, iga, on = ["Animal_ID", "q"], how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso, igd, on = ["Animal_ID", "q"], how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso, ige, on = ["Animal_ID", "q"], how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso, igm, on = ["Animal_ID", "q"], how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso, igg, on = ["Animal_ID", "q"], how = "inner")

# ## Beta isotypes

iga = pd.read_csv("../../analysis/results_tables/RDI_IGA_samp_SPL_isosum.tsv", sep = " ")
igd = pd.read_csv("../../analysis/results_tables/RDI_IGD_samp_SPL_isosum.tsv", sep = " ")
ige = pd.read_csv("../../analysis/results_tables/RDI_IGE_samp_SPL_isosum.tsv", sep = " ")
igm = pd.read_csv("../../analysis/results_tables/RDI_IGM_samp_SPL_isosum.tsv", sep = " ")
igg = pd.read_csv("../../analysis/results_tables/RDI_IGG_samp_SPL_isosum.tsv", sep = " ")
# We remove the ones that the RDI is == 0 because that is comparison of a sample against itself
iga = iga[iga["value"] != 0.0].drop_duplicates()
igd = igd[igd["value"] != 0.0].drop_duplicates()
ige = ige[ige["value"] != 0.0].drop_duplicates()
igm = igm[igm["value"] != 0.0].drop_duplicates()
igg = igg[igg["value"] != 0.0].drop_duplicates()
# We remove the ones that the RDI is duplicated because is the same comparison
iga = iga[~iga['col'].isin(['value'])]
igd = igd[~igd['col'].isin(['value'])]
ige = ige[~ige['col'].isin(['value'])]
igm = igm[~igm['col'].isin(['value'])]
igg = igg[~igg['col'].isin(['value'])]
iga = iga.drop(columns = "col")
igd = igd.drop(columns = "col")
ige = ige.drop(columns = "col")
igm = igm.drop(columns = "col")
igg = igg.drop(columns = "col")
iga = pd.DataFrame(iga.groupby(["row"])["value"].agg("min")).reset_index()
igd = pd.DataFrame(igd.groupby(["row"])["value"].agg("min")).reset_index()
ige = pd.DataFrame(ige.groupby(["row"])["value"].agg("min")).reset_index()
igm = pd.DataFrame(igm.groupby(["row"])["value"].agg("min")).reset_index()
igg = pd.DataFrame(igg.groupby(["row"])["value"].agg("min")).reset_index()
iga.columns = ["sample_id", "RDI_uniqueness_IGA"]
igd.columns = ["sample_id", "RDI_uniqueness_IGD"]
ige.columns = ["sample_id", "RDI_uniqueness_IGE"]
igm.columns = ["sample_id", "RDI_uniqueness_IGM"]
igg.columns = ["sample_id", "RDI_uniqueness_IGG"]

# {{{
metadata = metadata.drop(columns = ["Diet", "Age", "Illumina", "Illumina2", "Barcode"])
metadata.columns = ["Animal_ID", "sample_id"]
iga = pd.merge(metadata, iga, on = "sample_id", how = "inner").drop(columns = "sample_id")
igd = pd.merge(metadata, igd, on = "sample_id", how = "inner").drop(columns = "sample_id")
ige = pd.merge(metadata, ige, on = "sample_id", how = "inner").drop(columns = "sample_id")
igm = pd.merge(metadata, igm, on = "sample_id", how = "inner").drop(columns = "sample_id")
igg = pd.merge(metadata, igg, on = "sample_id", how = "inner").drop(columns = "sample_id")

tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso, iga, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso, igd, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso, ige, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso, igm, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso, igg, on = "Animal_ID", how = "inner")
# }}}

# ## Clonal abundance isotype

# {{{
igA = pd.read_csv("../../analysis/plots/clonal_abundance/" + 'clab_IGA_{}.tsv'.format("SPL"), sep = " ")
igM = pd.read_csv("../../analysis/plots/clonal_abundance/" + 'clab_IGM_{}.tsv'.format("SPL"), sep = " ")
igD = pd.read_csv("../../analysis/plots/clonal_abundance/" + 'clab_IGD_{}.tsv'.format("SPL"), sep = " ")
igE = pd.read_csv("../../analysis/plots/clonal_abundance/" + 'clab_IGE_{}.tsv'.format("SPL"), sep = " ")
igG12 = pd.read_csv("../../analysis/plots/clonal_abundance/" + 'clab_IGG12_{}.tsv'.format("SPL"), sep = " ")
igG3 = pd.read_csv("../../analysis/plots/clonal_abundance/" + 'clab_IGG3_{}.tsv'.format("SPL"), sep = " ")

igG = pd.concat([igG12, igG3])

igA["rank"] = igA["rank"].astype(float)
igA["p"] = igA["p"].astype(float)
p20_igA = igA[igA["rank"] <= 20]
p20_igA = pd.DataFrame(p20_igA.loc[:, ["sample_id", "p"]].groupby(['sample_id'])["p"].sum())
p20_igA.reset_index(inplace = True)
p20_igA = pd.merge(p20_igA, metadata, on = "sample_id", how = "inner").drop(columns = "sample_id")
p20_igA.columns = ["p20_IGA", "Animal_ID"]

igM["rank"] = igM["rank"].astype(float)
igM["p"] = igM["p"].astype(float)
p20_igM = igM[igM["rank"] <= 20]
p20_igM = pd.DataFrame(p20_igM.loc[:, ["sample_id", "p"]].groupby(['sample_id'])["p"].sum())
p20_igM.reset_index(inplace = True)
p20_igM = pd.merge(p20_igM, metadata, on = "sample_id", how = "inner").drop(columns = "sample_id")
p20_igM.columns = ["p20_IGM", "Animal_ID"]

igD["rank"] = igD["rank"].astype(float)
igD["p"] = igD["p"].astype(float)
p20_igD = igD[igD["rank"] <= 20]
p20_igD = pd.DataFrame(p20_igD.loc[:, ["sample_id", "p"]].groupby(['sample_id'])["p"].sum())
p20_igD.reset_index(inplace = True)
p20_igD = pd.merge(p20_igD, metadata, on = "sample_id", how = "inner").drop(columns = "sample_id")
p20_igD.columns = ["p20_IGD", "Animal_ID"]

igE["rank"] = igE["rank"].astype(float)
igE["p"] = igE["p"].astype(float)
p20_igE = igE[igE["rank"] <= 20]
p20_igE = pd.DataFrame(p20_igE.loc[:, ["sample_id", "p"]].groupby(['sample_id'])["p"].sum())
p20_igE.reset_index(inplace = True)
p20_igE = pd.merge(p20_igE, metadata, on = "sample_id", how = "inner").drop(columns = "sample_id")
p20_igE.columns = ["p20_IGE", "Animal_ID"]

igG["rank"] = igG["rank"].astype(float)
igG["p"] = igG["p"].astype(float)
p20_igG = igG[igG["rank"] <= 20]
p20_igG = pd.DataFrame(p20_igG.loc[:, ["sample_id", "p"]].groupby(['sample_id'])["p"].sum())
p20_igG.reset_index(inplace = True)
p20_igG = pd.merge(p20_igG, metadata, on = "sample_id", how = "inner").drop(columns = "sample_id")
p20_igG.columns = ["p20_IGG", "Animal_ID"]
# }}}

tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso, p20_igA, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso, p20_igM, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso, p20_igD, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso, p20_igE, on = "Animal_ID", how = "inner")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso, p20_igG, on = "Animal_ID", how = "inner")

# ## SHM isotypes

shm = pd.read_csv("../../analysis/results_tables/SHM_SPL_isotypes_sample_20220404.csv", sep = ",")
shm = pd.merge(shm, metadata, on = "sample_id", how = "inner").drop(columns = "sample_id")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso, shm, on = "Animal_ID", how = "inner")

# ## CDR3 isotypes

cdr3 = pd.read_csv("../../analysis/results_tables/cdr3_isotypes_20220404.csv", sep = ",").drop(columns = "Unnamed: 0")
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso = pd.merge(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso, cdr3, on = "Animal_ID", how = "inner")

tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso.columns


# # Save the dataframe with everything and tumors

# {{{
def replace01(df, pheno, s):
# Replace to 0 and 1
    df = df.set_index("Animal_ID")
    for e in pheno:
        l = []
        for t in df[e]:
            if t < 1:
                l.append(0)
            else:
                l.append(1)
        df[e] = l

    # 0 absence, 1 one organ, 2 two or more organs
    classif = []
    for e in df.iloc[:, 3:].sum(axis = 1):
        if e <= 1:
            classif.append(e)
        else:
            classif.append(2)
    df[s] = classif
    return(df)

smpheno = ['Eye', 'Connective_tissue', 'Liver', 'Kidney', 'Ovary',
       'Uterine', 'Lymph_node', 'Lung', 'Colorectal', 'Brain', 'Heart']
tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso = replace01(tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso, smpheno, "presAbsTumor")
# }}}

tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso.to_csv("../../analysis/results_tables/tumor_allBCRmetrics_{}.csv".format(new_day), sep = ";")

alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso = tum_alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso.drop(columns = ['Eye', 'Connective_tissue', 'Liver', 'Kidney', 'Ovary',
       'Uterine', 'Lymph_node', 'Lung', 'Colorectal', 'Brain', 'Heart', "presAbsTumor"])

# # Getting the pathology to do the same

# {{{
## Since this should be systemic comparisons, we are doing spleen vs el_diablo_F2
metadata = pd.read_csv("../../metadata/F2-cross-sectional-elDiablo.csv", sep = ";", index_col = "Animal_ID", nrows=288).fillna(0)
metadata = metadata.replace(1000, 10)

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

dic_diet = {"chronic_AL":"AL lifelong", "chronic_DR":"DR lifelong", "AL_DR52W":"AL_DR12M", 
            "AL_DR69W":"AL_DR16M", "AL_DR87W":"AL_DR20M", "AL_DR104W":"AL_DR24M"}
metadata["Treatment"] = metadata["Cohort"].map(dic_diet)

# Total tumour counts
df_col = list(metadata.columns)
# }}}

# {{{
tumour_col = []

for e in df_col:
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
tum_df = metadata.loc[:, tumour_col]
tum_df = tum_df.replace(2, 1)
# }}}

tum_df.columns = ["Cause_of_death", 'Autolysis_of_organs', 'All_organs_discoloured',
       'Meager_hardly_any_WAT', 'Crooked_back_bad_habitus',
       'Enlarged_spleen', 'Small_spleen', 'Granulated_spleen',
       'Spleen_discoloured_black', 'Greyish_greenish_blackish_WAT',
       'WAT_brown', 'Pancreas_knotty', 'Pancreas_glassy',
       'Pancreas_pale_discoloured', 'Pancreas_black',
       'Gastrointestinal_bleedings', 'Bloated_stomach', 'Fluid_filled_stomach',
       'Bloated_SI_and_colon', 'Discoloured_SI', 'Hepatomegalia',
       'Liver_lobes_discoloured_black', 'Liver_pale',
       'Liver_appears_granulated', 'Liver_yellow_Jaundice',
       'Uterine_cyst_left_horn', 'Uterine_cyst_right_horn',
       'Uterus_discoloured_black', 'Bleedings_into_uterus', 'Cyst_left_ovary',
       'Cyst_right_ovary', 'Enlarged_heart', 'Enlarged_kidneys',
       'Fully_filled_gall_bladder', 'Internal_bleeding_into_abdominal_cavity',
       'Internal_bleeding_into_pleural_cavity', 'WAT_cyst',
       'Kidney_spotty_or_granulated', 'Hind_leg_paralysis',
       'Agglutinated_vagina', 'Gasping_for_air_agonal_respiration',
       'Cold_to_the_touch', 'On_the_verge_of_dying', 'Critical_weight_loss',
       'enlarged_adrenal_gland', 'Bloody_poo', 'Black_gall_bladder',
       'Discoloured_cecum', 'Indigestion', 'Other_cysts',
       'WAT_in_pleural_cavity', 'lymph_vessels_close_to_spine_thickened',
       'Ill_kept_rumpled_furr', 'Uterus_thickened_and_hard',
       'Discoloured_kidneys', 'Kidneys_pale', 'Free_fluid_in_abdominal_cavity',
       'Enlarged_bladder', 'Open_wounds_on_tail', 'Enlarged_lungs',
       'Lungs_pale', "Treatment", 'Skull_bone_dissolved_due_to_eye_tumour',
       'Organs_dislocated_shifted_due_to_tumour']


def per_tumour_df(_tum):
    pert_ = pd.DataFrame()
    pert_["All_organs"] = _tum.loc[:, ['Gasping_for_air_agonal_respiration','Cold_to_the_touch',
                                       'On_the_verge_of_dying', 'Other_cysts', 
                                       'Organs_dislocated_shifted_due_to_tumour']].sum(axis = 1)
    pert_["Spleen"] = _tum.loc[:, ['Enlarged_spleen', 'Small_spleen', 'Granulated_spleen']].sum(axis = 1)
    pert_["WAT"] = _tum.loc[:, ['WAT_brown', 'WAT_cyst', 'WAT_in_pleural_cavity', 
                               'Meager_hardly_any_WAT']].sum(axis = 1)
    pert_["Pancreas"] = _tum.loc[:, ['Pancreas_knotty','Pancreas_pale_discoloured']].sum(axis = 1)
    pert_["Internal_bleeding"] = _tum.loc[:, ["Gastrointestinal_bleedings", "Fluid_filled_stomach", 
                                             "Bleedings_into_uterus", 'Internal_bleeding_into_abdominal_cavity', 
                                             'Internal_bleeding_into_pleural_cavity', 'Free_fluid_in_abdominal_cavity']].sum(axis = 1)
    pert_["Intestinal_tract"] = _tum.loc[:, ['Bloated_stomach', 'Bloated_SI_and_colon', 
                                            'Bloody_poo', 'Indigestion']].sum(axis = 1)
    pert_["Liver"] = _tum.loc[:, ['Hepatomegalia','Liver_pale',
                                  'Liver_appears_granulated','Liver_yellow_Jaundice']].sum(axis = 1)
    pert_["Kidneys"] = _tum.loc[:, ['Enlarged_kidneys', 'Kidney_spotty_or_granulated', 'Discoloured_kidneys', 
                                   'Kidneys_pale']].sum(axis = 1)
    pert_["Reproductive_tract"] = _tum.loc[:, ['Uterine_cyst_left_horn','Uterine_cyst_right_horn',
                                               'Cyst_left_ovary','Cyst_right_ovary', 
                                              'Agglutinated_vagina', 'Uterus_thickened_and_hard']].sum(axis = 1)
    pert_["gall_bladder"] = _tum.loc[:, ['Fully_filled_gall_bladder', 'Black_gall_bladder', 'Enlarged_bladder']].sum(axis = 1)
    pert_["Lungs"] = _tum.loc[:, ['Enlarged_lungs','Lungs_pale']].sum(axis = 1)
    pert_["Heart"] = _tum.loc[:, ["Enlarged_heart"]].sum(axis = 1)
    pert_["Physical_appearance"] = _tum.loc[:, ['Crooked_back_bad_habitus', 'Hind_leg_paralysis', 'Critical_weight_loss', 
                                               'Ill_kept_rumpled_furr', 'Open_wounds_on_tail']].sum(axis = 1)
    pert_["Lymph_vessels"] = _tum.loc[:, ['lymph_vessels_close_to_spine_thickened']].sum(axis = 1)
    pert_["Adrenal_gland"] = _tum.loc[:, ['enlarged_adrenal_gland']].sum(axis = 1)
    return(pert_)


path_df = per_tumour_df(tum_df)

# {{{
path_mer = pd.merge(path_df.reset_index(), alpha_p20_shm_cdr3_cs_beta_v_j_reliso_alphiso_betiso_p20iso_shmiso_cdr3iso, on = "Animal_ID", how = "inner")

smpheno = ['All_organs', 'Spleen', 'WAT', 'Pancreas',
       'Internal_bleeding', 'Intestinal_tract', 'Liver', 'Kidneys',
       'Reproductive_tract', 'gall_bladder', 'Lungs', 'Heart',
       'Physical_appearance', 'Lymph_vessels', 'Adrenal_gland']
path_mer = replace01(path_mer, smpheno, "presAbsPatho")
# }}}

path_mer.to_csv("../../analysis/results_tables/pathology_allBCRmetrics_{}.csv".format(new_day), sep = ";")




