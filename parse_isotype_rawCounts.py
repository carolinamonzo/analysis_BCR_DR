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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# {{{
import os
import datetime
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

new_day = datetime.datetime.now().strftime("%Y%m%d")
path = "../analysis/plots/isotype_raw/"

run_type = "dry"
#run_type = "wet"

palette2 = {"DR lifelong":"red", "AL lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"teal", "AL_DR20M":"gold", "AL_DR24M":"grey"}
# }}}

# {{{
## Reading metadata files
eD = pd.read_csv("../metadata/F2_tissue_collection-cross-sectional-pathology.csv", sep = ";", nrows = 288)

nam = pd.read_csv("../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 5", "Unnamed: 6", "Unnamed: 7"])
nam.columns = ["Animal_ID", "Illumina", "Illumina2", "Barcode", "index"]

# Format metadata cross sectional
tomod = eD["Cause_of_death"].tolist()
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

eD["month"] = mod
eD["month"] = eD["month"].astype(float)

# And treatment from metadata cross sectional
dic_diet = {"chronic_AL":"AL lifelong", "chronic_DR":"DR lifelong", "AL_DR52W":"AL_DR12M", 
            "AL_DR69W":"AL_DR16M", "AL_DR87W":"AL_DR20M", "AL_DR104W":"AL_DR24M"}
eD["Treatment"] = eD["Cohort"].map(dic_diet)

## Select only Ileum samples

nam = nam[nam['Illumina'].str.contains("I")]
nam["index"] = [sub.replace("_isotype.txt", "") for sub in list(nam["index"])]

## Merge into one metadata dataframe

metadata = pd.DataFrame()
metadata = pd.merge(eD.loc[:, ["Animal_ID", "Treatment", "month"]], nam.loc[:, ["Animal_ID", "index"]], on = "Animal_ID", how = "inner")
metadata.index = metadata["index"]
# }}}

# {{{
file_path = "/beegfs/group_lp/home/CMonzo/CM_IGseq/analysis/isotype_raw/"
fof = "isotype_rawCounts.fof"

## Read files with stats
stat_files = []

# get a list of files to parse
with open(file_path + fof, "r") as fof_list:
    for f in fof_list:
        stat_files.append(file_path + f.strip("\n"))

# Merge into one dataframe
df = pd.DataFrame()

df = pd.read_csv(stat_files[0], delim_whitespace = True)
df = df.reset_index()
co = []
co.append(df.columns[1])
co.append(df.columns[0])
df.columns = co

for fi in stat_files[1:]:
    t = pd.read_csv(fi, delim_whitespace = True)
    t = t.reset_index()
    co = []
    co.append(t.columns[1])
    co.append(t.columns[0])
    t.columns = co
    
    # Merge
    df = pd.merge(df, t, on = "index", how = "outer")

df = df.set_index("index").T

df.index = df.index.str[17:]
df.index = [sub.replace(".fasta", "") for sub in list(df.index)]
# }}}

# {{{
# Merge with metadata

df = pd.merge(df, metadata.drop(columns = ["index"]), how = "inner", left_index = True, right_index = True)

## Normalize df per total reads

df["Total"] = df.sum(axis = 1, numeric_only = True)

df["IGG"] = df["IGG12"] + df["IGG3"]

norm = pd.DataFrame()

for e in ["IGA", "IGD", "IGE", "IGG", "IGM"]:
    norm[e] = df[e]/df["Total"]
    
# Add missed columns
norm["Animal_ID"] = df["Animal_ID"]
norm["Treatment"] = df["Treatment"]
norm["month"] = df["month"]

norm.index = norm["Animal_ID"]
norm = norm.drop(columns = ["Animal_ID"])
# }}}



def violin_IG(data, text):
    fig, ax = plt.subplots(figsize = (10,6))
    ax = sns.violinplot(data = data, x = "month", y = text, hue = "Treatment", palette = palette2, hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"])
    ax = sns.swarmplot(data = data, x = "month", y = text, hue = "Treatment", hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"], color = ".25", dodge = True)
    
    ax.tick_params(axis = "x", labelsize = 18)
    ax.tick_params(axis = "y", labelsize = 18)
    ax.set_ylabel("Normalized {} amount".format(text), fontsize = 20)
    ax.set_xlabel("Month", fontsize = 20)
    plt.title(text, fontsize = 22)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4], bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0., fontsize = 18)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    plt.show()


violin_IG(norm, "IGA")
violin_IG(norm, "IGM")
violin_IG(norm, "IGG")
violin_IG(norm, "IGD")
violin_IG(norm, "IGE")


def violin_IG_treat(data, text):
    fig, ax = plt.subplots(figsize = (10,6))
    ax = sns.lineplot(data = data, x = "month", y = text, hue = "Treatment", palette = palette2, hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"], linewidth = 3.5)
    
    ax.tick_params(axis = "x", labelsize = 18)
    ax.tick_params(axis = "y", labelsize = 18)
    ax.set_ylabel("Normalized {} amount".format(text), fontsize = 20)
    ax.set_xlabel("Month", fontsize = 20)
    plt.title(text, fontsize = 22)
    
    ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0., fontsize = 18)
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    plt.show()


violin_IG_treat(norm, "IGA")
violin_IG_treat(norm, "IGM")
violin_IG_treat(norm, "IGG")
violin_IG_treat(norm, "IGD")
violin_IG_treat(norm, "IGE")


