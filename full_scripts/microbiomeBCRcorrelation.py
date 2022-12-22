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
path = "../../analysis/plots/path_BCR/"

run_type = "dry"
#run_type = "wet"

#organ = "ILE"
organ = "SPL"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", "AL_DR16M":"teal", "AL_DR20M":"gold"}

gr = ['5AL_lifelong', '16AL_lifelong', '20AL_lifelong', '24AL_lifelong', '5DR_lifelong', '20DR_lifelong', '24DR_lifelong','20AL_DR16M','24AL_DR16M', '24AL_DR20M']
co = ["#D3E0F1", "#A2C2DD", "#6897C6", "#3D64A8", "#ECBAA1", "#DD694D", "#AD1E22", "#779D45", "#416F6F", "#EDD853"]
palette3 = dict(zip(gr, co))
# }}}

# {{{
df = pd.read_csv("../../analysis/results_tables/allBCRmetrics_tumpathoscore_{}.csv".format(new_day), sep = ";")
df = df.drop(columns = "Unnamed: 0")
df = df.sort_values(by = "Animal_ID", ascending = False)
df = df.set_index("Animal_ID")
microbiome = pd.read_csv("../../metadata/norm-CLEAN_ASVs_counts_merged_20210615.tsv", sep = "\t")
otus = microbiome["OTU"].to_list()
microbiome = microbiome.set_index("OTU").T
microbiome.reset_index(inplace = True)
microbiome["Animal_ID"] = microbiome["index"]
microbiome = microbiome.drop(columns = "index")
midrobiome = microbiome[microbiome["Animal_ID"].isin(df.index.to_list())]
microbiome = microbiome.sort_values(by = "Animal_ID", ascending = False)
microbiome = microbiome.set_index("Animal_ID")


# }}}

df1 = df.drop(columns = ["biogroup", "Age", "Diet"])
df2 = microbiome
cordf = pd.concat([df1, df2], axis=1, keys=['df1', 'df2']).corr().loc['df2', 'df1'].dropna(axis = 0, how = "all")
cordf = cordf.sort_values(by = "p20", ascending = False)

sns.clustermap(data = cordf.fillna(0))


