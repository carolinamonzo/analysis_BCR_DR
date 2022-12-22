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
import scipy.stats as stats
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm
from IPython.display import display
import scikit_posthocs as sp

new_day = datetime.datetime.now().strftime("%Y%m%d")
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
path = "../../analysis/plots/QC/"

#organ = "ILE"
organ = "SPL"

run_type = "dry"
#run_type = "wet"

palette2 = {"DR_lifelong":"red", "AL_lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"tab:green", "AL_DR20M":"yellow", "AL_DR24M":"grey"}
palette3 = {"5M AL_lifelong":"#D3E0F1", "16M AL_lifelong":"#A2C2DD", "20M AL_lifelong":"#6897C6", "24M AL_lifelong":"#3D64A8", "5M DR_lifelong":"#ECBAA1", 
            "20M DR_lifelong":"#DD694D", "24M DR_lifelong":"#AD1E22","20M AL_DR16M":"#779D45", "24M AL_DR16M":"#416F6F", "24M AL_DR20M":"#EDD853"}

iso_palette = {"IGA":"tab:blue", "IGD":"tab:orange", "IGE":"tab:green", "IGG":"tab:red", "IGM":"tab:brown"}
# }}}

# {{{
# Read metadata so we can get age and diet from there
metadata = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";").drop(columns = ["Unnamed: 7", "Unnamed: 8", "Illumina", "Illumina2", "Barcode"])
metadata.columns = ["Mouse", "Diet", "Age", "sample_id"]

metadata['Diet'] = metadata['Diet'].str.replace(' ','_')

metadata["biogroup"] = metadata["Age"].astype(str) + metadata["Diet"].astype(str)
# }}}

# Read full file to get counts
df = pd.read_csv("../../analysis/created_germlines/reannotated_db/merged_changeo_{}_isosum.tsv".format(organ), sep = "\t")

# {{{
numseq_samp = []
numseq_seq = []

for samp in metadata["sample_id"].to_list():
    numseq_samp.append(samp)
    numseq_seq.append(len(df[df["sample_id"] == samp]["sequence_id"].to_list()))
# }}}

dfseq = pd.DataFrame({"sample_id":numseq_samp, "numseq":numseq_seq})
dfseq= dfseq[dfseq['numseq'] != 0]

dfseq = metadata.merge(dfseq, on = "sample_id", how = "inner")

sns.violinplot(data = dfseq, x = "Age", y = "numseq", hue = "Diet", palette = palette2)

print(stats.kruskal(dfseq["numseq"].to_list(), dfseq["biogroup"].to_list()))
sp.posthoc_mannwhitney(dfseq, val_col = "numseq", group_col = "biogroup")

dfseq.describe()


