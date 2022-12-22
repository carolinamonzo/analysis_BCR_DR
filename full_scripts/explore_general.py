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

new_day = datetime.datetime.now().strftime("%Y%m%d")
path = "../plots/exploratory_data/"

run_type = "dry"
#run_type = "wet"

organ = "ILE"
#organ = "SPL"

palette2 = {"DR lifelong":"red", "AL lifelong":"dodgerblue", 
            "AL_DR12M":"magenta", "AL_DR16M":"teal", "AL_DR20M":"gold", "AL_DR24M":"grey"}
# }}}

# Reading metadata file
eD = pd.read_csv("../../metadata/SampleSheet_IGSeq.csv", sep = ";", usecols = ["Mouse", "Diet", 
                "Age", "Illumina", "Illumina2", "Barcode", "File"])
eD.head()

# Read Spleen hierarchical clustering
spl = pd.read_csv("../productive_vdj/SPL_clustering_table.tsv", sep = "\t")
# Read Ileum hierarchical clustering
ile = pd.read_csv("../productive_vdj/ILE_clustering_table.tsv", sep = "\t")

ilem = pd.merge(eD, ile, on = ["File"], how = "inner")
splm = pd.merge(eD, spl, on = ["File"], how = "inner")


def plot_threshold_method(df, text):
    fig, ax = plt.subplots(figsize = (10,6))
    ax = sns.swarmplot(data = df, x = "Age", y = "Method", hue = "Diet", palette = palette2, 
                       hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"], size = 10, dodge = True, alpha = 0.7, linewidth = 1.0)
    ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0., fontsize = 18)
    ax.tick_params(axis = "x", labelsize = 18)
    ax.tick_params(axis = "y", labelsize = 18)
    ax.set_xlabel("Age (M)", fontsize = 20)
    ax.set_ylabel("")
    plt.title("{} - Clonal threshold method".format(text), fontsize = 22)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    plt.show()


plot_threshold_method(ilem, "Ileum")
plot_threshold_method(splm, "Spleen")


# {{{
def countplot(x, hue, **kwargs):
    sns.countplot(x=x, hue=hue, **kwargs)
    
def plot_bars_method(df, text):    
    
    grid = sns.FacetGrid(data=df,row='Method',size=4,aspect=1)
    fig = grid.map(countplot,'Age','Diet',palette=palette2,order=[5, 16, 20, 24], 
                   hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"])
    fig.add_legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0., fontsize = 18)

    fig.axes[0,0].set_ylabel("Spectral clustering")
    fig.axes[1,0].set_ylabel("Hierarchical clustering")

    for ax in fig.axes.flat:
        if ax.get_title():
            ax.set_title("")
        ax.tick_params(axis = "x", labelsize = 18)
        ax.tick_params(axis = "y", labelsize = 18)
        ax.set_xlabel(ax.get_xlabel(), fontsize=16)
        ax.set_ylabel(ax.get_ylabel(), fontsize=16)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    
    if run_type != "dry":
        plt.savefig("{}{}_thresholdMethod_{}.pdf".format(path, text, new_day))
    else:
        plt.show()
# }}}


plot_bars_method(ilem, "Ileum")
plot_bars_method(splm, "Spleen")

prop = pd.DataFrame({">0":[2, 43], "=0":[48, 6], "Organ":["Ileum", "Spleen"]})


# {{{
def stacked_percent_df(df, diets):
    r = [0,1]

    # From raw value to percentage
    totals = [i+j for i,j in zip(df['=0'], df['>0'])]
    Bars0 = [i / j * 100 for i,j in zip(df['=0'], totals)]
    Bars1 = [i / j * 100 for i,j in zip(df['>0'], totals)]


    # Make dataframe with data, percentage of mice with x number of tumours per diet
    st = pd.DataFrame([Bars0, Bars1], columns = diets)
    st.index = ["=0", ">0",]
    return(st)

def make_staked_plot(df):
    fig, ax = plt.subplots(figsize=(6,4))
    df.T.plot(kind = "bar", stacked = True, ax = ax, color = ["blue", "orange"], width = 0.85, edgecolor = "black", linewidth = 1.5)
    plt.xticks(rotation = 0)
    ax.tick_params(axis = "x", labelsize=20)
    ax.tick_params(axis = "y", labelsize=20)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
    plt.title("Proportion of samples with novel alleles \n", fontsize = 18)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.tight_layout()
    
    if run_type != "dry":
        plt.savefig("{}proportionMice_novelAlleles_{}.pdf".format(path, new_day))
    else:
        plt.show()
# }}}


st = stacked_percent_df(prop.drop(columns = ["Organ"]), ["Ileum", "Spleen"])
make_staked_plot(st)

# {{{
# Read novel alleles file
nal = pd.read_csv("../productive_vdj/novel_alleles/novel_alleles_all.tsv", sep = " ")
nal.reset_index(inplace = True)
nal = nal.drop(columns = ["germline_call", "note", "level_0", "level_2"])
nal = nal.rename(columns={"level_1":"germline_call"})

# Get diet and age info
nada = pd.merge(eD.loc[:, ["Diet", "Age", "File"]], nal, on = ["File"], how = "inner", right_index = True)
# }}}

# We want to plot for each sample the number of alleles...
t = pd.pivot_table(nada, values = "germline_call", index = "File", columns = ["Diet", "Age"], 
               aggfunc = "count")
t = t.T.reset_index()
t = pd.melt(t, id_vars = ["Diet", "Age"], value_vars = ['S1_S1_001_PRIMER-S10_R1', 'S1_S1_001_PRIMER-S1_R1',
       'S1_S1_001_PRIMER-S2_R1', 'S1_S1_001_PRIMER-S3_R1','S1_S1_001_PRIMER-S4_R1', 'S1_S1_001_PRIMER-S5_R1',
       'S1_S1_001_PRIMER-S6_R1', 'S1_S1_001_PRIMER-S8_R1','S1_S1_001_PRIMER-S9_R1', 'S2_S2_001_PRIMER-S10_R1',
       'S2_S2_001_PRIMER-S1_R1', 'S2_S2_001_PRIMER-S2_R1','S2_S2_001_PRIMER-S3_R1', 'S2_S2_001_PRIMER-S4_R1',
       'S2_S2_001_PRIMER-S5_R1', 'S2_S2_001_PRIMER-S6_R1','S2_S2_001_PRIMER-S8_R1', 'S2_S2_001_PRIMER-S9_R1',
       'S3_S3_001_PRIMER-S10_R1', 'S3_S3_001_PRIMER-S2_R1','S3_S3_001_PRIMER-S3_R1', 'S3_S3_001_PRIMER-S4_R1',
       'S3_S3_001_PRIMER-S5_R1', 'S3_S3_001_PRIMER-S6_R1','S3_S3_001_PRIMER-S7_R1', 'S3_S3_001_PRIMER-S8_R1',
       'S3_S3_001_PRIMER-S9_R1', 'S4_S4_001_PRIMER-S10_R1','S4_S4_001_PRIMER-S1_R1', 'S4_S4_001_PRIMER-S3_R1',
       'S4_S4_001_PRIMER-S4_R1', 'S4_S4_001_PRIMER-S5_R1','S4_S4_001_PRIMER-S6_R1',
       'S4_S4_001_PRIMER-S8_R1', 'S4_S4_001_PRIMER-S9_R1','S5_S5_001_PRIMER-S10_R1', 'S5_S5_001_PRIMER-S1_R1',
       'S5_S5_001_PRIMER-S2_R1', 'S5_S5_001_PRIMER-S4_R1','S5_S5_001_PRIMER-S5_R1', 'S5_S5_001_PRIMER-S7_R1',
       'S5_S5_001_PRIMER-S8_R1', 'S5_S5_001_PRIMER-S9_R1']).dropna().reset_index()

# Add also files that had 0 novel alleles so we have accurate number of files studied
t = pd.concat([pd.DataFrame({"index":[0, 0, 0, 0, 0, 0],
              "Diet":["AL lifelong", "AL lifelong", "DR lifelong", "DR lifelong", "AL lifelong", "AL_DR20M"], 
              "Age":[24, 24, 24, 20, 24, 24], 
              "File":["S1_S1_001_PRIMER-S7_R1", "S3_S3_001_PRIMER-S1_R1", "S4_S4_001_PRIMER-S2_R1", "S5_S5_001_PRIMER-S3_R1", 
                      "S5_S5_001_PRIMER-S6_R1", "S4_S4_001_PRIMER-S7_R1"], 
              "value":[0, 0, 0, 0, 0, 0]}), t], ignore_index=True)

# {{{
fig, ax = plt.subplots(figsize = (8,4))
ax = sns.boxplot(data = t, x = "Age", y = "value", hue = "Diet", palette = palette2, 
           hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"], 
                showfliers = False)
ax = sns.swarmplot(x='Age', y='value', data=t, color=".25", hue = "Diet", dodge = True, 
                  hue_order = ["DR lifelong", "AL_DR16M", "AL_DR20M", "AL lifelong"])

ax.tick_params(axis = "x", labelsize = 18)
ax.tick_params(axis = "y", labelsize = 18)
ax.set_xlabel("Age (M)", fontsize = 20)
ax.set_ylabel("New alleles count", fontsize = 20)
plt.title("Number of new alleles per sample\n", fontsize = 22)
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles[:4], labels[:4], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 16)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for line in leg.get_lines():
    line.set_linewidth(3.5)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}count_novelAlleles_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

# {{{
# Check the diet effect, there are no real differences

list_ages = [5, 20, 24]

for d in list_ages:
    mod = ols('value ~ Diet', data = t[t["Age"] == d]).fit()
    aov_table = sm.stats.anova_lm(mod, typ=2)
        # If we want to get the degrees of freedom F statistic and sum of squares, print
        #print(aov_table)

        # Pairwise t-tests to get p-value
    pair_t = mod.t_test_pairwise('Diet', method = 'Bonferroni', alpha = 0.05)
    print("ANOVA Diet comparison from Age: {}".format(d))
    display(pd.DataFrame(pair_t.result_frame))

# }}}

# {{{
list_diets = ["AL_DR16M", "DR lifelong", "AL_DR20M", "AL lifelong"]

for d in list_diets:
    temp = t[t["Diet"] == d].copy()
    if d == "AL_DR16M":
        temp1 = t[t["Diet"] == "AL lifelong"]
        temp = pd.concat([temp, temp1[temp1["Age"] == 5], temp1[temp1["Age"] == 16]], 
                         ignore_index=True)
    elif d == "AL_DR20M":
        temp1 = t[t["Diet"] == "AL lifelong"]
        temp = pd.concat([temp, temp1[temp1["Age"] == 5], temp1[temp1["Age"] == 16],
                          temp1[temp1["Age"] == 20]], ignore_index=True)
    # Run statistics
    temp["Age"] = temp["Age"].astype(str)
    mod = ols('value ~ Age', data = temp).fit()
    aov_table = sm.stats.anova_lm(mod, typ=2)
    print("ANOVA Age comparison from Diet: {}".format(d))
    display(pd.DataFrame(aov_table))

    try:
        pair_t = mod.t_test_pairwise('Age', alpha = 0.05)
        display(pd.DataFrame(pair_t.result_frame))
    except:
        continue
# }}}

# #### It seems like the number of novel alleles per sample is not influenced by age or diet unfortunately  
#   
# #### It only shows that from young to super old AL there is a significant decrease

## Checking the evidence of assignment
# Read file
eal = pd.read_csv("../productive_vdj/novel_alleles/evidence_assignment.tsv", sep = " ")
# Get diet and age info
eal = pd.merge(eD.loc[:, ["Diet", "Age", "File"]], eal, on = ["File"], how = "inner", right_index = True)
# Long the before and after
eal_l = pd.melt(eal, id_vars = ["Diet", "Age", "File", "Call"], value_vars = ["Before", "After"]).dropna()

# {{{
fig, ax = plt.subplots(figsize = (8,4))
sns.pointplot(data = eal_l[eal_l["Call"] == "NotInGenotype"], x = "variable", 
              y = "value", hue = "Diet", palette = palette2)

ax.tick_params(axis = "x", labelsize = 18)
ax.tick_params(axis = "y", labelsize = 18)
ax.set_xlabel("")
ax.set_ylabel("Alleles not in genotype", fontsize = 20)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}Not_inGenotype_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

# {{{
fig, ax = plt.subplots(figsize = (8,4))
sns.pointplot(data = eal_l[eal_l["Call"] == "Ambiguous"], x = "variable", y = "value", hue = "Diet", palette = palette2)
ax.tick_params(axis = "x", labelsize = 18)
ax.tick_params(axis = "y", labelsize = 18)
ax.set_xlabel("")
ax.set_ylabel("Ambiguous genotype calls", fontsize = 20)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}AmbiguousCalls_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}
# ### Now we want to look for alleles that are shared between our mice


eal_l[eal_l["variable"] == "After"]["value"].mean()

eal_l[eal_l["variable"] == "Before"]["value"].mean()

0.2522965178380688*100

# {{{
# Get the count of alleles found in each file

matrix_nada = pd.merge(nada.groupby(["File", 'polymorphism_call']).size().unstack(fill_value=0), nada.loc[:, ["Diet", "Age", "File"]].drop_duplicates(), 
                       on = "File", how = "inner")
# sort the matrix according to age and diet
mat2 = matrix_nada.sort_values(by = ["Age", "Diet"])
mat = mat2.drop(columns = ["Age", "Diet"])
# Make a stable list of files so they keep the order
matl = mat["File"].unique().tolist()
mat = mat.set_index("File")
mat.index = pd.CategoricalIndex(mat.index, categories = matl)
mat.sort_index(level=0, inplace = True)
# }}}

# {{{
# Colors for the ages
batches2 = pd.Series(list(mat2["Age"]), index = list(mat2["File"]))
lut = dict(zip(batches2.unique(), ["lavenderblush", "pink", "hotpink", "magenta"]))
row_colors = batches2.map(lut)

# color per diet
batches3 = pd.Series(list(mat2["Diet"]), index = list(mat2["File"]))
lut2 = dict(zip(batches3.unique(), ["dodgerblue", "red", "teal", "gold"]))

row_colors2 = batches3.map(lut2)
g = sns.clustermap(mat, cmap="Blues", cbar_kws = {'orientation':'horizontal'},
            figsize = (16,35), col_cluster = False, row_cluster = False, row_colors=[row_colors, row_colors2])

g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")

matplotlib.rcParams['pdf.fonttype'] = 42

if run_type != "dry":
    plt.savefig("{}Heatmap_SharedNovel_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

# {{{
# Make legends in a separate file, we will put them in on inkscape
g = plt.figure()
gs2 = matplotlib.gridspec.GridSpec(2,5)
# Make Legend of Diets

ax2 = g.add_subplot(gs2[1, 0], frameon = False)
ax2.axis("off")
blue_patch = mpatches.Patch(color="#08306b", label="Novel allele")
ax2.legend(handles=[blue_patch], ncol =1, 
           loc='upper right')

# Make legend of Ages
#sub3 = fig.add_subplot(255, frameon = False) # {20: 'tab:pink', 5: 'tab:gray', 16: 'tab:olive', 24: 'tab:cyan'}
ax3 = g.add_subplot(gs2[1,2], frameon = False)
ax3.axis("off")
blue_patch = mpatches.Patch(color="lavenderblush", label="5M")
orange_patch = mpatches.Patch(color="pink", label="16M")
green_patch = mpatches.Patch(color="hotpink", label="20M")
red_patch = mpatches.Patch(color="magenta", label="24M")
ax3.legend(handles=[blue_patch, orange_patch, green_patch, red_patch], ncol = 4, 
           loc='center')

ax4 = g.add_subplot(gs2[1,3], frameon = False)
ax4.axis("off")
blue_patch = mpatches.Patch(color="red", label="DR lifelong")
orange_patch = mpatches.Patch(color="teal", label="AL_DR16M")
green_patch = mpatches.Patch(color="gold", label="AL_DR20M")
red_patch = mpatches.Patch(color="dodgerblue", label="AL lifelong")
ax4.legend(handles=[blue_patch, orange_patch, green_patch, red_patch], ncol = 4, 
           loc='lower right')
plt.tight_layout()

matplotlib.rcParams['pdf.fonttype'] = 42

if run_type != "dry":
    plt.savefig("{}Heatmap_labels_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

# {{{
# Plot counts of per allele
ca = mat.T
ca_df = pd.DataFrame(test.sum(axis = 1))

fig, ax = plt.subplots(figsize = (12, 3))
g = sns.barplot(ca_df.index, ca_df[0], color = "#08306b")

g.set(xticklabels=[])
#ax.tick_params(axis = "x", labelsize = 0)
ax.tick_params(axis = "y", labelsize = 18)
ax.set_xlabel("")
ax.set_ylabel("Number of mice", fontsize = 20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

matplotlib.rcParams['pdf.fonttype'] = 42
plt.tight_layout()

if run_type != "dry":
    plt.savefig("{}MiceNum_perNovelAllele_{}.pdf".format(path, new_day))
else:
    plt.show()
# }}}

# ## Checking mice that had polymorphisms in the ileum if they have them in the spleen

# {{{
# Mouse SGRO-1137 (I6_S5 ileum)
nada[nada["File"] == "S1_S1_001_PRIMER-S5_R1"]
# In Ileum we also had the IGHV5-9*04_A163T_C164A_G172A_T222C_A251G

# But ileum also had two more: IGHV1-5*01_A73G_A87C_A159T_G162C_C249T and IGHV1-82*01_C84T
# }}}

# {{{
# Mouse SGRO-0975 (I8_S8 ileum)
nada[nada["File"] == "S2_S2_001_PRIMER-S8_R1"]

# In ileum we had two alleles, but none are found here: IGHV1-5*01_A73G_A87C_A159T_G162C_C249T and IGHV1-69*01_T226G
# }}}




