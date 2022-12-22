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
import numpy as np

path = "../../analysis/created_germlines/reannotated_db/"
# }}}

# {{{
## Formatting to get isotypes 

df = pd.read_csv(path + "merged_changeo_SPL.tsv", sep = "\t")
df = df.assign(sum_iso=df["isotype"].str[:3])
df.to_csv(path + "merged_changeo_SPL_isosum.tsv", sep = "\t", index=False)
# }}}

## Same with other tissue
df = pd.DataFrame()
df = pd.read_csv(path + "merged_changeo_ILE.tsv", sep = "\t")
df = df.assign(sum_iso=df["isotype"].str[:3])
df.to_csv(path + "merged_changeo_ILE_isosum.tsv", sep = "\t", index=False)
