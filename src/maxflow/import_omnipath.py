#%%


import omnipath as op
import pandas as pd



omni_interactions = op.interactions.OmniPath.get()

omni_interactions.to_csv('/Users/charliebarker/Desktop/Melanoma_Resistance/data/omnipath.csv')

# %%
