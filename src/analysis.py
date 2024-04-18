import tellurium as te
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import os

print(os.getcwd())

######### USER-INPUT SECTION
MODEL_FILE_PATH = 'data/interim/Antimony/Simplified_Teusink_yeast.ant'

PT_LVL = '3' # must be str


######### END OF USER-INPUT SECTION

model = te.loada(MODEL_FILE_PATH)

# Establish labels for metabolite and reaction names
m_labels = [m for m in model.getFloatingSpeciesIds()]
r_labels = [r for r in model.getReactionIds()]

ex_labels = np.array([['$\epsilon_{' + '{0},{1}'.format(rlabel, mlabel) + '}$'
                    for mlabel in m_labels] for rlabel in r_labels]).flatten()

ex_file_labels = np.array([['E_' + '{0},{1}'.format(rlabel, mlabel)
                    for mlabel in m_labels] for rlabel in r_labels]).flatten()

# load the predicted Ex dataset
# then reshape into 1000 rows, 176 columns
# plot all data in each column as a scatter plot
pt3_dfs = []
for i in range(10):
    # df3 = pd.read_csv(f'data/interim/generated_data/simplTeusink-noReg_3x_{i}/output-allData/{PT_LVL}x_PredictedExs.csv', index_col=0).values
    df3 = pd.read_csv(f'data/interim/generated_data/simplTeusink-noReg_3x_{i}/output-allData/0_PredictedExs.csv', index_col=0).values
    pt3_dfs.append(pd.DataFrame(df3.reshape(1000,176), columns=ex_labels))

# make a plotting method
def plot10distr(colNo):
    only_iter = []
    for df in pt3_dfs:
        only_iter.append(df.loc[:,ex_labels[colNo]])
        # sns.swarmplot(data=df, x=ex_labels[colNo], size=3, alpha=0.8)
    only_iter_df = pd.concat(only_iter, axis=1)
    only_iter_df.columns = range(10)
    
    # sns.boxplot(data=only_iter_df)
    sns.violinplot(data=only_iter_df, palette = "Spectral", inner_kws=dict(box_width=15, whis_width=2, color=".8"))
    plt.xlabel('iterations')
    plt.ylabel('elasticity coefficient value')
    plt.title(f'10 iterations of BMCA runs for {ex_labels[colNo]}\n at perturbation {PT_LVL}x')
    plt.grid(axis='y')
    plt.savefig(f'data/results/{ex_file_labels[colNo]}')
    plt.close()
    
for i in range(len(ex_labels)):
    plot10distr(i)

# examples of dimension wrangling
"""a = np.linspace(0,175, 176).reshape((1,16,11))
b = np.tile(a,(1000,1))
c = b.reshape(1000,16,11)

c[0,:,:]

dim = c.shape
c.reshape(1000, 176)"""