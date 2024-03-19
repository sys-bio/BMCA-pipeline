import sys

from . import antemll 



try: 
    os.mkdir(output_dir)
except FileExistsError:
    pass

ant = None
#'../../../data/interim/Antimony/JSexample22_reg1.ant'  
data_file = None
# '../../../data/interim/generated_data/JSexample22-reg1/JSexample22_reg1_1.01.csv'
output_dir = None # sys.argv[3] + "/output/"
# '../../../data/interim/generated_data/JSexample22-reg1'
crispri = None # 0 or 1 for yes or no

#################################################
#################################################

import tellurium as te
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import aesara.tensor as at
import aesara
floatX = aesara.config.floatX
import pymc as pm
import numpy as np

import cobra
import os

from scipy import stats
import scipy as sp

os.chdir('../../..')
from src import util
import emll
from emll.aesara_utils import LeastSquaresSolve
os.chdir('notebooks/topologyB/all_data/')

###
"""
For crispri, 
    A: 0.1x
    B: 0.2x
    C: 0.3x
    D: 0.4x
    E: 0.5x
    F: 1.01x

For CRISPRa, 
    A: 1.01x
    B: 1.5x
    C: 3x
    D: 5x
    E: 7x
    F: 10x

"""

#########################################################
##########################################################

def runBMCA(ant, data_file, output_dir): 

    r = te.loada(ant)
    r.conservedMoietyAnalysis = True
    r.steadyState()

        data_file_A = data_file + '_0.1.csv'
        data_file_B = data_file + '_0.2.csv'
        data_file_C = data_file + '_0.3.csv'
        data_file_D = data_file + '_0.4.csv'
        data_file_E = data_file + '_0.5.csv'
        data_file_F = data_file + '_1.01.csv'
        data_file_G = data_file + '_1.5.csv'
        data_file_H = data_file + '_3.csv'
        data_file_I = data_file + '_5.csv'
        data_file_J = data_file + '_7.csv'
        data_file_K = data_file + '_10.csv'
        pt_labels = ['0.1x', '0.2x', '0.3x', '0.4x', '0.5x','1.01x', 
                     '1.5x', '3x', '5x', '7x', '10x']

    BMCA_obj_A = antemll.BMCA(ant, data_file_A)
    BMCA_obj_B = antemll.BMCA(ant, data_file_B)
    BMCA_obj_C = antemll.BMCA(ant, data_file_C)
    BMCA_obj_D = antemll.BMCA(ant, data_file_D)
    BMCA_obj_E = antemll.BMCA(ant, data_file_E)
    BMCA_obj_F = antemll.BMCA(ant, data_file_F)
    BMCA_obj_G = antemll.BMCA(ant, data_file_G)
    BMCA_obj_H = antemll.BMCA(ant, data_file_H)
    BMCA_obj_I = antemll.BMCA(ant, data_file_I)
    BMCA_obj_J = antemll.BMCA(ant, data_file_J)
    BMCA_obj_K = antemll.BMCA(ant, data_file_K)

    ### Run BMCA

    traceA = util.run_ADVI(BMCA_obj_A)
    traceB = util.run_ADVI(BMCA_obj_B)
    traceC = util.run_ADVI(BMCA_obj_C)
    traceD = util.run_ADVI(BMCA_obj_D)
    traceE = util.run_ADVI(BMCA_obj_E)
    traceF = util.run_ADVI(BMCA_obj_F)
    traceG = util.run_ADVI(BMCA_obj_G)
    traceH = util.run_ADVI(BMCA_obj_H)
    traceI = util.run_ADVI(BMCA_obj_I)
    traceJ = util.run_ADVI(BMCA_obj_J)
    traceK = util.run_ADVI(BMCA_obj_K)

    ExTrace_A = (traceA['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_B = (traceB['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_C = (traceC['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_D = (traceD['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_E = (traceE['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_F = (traceF['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_G = (traceG['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_H = (traceH['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_I = (traceI['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_J = (traceJ['posterior']['Ex']).to_numpy().squeeze()
    ExTrace_K = (traceK['posterior']['Ex']).to_numpy().squeeze()

    medEx_A = np.median(ExTrace_A, axis=0)
    medEx_B = np.median(ExTrace_B, axis=0)
    medEx_C = np.median(ExTrace_C, axis=0)
    medEx_D = np.median(ExTrace_D, axis=0)
    medEx_E = np.median(ExTrace_E, axis=0)
    medEx_F = np.median(ExTrace_F, axis=0)
    medEx_G = np.median(ExTrace_G, axis=0)
    medEx_H = np.median(ExTrace_H, axis=0)
    medEx_I = np.median(ExTrace_I, axis=0)
    medEx_J = np.median(ExTrace_J, axis=0)
    medEx_K = np.median(ExTrace_K, axis=0)

    medExs = [medEx_A, medEx_B, medEx_C, medEx_D, medEx_E, medEx_F, 
              medEx_G, medEx_H, medEx_I, medEx_J, medEx_K]

    ### Compare elasticity results with ground truth values

    gtE = pd.DataFrame(r.getScaledElasticityMatrix(), index=r.getReactionIds(), columns=r.getFloatingSpeciesIds())
    gtE['pt'] = ['gt']* len(gtE)
    gtE.index.name = 'reactions'
    gtE.reset_index(inplace=True)
    gtE.set_index(['pt', 'reactions'], inplace=True)
    medEx_pile = [gtE]

    for i, lvl in enumerate(pt_labels):
        mdEx = pd.DataFrame(medExs[i], index=r.getReactionIds(), columns=r.getFloatingSpeciesIds())  
        mdEx_by_pt = pd.concat([mdEx], keys=[lvl], names=['perturbation'])
        medEx_pile.append(mdEx_by_pt)

    medEx_df = pd.concat(medEx_pile)
    medEx_df.to_csv('kataraTEST.csv')

    ### CALCULATING FCCs
    gtFCC = pd.DataFrame(r.getScaledFluxControlCoefficientMatrix(), index=r.getReactionIds(), columns=r.getReactionIds())

    # we are trying to add the ground truth to the larger
    # df of FCC predictions
    gtFCC.index.name = 'pt_rxn'
    gtFCC.reset_index(inplace=True)
    gtFCC['pt']=['gt']*len(gtFCC)
    gtFCC.set_index(['pt', 'pt_rxn'], inplace=True)

    postFCCs = []

    for i, ExTrace in enumerate(ExTraces):
        postFCCs.append(util.estimate_CCs(BMCA_obj[i], ExTrace))

    
    postFCC1 = util.estimate_CCs(BMCA_obj_A, ExTrace_A)
    postFCC15 = util.estimate_CCs(BMCA_obj_B, ExTrace_B)
    postFCC3 = util.estimate_CCs(BMCA_obj_C, ExTrace_C)
    postFCC5 = util.estimate_CCs(BMCA_obj_D, ExTrace_D)
    postFCC7 = util.estimate_CCs(BMCA_obj_E, ExTrace_E)
    postFCC10 = util.estimate_CCs(BMCA_obj_F, ExTrace_F)

    postFCCs = [postFCC1, postFCC15, postFCC3, postFCC5, postFCC7, postFCC10]

    postFCCdfs = pd.concat([util.append_FCC_df(postFCCs[i], pt_labels[i], r) for i in range(len(postFCCs))])
    prdFCCs = pd.pivot_table(postFCCdfs, index=['pt_str','pt_rxn'], aggfunc='median', sort=False)
    prdFCCmeds = pd.concat([gtFCC, prdFCCs])

    ## Evaluating FCC ranking
        
    gtFCC=pd.DataFrame(r.getScaledFluxControlCoefficientMatrix(), columns=r.getReactionIds(), index=r.getReactionIds()).abs()
    scores = []
    for pt_level in postFCCs:
        postFCC_med=pd.DataFrame(np.median(pt_level, axis=0), columns=r.getReactionIds(), index=r.getReactionIds()).abs()
        # m1 = gtFCC.index.values[:, None] == gtFCC.columns.values
        postFCC_med = pd.DataFrame(np.select([m1], [float('Nan')], postFCC_med), columns=gtFCC.columns, index=gtFCC.index)
        postFCC_med_rankings= postFCC_med.rank(axis=1, ascending=False, na_option='keep')

        m1 = postFCC_med_rankings.isin([1.0])  
        m2 = postFCC_med_rankings.isin([2.0])  
        m3 = postFCC_med_rankings.isin([3.0])  
        a = m1.mul(r.getReactionIds()).apply(lambda x: [i for i in x if i], axis=1)
        b = m2.mul(r.getReactionIds()).apply(lambda x: [i for i in x if i], axis=1)
        c = m3.mul(r.getReactionIds()).apply(lambda x: [i for i in x if i], axis=1)

        prdRanks = pd.concat([a,b,c], axis=1)
        prdRanks['topThree'] = prdRanks[0] + prdRanks[1] + prdRanks[2]

        scores.append([len([i for i in prdRanks['topThree'][rxn] if i in trueRanks['topThree'][rxn]]) for rxn in r.getReactionIds()])

    topThreeCheckdf = pd.DataFrame(scores, columns=r.getReactionIds(), index=pt_labels).T
    topThreeCheckdf.to_csv(output_dir + 'top3breakdown.csv')
    # topThreeCheckdf.style.background_gradient(cmap='RdYlBu', axis=None)

    summary = (topThreeCheckdf.sum(axis=0)/(len(r.getReactionIds())*3)).round(3)
    summary.to_csv(output_dir + 'top3summary.csv')
    ### SAVE this into output file

