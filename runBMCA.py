from src import antemll 

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

from src import util

##########################################################
##########################################################

def runBMCA(ant, data_file, output_dir, n_iter=45000, n_samp=3, omit=None): 

    if omit:
        output_dir = output_dir + f'/output-omit{omit}/'
    else: 
        output_dir = output_dir + '/output-allData/'
    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass

    r = te.loada(ant)
    r.conservedMoietyAnalysis = True
    r.steadyState()

    enzymes = ['e_' + i for i in r.getReactionIds()]
    internal = r.getFloatingSpeciesIds()
    external = r.getBoundarySpeciesIds()
    fluxes = ['v_' + i for i in r.getReactionIds()]

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

    data_files = [data_file_A, data_file_B, data_file_C, data_file_D, 
                  data_file_E, data_file_F, data_file_G, data_file_H, 
                  data_file_I, data_file_J, data_file_K]
    
    data = []
    if omit is None: 
        data = data_files
    if omit == 'fluxes':
        for file in data_files:
            data.append(pd.read_csv(file)[enzymes+internal+external])
        v_star = pd.read_csv(file)[fluxes].iloc[0].values
    elif omit == 'enzymes':
        for file in data_files:
            data.append(pd.read_csv(file)[fluxes+internal+external])
    elif omit == 'internal':
        for file in data_files:
            data.append(pd.read_csv(file)[fluxes+enzymes+external])
    elif omit == 'external':
        for file in data_files:
            data.append(pd.read_csv(file)[fluxes+enzymes+internal])
        
    BMCA_objs = []
    for i in data_files:
        BMCA_objs.append(antemll.antemll(ant, i))

    ### Run BMCA
    traces = []
    if omit is None:
        for i in BMCA_objs:
            traces.append(util.run_ADVI(i, output_dir, n_iter, n_samp=n_samp))

    elif omit == 'fluxes':
        for i in BMCA_objs:
            traces.append(util.runBayesInf_fluxes(i, r, data[0], output_dir, 
                                                  n_iter, n_samp=n_samp))
    elif omit == 'enzymes':
        for i in BMCA_objs:
            traces.append(util.runBayesInf_enzymes(i, r, data[0], output_dir, 
                                                  n_iter, n_samp=n_samp))
    elif omit == 'internal':
        for i in BMCA_objs:
            traces.append(util.runBayesInf_internal(i, r, data[0], output_dir, 
                                                  n_iter, n_samp=n_samp))
    elif omit == 'external':
        for i in BMCA_objs:
            traces.append(util.runBayesInf_external(i, r, data[0], output_dir, 
                                                  n_iter, n_samp=n_samp))

    ExTraces = []
    if n_samp == 1:
        for i in traces:
            ExTraces.append((i['posterior']['Ex']).to_numpy().squeeze())
    elif n_samp > 1: 
        for i in traces:
            Ex_samples = []
            for ii in range(n_samp):
                Ex_samples.append((i[ii]['posterior']['Ex']).to_numpy().squeeze())                
            trace = np.concatenate(Ex_samples)
            ExTraces.append(trace)
    
    medExs =[]
    for i in ExTraces:
        medExs.append(np.median(i, axis=0))

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
    medEx_df.to_csv(output_dir + 'medianPredictedExs.csv')

    
    ### CALCULATING FCCs
    gtFCC = pd.DataFrame(r.getScaledFluxControlCoefficientMatrix(), index=r.getReactionIds(), columns=r.getReactionIds())

    # we are trying to add the ground truth to the larger
    # df of FCC predictions
    gtFCC.index.name = 'pt_rxn'
    gtFCC.reset_index(inplace=True)
    gtFCC['pt']=['gt']*len(gtFCC)
    gtFCC.set_index(['pt', 'pt_rxn'], inplace=True)
    
    postFCCs = []
    
    if omit == 'enzymes':
        EtTraces = []
        et_samples = []
        for i in traces:
            for ii in range(n_samp):
                et_samples.append((i[ii]['posterior']['e_t']).to_numpy().squeeze())
            trace = np.concatenate(et_samples)
            EtTraces.append(trace)
        medEts =[]
        for i in EtTraces:
            medEts.append(np.median(i, axis=0).transpose())
        
        for i, med_et in enumerate(medEts):
            BMCA_objs[i].vn[BMCA_objs[i].vn == 0] = 1e-6
            a = np.diag(med_et / BMCA_objs[i].vn.values)    
            postFCCs.append(util.estimate_CCs(BMCA_objs[i], ExTraces[i], n_samp, a))

    elif omit == 'fluxes':
        vtTraces = []
        vt_samples = []
        for i in traces:
            for ii in range(n_samp):
                vt_samples.append((i[ii]['posterior']['v_t']).to_numpy().squeeze())
            trace = np.concatenate(vt_samples)
            vtTraces.append(trace)
        medvts =[]
        for i in vtTraces:
            medvts.append(np.median(i, axis=0).transpose())
        
        for i, med_vt in enumerate(medvts):
            BMCA_objs[i].vn[BMCA_objs[i].vn == 0] = 1e-6
            a = np.diag(BMCA_objs[i].en.values / med_vt)# BMCA_obj.vn.values)
            postFCCs.append(util.estimate_CCs(BMCA_objs[i], ExTraces[i], n_samp, a))

    else:
        print('else')
        for i, BMCA_obj in enumerate(BMCA_objs):
            BMCA_obj.vn[BMCA_obj.vn == 0] = 1e-6
            a = np.diag(BMCA_obj.en.values / BMCA_obj.vn.values)
            postFCCs.append(util.estimate_CCs(BMCA_obj, ExTraces[i], n_samp, a))
    
    postFCCdfs = pd.concat([util.append_FCC_df(postFCCs[i], pt_labels[i], r) for i in range(len(postFCCs))])
    prdFCCs = pd.pivot_table(postFCCdfs, index=['pt_str','pt_rxn'], aggfunc='median', sort=False)
    prdFCCs.to_csv(output_dir + 'predictedFCCs.csv')   
    prdFCCmeds = pd.concat([gtFCC, prdFCCs])
    prdFCCmeds.to_csv(output_dir + 'predictedFCCmedians.csv')   
  
    ## Evaluating FCC ranking
        
    gtFCC=pd.DataFrame(r.getScaledFluxControlCoefficientMatrix(), columns=r.getReactionIds(), index=r.getReactionIds()).abs()
    
    m1 = gtFCC.index.values[:, None] == gtFCC.columns.values
    gtFCC = pd.DataFrame(np.select([m1], [float('Nan')], gtFCC), columns=gtFCC.columns, index=gtFCC.index)
    gtFCC_rankings= gtFCC.rank(axis=1, ascending=False, na_option='keep')
    m1 = gtFCC_rankings.isin([1.0])  
    m2 = gtFCC_rankings.isin([2.0])  
    m3 = gtFCC_rankings.isin([3.0])  
    a = m1.mul(r.getReactionIds()).apply(lambda x: [i for i in x if i], axis=1)
    b = m2.mul(r.getReactionIds()).apply(lambda x: [i for i in x if i], axis=1)
    c = m3.mul(r.getReactionIds()).apply(lambda x: [i for i in x if i], axis=1)

    trueRanks = pd.concat([a,b,c], axis=1)
    trueRanks['topThree'] = trueRanks[0] + trueRanks[1] + trueRanks[2]
    
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

