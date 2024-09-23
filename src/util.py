# from PSJ's emll/py

import numpy as np
import scipy as sp
from scipy.stats import spearmanr

import tellurium as te
import cobra
import aesara
import aesara.tensor as at
import csv
import numpy as np
import libsbml
import os
import random
import logging
logging.getLogger("cobra").setLevel(logging.ERROR)

import pandas as pd

from emll.aesara_utils import LeastSquaresSolve

import aesara.tensor as pt
import pymc as pm

# plotting
import matplotlib.pyplot as plt
import seaborn as sns
import arviz as az

def generate_data(model_file, perturbation_levels, data_folder, concurrent=1):
    """
    Takes in SBML model_file and creates perturbation data files. 
    Deposits perturbation data files into data_folder

    data_folder='data\interim\generated_data'
    """
    model_name = model_file.split('.ant')[0]
    model_name = model_name.split('/')[-1]
    
    for pl in perturbation_levels:
        
        r = te.loada(model_file)
        r.conservedMoietyAnalysis = True
        
        exMet = r.getBoundarySpeciesIds()
        inMet = r.getFloatingSpeciesIds()
        fluxIDs = ['v_' + i for i in r.getReactionIds()]
        e_list = [i for i in r.getGlobalParameterIds() if 'e_' in i]   
        
        pertLevel = pl #/100 
        perturbation_level = [pertLevel]# [1 - pertLevel, 1 + pertLevel]
        # get ride of 0 values
        # here

        header = e_list + exMet + inMet + fluxIDs        
        cont = True
        
        with open(data_folder + f'{model_name}_{pl}.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
    
            try: # base case    
                spConc = list(r.simulate(0,1000000)[-1])[1:]
                # r.conservedMoietyAnalysis = True
                # r.steadyState()

                enzymes = [r.getValue(e) for e in e_list]
                exMet_values = [r.getValue(m) for m in exMet]
                fluxes = list(r.getReactionRates())

                writer.writerow(enzymes + exMet_values + spConc + fluxes)
            except:
                print('failed at', params, level) # pass #    
                cont=False
            
            # perturbed enzyme cases
            if cont:
                if concurrent==1:
                    for params in e_list:
                        for level in perturbation_level:
                            try: 
                                r.resetToOrigin()
                                r.setValue(params, level*r.getValue(params))
                                
                                spConc = list(r.simulate(0,1000000)[-1])[1:]
                                # r.steadyState()
                                enzymes = [r.getValue(e) for e in e_list]
                                exMet_values = [r.getValue(m) for m in exMet]
                                fluxes = list(r.getReactionRates())
                                
                                writer.writerow(enzymes + exMet_values + spConc + fluxes)
                            except:
                                pass
                

def ant_to_cobra(antimony_file):
    """
    This method takes in an Antimony file and converts it to a Cobra-
    friendly format by removing all boundary species and replacing all
    rate laws with a constant. Both an Antimony file and an SBML are produced.
    
    """
    output_name = antimony_file.split('/')[-1]
    output_name = output_name.split('.')[0]

    # load normal Antimony file
    with open(antimony_file,'r') as file:
        lines = file.readlines()
    section_indices = []
    c_index = -1
    for i, line in enumerate(lines): 
        if '//' in line:
            section_indices.append(i)
        if '// Compartment' in line and c_index == -1: 
            c_index = i

    next_section = section_indices.index(c_index) + 1
    with open(antimony_file,'r') as file:
        lines = file.readlines()[c_index: section_indices[next_section]]
    
    for line in lines:
        line = line.strip()
        if '$' not in line: 
            with open(f'{output_name}_cobra.ant', 'a') as f:
                f.write(line + '\n')

        else: 
            no_bd_sp = [i for i in line.split(',') if '$' not in i]
            no_bd_sp = [i for i in no_bd_sp if i!='']

            line = ','.join(no_bd_sp)
            if line != '' and line[-1] != ';': 
                line += ';'
            with open(f'{output_name}_cobra.ant', 'a') as f:
                f.write(line + '\n')
    
    r = te.loada(antimony_file)
    doc = libsbml.readSBMLFromString(r.getSBML())
    model = doc.getModel()

    bd_sp = r.getBoundarySpeciesIds()

    reactants_list=[]
    products_list=[]

    for n in range(len(r.getReactionIds())): 
        rxn = model.getReaction(n)
        reactants = []
        products = []

        for reactant in range(rxn.getNumReactants()):   
            stoich = rxn.getReactant(reactant).getStoichiometry()
            if stoich == 1: 
                reactants.append(rxn.getReactant(reactant).species)    
            else:
                reactants.append(str(stoich) + ' ' + rxn.getReactant(reactant).species)
        reactants_list.append([i for i in reactants if i not in bd_sp])
        
        for product in range(rxn.getNumProducts()):
            stoich = rxn.getProduct(product).getStoichiometry()
            if stoich == 1: 
                products.append(rxn.getProduct(product).species)    
            else:
                products.append(str(stoich) + ' ' + rxn.getProduct(product).species)
        products_list.append([i for i in products if i not in bd_sp])
        
    for i in range(len(reactants_list)):
        r1 = ' + '.join(reactants_list[i])
        p1 = ' + '.join(products_list[i])
        with open(f'{output_name}_cobra.ant', 'a') as f:
            f.write(r.getReactionIds()[i]+ ': ' + r1 + ' -> ' + p1 + '; 1;\n')

    with open(f'{output_name}_cobra.ant', 'a') as f:
        f.write('\n')

    for sp in r.getFloatingSpeciesIds():
        with open(f'{output_name}_cobra.ant', 'a') as f:
            f.write(sp + ' = 1;\n')



def initialize_elasticity(ela_matrix, name=None, b=0.01, alpha=5, sd=1,
                          m_compartments=None, r_compartments=None):
    """ Initialize the elasticity matrix, adjusting priors to account for
    reaction stoichiometry. Uses `SkewNormal(mu=0, sd=sd, alpha=sign*alpha)`
    for reactions in which a metabolite participates, and a `Laplace(mu=0,
    b=b)` for off-target regulation. 

    Also accepts compartments for metabolites and reactions. If given,
    metabolites are only given regulatory priors if they come from the same
    compartment as the reaction.
    
    Parameters
    ==========

    ela_matrix : np.ndarray
        A (nr x nm) elasticity matrix for the given reactions and metabolites
    name : string
        A name to be used for the returned pymc probabilities
    b : float
        Hyperprior to use for the Laplace distributions on regulatory interactions
    alpha : float
        Hyperprior to use for the SkewNormal distributions. As alpha ->
        infinity, these priors begin to resemble half-normal distributions.
    sd : float
        Scale parameter for the SkewNormal distribution.
    m_compartments : list
        Compartments of metabolites. If None, use a densely connected
        regulatory prior.
    r_compartments : list
        Compartments of reactions

    Returns
    =======

    E : pymc matrix
        constructed elasticity matrix

    """
    
    if name is None:
        name = 'ex'

    if m_compartments is not None:
        assert r_compartments is not None, \
            "reaction and metabolite compartments must both be given"

        regulation_array = np.array(
            [[a in b for a in m_compartments]
              for b in r_compartments]).flatten()
        
    else:
        # If compartment information is not given, assume all metabolites and
        # reactions are in the same compartment
        regulation_array = np.array([True] * (ela_matrix.shape[1] * ela_matrix.shape[0]))


    # Find where the guessed E matrix has zero entries
    e_flat = ela_matrix.flatten()
    nonzero_inds = np.where(e_flat != 0)[0]
    offtarget_inds = np.where(e_flat == 0)[0]
    e_sign = np.sign(e_flat[nonzero_inds])

    # For the zero entries, determine whether regulation is feasible based on
    # the compartment comparison
    offtarget_reg = regulation_array[offtarget_inds]
    reg_inds = offtarget_inds[offtarget_reg]
    zero_inds = offtarget_inds[~offtarget_reg]

    num_nonzero = len(nonzero_inds)
    num_regulations = len(reg_inds)
    num_zeros = len(zero_inds)
    
    # Get an index vector that 'unrolls' a stacked [kinetic, capacity, zero]
    # vector into the correct order
    flat_indexer = np.hstack([nonzero_inds, reg_inds, zero_inds]).argsort()
        
    if alpha is not None:
        e_kin_entries = pm.SkewNormal(
            name + '_kinetic_entries', sigma=sd, alpha=alpha, shape=num_nonzero,
            initval= 0.1 + np.abs(np.random.randn(num_nonzero)))
    else:
        e_kin_entries = pm.HalfNormal(
            name + '_kinetic_entries', sigma=sd, shape=num_nonzero,
            initval= 0.1 + np.abs(np.random.randn(num_nonzero)))
    
    e_cap_entries = pm.Laplace(
        name + '_capacity_entries', mu=0, b=b, shape=num_regulations,
        initval=b * np.random.randn(num_regulations))
    
    flat_e_entries = pt.concatenate(
        [e_kin_entries * e_sign,  # kinetic entries
         e_cap_entries,           # capacity entries
         pt.zeros(num_zeros)])     # different compartments
        
    E = flat_e_entries[flat_indexer].reshape(ela_matrix.shape)
    
    return E

def elasticity_to_CCC(BMCA, scaledE=None):

    if scaledE is None:
        scaledE = BMCA.Ex

    r = te.loada(BMCA.model_file)
    r.conservedMoietyAnalysis = True
    link_matrix = r.getLinkMatrix()
    Nr = r.getReducedStoichiometryMatrix()

    ##### this line needs to be workshopped
    unscaledE = np.linalg.inv(np.diag(1/BMCA.v_star)) @ scaledE @ np.linalg.inv(np.diag(BMCA.x_star))

    invJac = np.linalg.inv(-Nr@unscaledE@link_matrix)
    idMat = np.identity(len(BMCA.v_star))

    # unscaled concentration and flux control coefficients, respectively
    Cx = link_matrix@invJac@Nr 
    CJ = np.matmul(unscaledE, Cx) + idMat # unscaled FCC

    # scaled concentration and flux control coefficients, respectively  
    ##### these two lines need to be workshopped
    CxS = np.diag(1/BMCA.x_star) @ Cx @ np.diag(BMCA.v_star)
    CJS = np.diag(1/BMCA.v_star) @ CJ @ np.diag(BMCA.v_star)

    return CxS, CJS

def estimate_CCs(BMCA_obj, Ex):
    BMCA_obj.vn[BMCA_obj.vn == 0] = 1e-6
    
    a = np.diag(BMCA_obj.en.values / BMCA_obj.vn.values)
    a = np.diag(a)
    a = a[np.newaxis,:].repeat(1000, axis=0)

    Ex_ss = a @ Ex
    As = BMCA_obj.N @ np.diag(BMCA_obj.v_star) @ Ex_ss
    bs = BMCA_obj.N @ np.diag(BMCA_obj.v_star)
    bs = bs[np.newaxis, :].repeat(1000, axis=0)
    
    As = at.as_tensor_variable(As)
    bs = at.as_tensor_variable(bs)

    def solve_aesara(A, b):
        os.chdir('..')
        from emll.aesara_utils import LeastSquaresSolve
        os.chdir('src')
        rsolve_op = LeastSquaresSolve()
        return rsolve_op(A, b).squeeze()

    CCC, _ = aesara.scan(lambda A, b: solve_aesara(A, b),
                        sequences=[As, bs], strict=True)

    identity = np.eye(len(BMCA_obj.N.T))
    identity = identity[np.newaxis,:].repeat(1000, axis=0)
    
    FCC = (Ex_ss @ CCC.eval()) + identity
    
    return CCC.eval(), FCC

def calculate_e_hat(BMCA_obj, v_hat_obs, x_terms, y_terms): 
    one_n = np.ones([len(x_terms.eval()),len(BMCA_obj.en)])
    product = (v_hat_obs * (one_n + x_terms + y_terms)).eval()
    product[product == 0 ] = 1E-6

    return aesara.tensor.reciprocal(product)

def plot_elbo(approx, output_dir, n_iter, x_start=0):
    with sns.plotting_context('notebook', font_scale=1.2):

        fig = plt.figure(figsize=(5,4))
        plt.plot(approx.hist + 30, '.', rasterized=True, ms=1)
        # plt.ylim([-1E1, 1E3])
        plt.xlim([x_start, n_iter])
        sns.despine(trim=True, offset=10)

        plt.ylabel('-ELBO')
        plt.xlabel('Iteration')
        plt.title('in vitro ADVI convergence')
        plt.tight_layout()

        plt.savefig(output_dir + 'elbo_convergence.svg', transparent=True, dpi=200)

def sample_draw(approx, n_samp):
    if n_samp > 1:
        samples = []
        for i in range(n_samp): 
            samples.append(approx.sample(draws=1000, random_seed=i))
        return samples
    else:
        return approx.sample(draws=1000, random_seed=n_samp)

def run_ADVI(BMCA_obj, output_dir, n_iter, n_samp=1):
    with pm.Model() as pymc_model:
        
        # Initialize elasticities
        Ex_t = pm.Deterministic('Ex', initialize_elasticity(BMCA_obj.Ex.to_numpy(), name='Ex'))
        Ey_t = pm.Deterministic('Ey', initialize_elasticity(BMCA_obj.Ey.to_numpy(), name='Ey'))
        e_obs = pm.Normal('e_obs', mu=1, sigma=1, observed=BMCA_obj.en.T)
        chi_obs = pm.Normal('chi_obs', mu=0, sigma=10, observed=BMCA_obj.xn.T)
        y_obs = pm.Normal('y_obs', mu=0, sigma=10, observed=BMCA_obj.yn.T)
        likelihood = pm.Deterministic('vn', e_obs * (np.ones(BMCA_obj.en.T.shape) + pm.math.dot(Ex_t,chi_obs) + pm.math.dot(Ey_t,y_obs)))
        v_hat_obs = pm.Normal('v_hat_obs', mu=likelihood, sigma=0.1, observed=BMCA_obj.vn.squeeze().T)
    
    with pymc_model:
        advi = pm.ADVI()
        tracker = pm.callbacks.Tracker(
            mean = advi.approx.mean.eval,
            std = advi.approx.std.eval
        )
        approx = advi.fit(
            n=n_iter, 
            callbacks = [tracker],
            obj_optimizer=pm.adagrad_window(learning_rate=5E-3), 
            total_grad_norm_constraint=0.7,
            obj_n_mc=1)
    
    plot_elbo(approx, output_dir, n_iter)
    return sample_draw(approx, n_samp)

def runBayesInf_enzymes(BMCA_obj, r, data, output_dir, n_iter, n_samp=1):
    enzymes = ['e_' + i for i in r.getReactionIds()]
        
    known_e_inds = []
    omitted_e_inds = []
    for i, e in enumerate(enzymes):
        if e in data.columns:
            known_e_inds.append(i)
        else: 
            omitted_e_inds.append(i)
    e_inds = np.hstack([known_e_inds, omitted_e_inds]).argsort()

    with pm.Model() as pymc_model:

        # Initialize elasticities
        Ex_t = pm.Deterministic('Ex', initialize_elasticity(BMCA_obj.Ex.to_numpy(), name='Ex'))
        Ey_t = pm.Deterministic('Ey', initialize_elasticity(BMCA_obj.Ey.to_numpy(), name='Ey'))
        
        #Protein Expression Priors
        e_measured = pm.Normal('e_measured', mu=1, sigma=0.1, observed=BMCA_obj.en.T)
        e_unmeasured = pm.Normal('e_unmeasured', mu=1, sigma=0.1, shape=(len(omitted_e_inds), len(BMCA_obj.en)))
        e_t = at.concatenate([e_measured, e_unmeasured], axis=0)[e_inds, :]
        pm.Deterministic('e_t', e_t)
        
        chi_t = pm.Normal('chi_t', mu=0, sigma=0.5, observed=BMCA_obj.xn.T)
        y_t = pm.Normal('y_t', mu=0, sigma=0.5, observed=BMCA_obj.yn.T)
        
        likelihood = pm.Deterministic('vn', e_t * (np.ones((len(e_inds), len(BMCA_obj.en))) + pm.math.dot(Ex_t,chi_t) + pm.math.dot(Ey_t,y_t)))
        v_hat_obs = pm.Normal('v_hat_obs', mu=likelihood, sigma=0.1, observed=BMCA_obj.vn.squeeze().T)

        advi = pm.ADVI()
        tracker = pm.callbacks.Tracker(
            mean = advi.approx.mean.eval,
            std = advi.approx.std.eval
        )
    
        approx = advi.fit(
            n=n_iter, 
            callbacks = [tracker],
            obj_optimizer=pm.adagrad_window(learning_rate=1E-1), 
            total_grad_norm_constraint=0.7,
            obj_n_mc=1)

    plot_elbo(approx, output_dir, n_iter)
    return sample_draw(approx, n_samp)

def runBayesInf_fluxes(BMCA_obj, r, data, output_dir, n_iter, n_samp=1):
    flux = ['v_' + i for i in r.getReactionIds()]
        
    known_v_inds = []
    omitted_v_inds = []
    for i, v in enumerate(flux):
        if v in data.columns:
            known_v_inds.append(i)
        else: 
            omitted_v_inds.append(i)
    v_inds = np.hstack([known_v_inds, omitted_v_inds]).argsort()

    with pm.Model() as pymc_model:

        # Initialize elasticities
        Ex_t = pm.Deterministic('Ex', initialize_elasticity(BMCA_obj.Ex.to_numpy(), name='Ex'))
        Ey_t = pm.Deterministic('Ey', initialize_elasticity(BMCA_obj.Ey.to_numpy(), name='Ey'))
        
        # flux priors
        v_measured = pm.Normal('v_measured', mu=0, sigma=0.1, observed=BMCA_obj.vn.T)
        v_unmeasured = pm.Normal('v_unmeasured', mu=0, sigma=1, shape=(len(omitted_v_inds), len(BMCA_obj.vn)))

        v_t = at.concatenate([v_measured, v_unmeasured], axis=0)[v_inds, :]
        pm.Deterministic('v_t', v_t)

        chi_t = pm.Normal('chi_t', mu=0, sigma=0.5, observed=BMCA_obj.xn.T)
        y_t = pm.Normal('y_t', mu=0, sigma=0.5, observed=BMCA_obj.yn.T)

        #### NEED TO ADD fitting equation here
        e_ss = calculate_e_hat(BMCA_obj, v_t, Ex_t@chi_t, Ey_t@y_t)
        e_t = pm.Normal('e_t', mu=e_ss, sigma=1, observed=BMCA_obj.en.squeeze().T)

        advi = pm.ADVI()
        tracker = pm.callbacks.Tracker(
            mean = advi.approx.mean.eval,
            std = advi.approx.std.eval
        )
        approx = advi.fit(
            n=n_iter, 
            callbacks = [tracker],
            obj_optimizer=pm.adagrad_window(learning_rate=5E-3), 
            total_grad_norm_constraint=0.7,
            obj_n_mc=1)
        
    plot_elbo(approx, output_dir, n_iter)
    return sample_draw(approx, n_samp)

def runBayesInf_internal(BMCA_obj, r, data, output_dir, n_iter, n_samp=1):
    known_chi_inds = []
    omitted_chi_inds = []
    for i, sp in enumerate(r.getFloatingSpeciesIds()):
        if sp in data.columns:
            known_chi_inds.append(i)
        else: 
            omitted_chi_inds.append(i)
    chi_inds = np.hstack([known_chi_inds, omitted_chi_inds]).argsort()
    
    with pm.Model() as pymc_model:
    
        # Initialize elasticities
        Ex_t = pm.Deterministic('Ex', initialize_elasticity(BMCA_obj.Ex.to_numpy(), name='Ex'))
        Ey_t = pm.Deterministic('Ey', initialize_elasticity(BMCA_obj.Ey.to_numpy(), name='Ey'))
        
        chi_measured = pm.Normal('chi_measured', mu=0, sigma=0.1, observed=BMCA_obj.xn.T)
        chi_unmeasured = pm.Normal('chi_unmeasured', mu=0, sigma=10, shape=(len(omitted_chi_inds), len(BMCA_obj.xn)))

        chi_t = at.concatenate([chi_measured, chi_unmeasured], axis=0)[chi_inds, :]
        # supposedly chi_t would be in the order listed in ss tellurium

        pm.Deterministic('chi_t', chi_t)

        e_t = pm.Normal('e_t', mu=1, sigma=1, observed=BMCA_obj.en.T) # e_hat?
        y_t = pm.Normal('y_t', mu=0, sigma=10, observed=BMCA_obj.yn.T) # yn?

        likelihood = pm.Deterministic('vn', e_t * (np.ones(BMCA_obj.en.T.shape) + pm.math.dot(Ex_t,chi_t) + pm.math.dot(Ey_t,y_t)))
        v_hat_obs = pm.Normal('v_hat_obs', mu=likelihood, sigma=0.1, observed=BMCA_obj.vn.squeeze().T)

        advi = pm.ADVI()
        tracker = pm.callbacks.Tracker(
            mean = advi.approx.mean.eval,
            std = advi.approx.std.eval
        )
        approx = advi.fit(
            n=n_iter, 
            callbacks = [tracker],
            obj_optimizer=pm.adagrad_window(learning_rate=5E-3), 
            total_grad_norm_constraint=0.7,
            obj_n_mc=1)
        
    plot_elbo(approx, output_dir, n_iter)
    return sample_draw(approx, n_samp)

def runBayesInf_external(BMCA_obj, r, data, output_dir, n_iter, n_samp=1):
    external = r.getBoundarySpeciesIds()
    
    known_y_inds = []
    omitted_y_inds = []
    for i, y in enumerate(external):
        if y in data.columns:
            known_y_inds.append(i)
        else: 
            omitted_y_inds.append(i)
    y_inds = np.hstack([known_y_inds, omitted_y_inds]).argsort()

    with pm.Model() as pymc_model:

         # Initialize elasticities
        Ex_t = pm.Deterministic('Ex', initialize_elasticity(BMCA_obj.Ex.to_numpy(), name='Ex'))
        Ey_t = pm.Deterministic('Ey', initialize_elasticity(BMCA_obj.Ey.to_numpy(), name='Ey'))
        
        # flux priors
        y_measured = pm.Normal('y_measured', mu=0, sigma=0.1, observed=BMCA_obj.vn.T)
        y_unmeasured = pm.Normal('y_unmeasured', mu=0, sigma=0.1, shape=(len(omitted_y_inds), len(BMCA_obj.vn)))

        y_t = at.concatenate([y_measured, y_unmeasured], axis=0)[y_inds, :]
        pm.Deterministic('y_t', y_t)
        
        chi_t = pm.Normal('chi_t', mu=0, sigma=0.1, observed=BMCA_obj.xn.T)
        e_t = pm.Normal('e_t', mu=1, sigma=0.1, observed=BMCA_obj.en.squeeze().T)

        likelihood = pm.Deterministic('vn', e_t * (np.ones(BMCA_obj.en.T.shape) + pm.math.dot(Ex_t,chi_t) + pm.math.dot(Ey_t,y_t)))
        v_hat_obs = pm.Normal('v_hat_obs', mu=likelihood, sigma=0.1, observed=BMCA_obj.vn.squeeze().T)

        advi = pm.ADVI()
        tracker = pm.callbacks.Tracker(
            mean = advi.approx.mean.eval,
            std = advi.approx.std.eval
        )
        approx = advi.fit(
            n=n_iter, 
            callbacks = [tracker],
            obj_optimizer=pm.adagrad_window(learning_rate=1E-1), 
            total_grad_norm_constraint=0.7,
            obj_n_mc=1)
        
    plot_elbo(approx, output_dir, n_iter)
    return sample_draw(approx, n_samp)


def estimate_CCs(BMCA_obj, Ex, n_samp, a):
    a = np.diag(a)
    a = a[np.newaxis,:].repeat(n_samp * 1000, axis=0)

    Ex_ss = a @ Ex
    As = BMCA_obj.N @ np.diag(BMCA_obj.v_star) @ Ex_ss
    bs = BMCA_obj.N @ np.diag(BMCA_obj.v_star)
    bs = bs[np.newaxis, :].repeat(n_samp * 1000, axis=0)
    
    As = at.as_tensor_variable(As)
    bs = at.as_tensor_variable(bs)

    def solve_aesara(A, b):
        rsolve_op = LeastSquaresSolve()
        return rsolve_op(A, b).squeeze()

    CCC, _ = aesara.scan(lambda A, b: solve_aesara(A, b),
                        sequences=[As, bs], strict=True)

    identity = np.eye(len(BMCA_obj.N.T))
    identity = identity[np.newaxis,:].repeat(n_samp * 1000, axis=0)
    
    FCC = (Ex_ss @ CCC.eval()) + identity
    
    # return CCC.eval(), FCC
    return FCC


def append_FCC_df(postFCC, label, r):
    dfs=[]
    
    for idx, rxn in enumerate(r.getReactionIds()):
        # negativity applied here
        df = -pd.DataFrame(postFCC[:,idx,:], columns=r.getReactionIds())
        df['pt_rxn']=[rxn]*len(df)
        dfs.append(df)
    
    w = pd.concat(dfs)
    w['pt_str']=[label]*len(w)
    return w


def assign_teams(data_path, noise_bound, n_noisy, fluxes):
    all_data = pd.read_csv(data_path)
    
    all_players = set(all_data[fluxes].columns)
    team1 = set(random.sample(all_data[fluxes].columns.tolist(), n_noisy))
    team2 = all_players - team1

    # first, decide on how much noise you want to add
    noise_bound = noise_bound/100
    noise = 1 + np.random.uniform(-noise_bound, noise_bound, [n_noisy,1])

    all_data.loc[:, all_data.columns.isin(list(team1))] *= noise.squeeze()

    return team1, team2, all_data


def calculate_noisy_fba(team1, team2, all_data, cobra_file):
    fba_model = cobra.io.read_sbml_model(cobra_file)
    sample_flux_responses = []
    for sample in range(len(all_data)):
        for t in team1: # for each column
            rxn = t.split('_')[1]
            # set the fba bounds according to the the noise generator
            try: 
                fba_model.reactions.get_by_id(rxn).upper_bound = all_data[t][sample]
                fba_model.reactions.get_by_id(rxn).lower_bound = all_data[t][sample]
            except: 
                fba_model.reactions.get_by_id(rxn).lower_bound = all_data[t][sample]
                fba_model.reactions.get_by_id(rxn).upper_bound = all_data[t][sample]

        # find the max of team 1 fluxes and make that the upper and lower 
        # bounds for team 2 reactions. 
        dbl_bd = all_data.loc[:, all_data.columns.isin(list(team1))].max(axis=1)[sample] * 2
        for t in team2: # for each column
            rxn = t.split('_')[1]
            try: 
                fba_model.reactions.get_by_id(rxn).upper_bound = dbl_bd
                fba_model.reactions.get_by_id(rxn).lower_bound = -dbl_bd
            except: 
                fba_model.reactions.get_by_id(rxn).lower_bound = dbl_bd
                fba_model.reactions.get_by_id(rxn).upper_bound = -dbl_bd

        fba_model.reactions.vGLT.upper_bound = 5
        fba_model.reactions.vSUC.lower_bound = 0.01#0.05
        fba_model.reactions.vGLYCO.lower_bound = 0.01#0.05
        fba_model.reactions.vTreha.lower_bound = 0.01#0.05

        ## run FBA
        fba_model.objective = fba_model.reactions.vADH
        # add to sample_flux_responses
        sample_flux_responses.append(fba_model.optimize().fluxes)
        # reset the fba model
        fba_model = cobra.io.read_sbml_model(cobra_file)
    # return sample_flux_responses
    b = pd.concat(sample_flux_responses, axis=1).T.reset_index().drop(labels='index', axis=1)
    b.columns = ['v_'+i for i in b.columns]
    return b

def make_FBA_test_data(ant_file, sbml_path, data_path, n_noisy, noise_bound, pt):
    """
    ant_file = '../data/interim/Antimony/Simplified_Teusink_yeast.ant'
    sbml_path = "../data/interim/sbml/Simplified_Teusink_yeast_cobra.xml" ## load the fba model
    data_path="../data/interim/generated_data/simplTeusink-noReg/Simplified_Teusink_yeast_3.csv"

    n_noisy = 5 #(int) number of species that will have random noise added
    noise_bound=2 #(int) maximum amount of noise that will be applied to species (in %)
    """
    
    fba_model = cobra.io.read_sbml_model(sbml_path)
    N = cobra.util.array.create_stoichiometric_matrix(fba_model)

    r = te.loada(ant_file)
    enzymes = ['e_' + i for i in r.getReactionIds()]
    internal = r.getFloatingSpeciesIds()
    external = r.getBoundarySpeciesIds()
    fluxes = ['v_' + i for i in r.getReactionIds()]
    
    i = 0
    while i < 100:
        diagnosis=[]
        fba_model = cobra.io.read_sbml_model(sbml_path)
        team1, team2, all_data = assign_teams(data_path, noise_bound, n_noisy, fluxes)
        a = calculate_noisy_fba(team1, team2, all_data, sbml_path)
        for i in range(len(a)):
            if np.all(np.isclose((N@a.iloc[i]),0)):
                diagnosis.append(True)
            else:
                diagnosis.append(False)
        if all(diagnosis) == True:
            print('complete')
            c = pd.concat([all_data[enzymes+internal+external],a], axis=1)
            c.to_csv(f'SimplTeusink_fba-noise_pt{pt}_sp{n_noisy}-max{noise_bound}.csv', index=False)
            break
        else:
            i+=1

    if i==100:
        print('could not complete in 100 iterations')


def bootstrap_spearman(x, y, num_bootstrap=1000, alpha=0.05):
    n = len(x)
    corr_list = []

    # Original Spearman correlation
    corr_original, p_value = spearmanr(x, y)

    for _ in range(num_bootstrap):
        # Generate bootstrap samples
        indices = np.random.randint(0, n, n)
        x_bootstrap = [x[i] for i in indices]
        y_bootstrap = [y[i] for i in indices]

        # Calculate Spearman correlation for the bootstrap sample
        corr, _ = spearmanr(x_bootstrap, y_bootstrap)
        corr_list.append(corr)

    # Convert to numpy array for convenience
    corr_list = np.array(corr_list)
    
    # Calculate the confidence intervals
    lower_bound = np.percentile(corr_list, (alpha/2) * 100)
    upper_bound = np.percentile(corr_list, (1 - alpha/2) * 100)
    
    return corr_original, p_value, lower_bound, upper_bound


def get_az_summary(t): 
    Ex_mean = az.summary(t)['mean'].reset_index()
    Ex_mean.columns = ['elasticity', 'mean']
    Ex_mean = Ex_mean[Ex_mean.elasticity.str.contains("Ex\[")]['mean'].values.flatten().reshape((-1,1))
    return np.mean(Ex_mean, axis=1)


def get_az_mean(traces): 
    trace_means = []
    for t in traces: 
        Ex_mean = az.summary(t)['mean'].reset_index()
        Ex_mean.columns = ['elasticity', 'mean']
        Ex_mean = Ex_mean[Ex_mean.elasticity.str.contains("Ex\[")]['mean'].values.flatten().reshape((-1,1))
        trace_means.append(Ex_mean)
    Ex = np.concatenate(trace_means, axis=1)
    return np.mean(Ex, axis=1)

def calculate_slope(x,y): 
    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x,y)
    return slope, intercept, r_value


def plt_spr_scatter(dataframe, title):
    plt = sns.scatterplot(data=dataframe, x='perturbation', y="r", hue='omit', style='omit', s=100, alpha=0.8, zorder=100)
    plt.grid(True, which='both', axis='both', zorder=0)
    plt.set(xscale='log')
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.set_xlabel('fold change in enzyme concentration', size=14)
    plt.set_ylabel('Spearman rank coefficient, $\it{r}$', size=14)
    plt.tick_params(axis='both', which='major', labelsize=13)
    plt.set_title(title, size=20)
    