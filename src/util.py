# from PSJ's emll/util.py

import numpy as np
import scipy as sp
import tellurium as te
import os
import re
import csv
import numpy as np
import libsbml


def generate_data(model_file, perturbation_levels, data_folder):
    """
    Takes in SBML model_file and creates perturbation data files. 
    Deposits perturbation data files into data_folder
    """
    model_name = model_file.split('.')[0]
    model_name = model_name.split('/')[-1]
    for pl in perturbation_levels:
        
        r = te.loada(model_file)
        r.conservedMoietyAnalysis = True
        
        exMet = r.getBoundarySpeciesIds()
        inMet = r.getFloatingSpeciesIds()
        fluxIDs = ['v_' + i for i in r.getReactionIds()]
        e_list = [i for i in r.getGlobalParameterIds() if 'e_' in i]   
        
        pertLevel = pl/100 
        perturbation_level = [1 - pertLevel, 1 + pertLevel]
        
        header = e_list + exMet + inMet + fluxIDs        
        
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
                """
                # perturbed enzyme cases
                for params in e_list:
                    for level in perturbation_level:
                        r.resetToOrigin()
                        r.setValue(params, level*r.getValue(params))
                        
                        spConc = list(r.simulate(0,1000000)[-1])[1:]
                        # r.steadyState()
                        enzymes = [r.getValue(e) for e in e_list]
                        exMet_values = [r.getValue(m) for m in exMet]
                        fluxes = list(r.getReactionRates())
                        
                        writer.writerow(enzymes + exMet_values + spConc + fluxes)
                """
                # perturbed boundary species cases
                for params in ['GLCo']:
                    for level in perturbation_level:
                        r.resetToOrigin()
                        r.setValue(params, level*r.getValue(params))
                        
                        spConc = list(r.simulate(0,1000000)[-1])[1:]
                        # r.steadyState()
                        enzymes = [r.getValue(e) for e in e_list]
                        exMet_values = [r.getValue(m) for m in exMet]
                        fluxes = list(r.getReactionRates())
                        
                        writer.writerow(enzymes + exMet_values + spConc + fluxes)        
            except:
                print('failed at', params, level) # pass #



def ant_to_cobra(antimony_file):
    """
    This method takes in an Antimony file and converts it to a Cobra-
    friendly format by removing all boundary species and replacing all
    rate laws with a constant. The returned file is in SBML
    
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

import aesara.tensor as pt
import pymc as pm

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
