# handy-dandy
import os
import sys
from tqdm import tqdm
import winsound
duration = 1000  # milliseconds
freq = 440  # Hz

# arrays/dataframes
import numpy as np
np.random.seed(0)
np.set_printoptions(threshold=sys.maxsize)

import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

from csv import writer

# math/stats
import scipy
import scipy.stats

# biochemical pathway simulators
import tellurium as te
import libsbml
import cobra

# plotting
import matplotlib.pyplot as plt

class BMCA():

    def __init__(self, model_file, data_file, desired_product=None, bd_exclude=[]):
                
        self.N, self.v_star, self.en, self.xn, self.yn, self.vn = \
            BMCA.load_model_data(model_file, data_file, desired_product, bd_exclude)
        """
        s = te.loada("models/MODEL1303260011_cobra.ant")
        with open("temp.txt", "w") as f:
            f.write(s.getSBML())
        model = cobra.io.read_sbml_model("temp.txt")
        os.remove("temp.txt") """

        self.Ex = BMCA.create_Visser_elasticity_matrix(model_file)
        self.Ey = BMCA.create_Visser_elasticity_matrix(model_file, Ex=False)


    def load_model_data(model_file, data_file, desired_product, bd_exclude):
        # a function, because a method takes in 'self'
        """
        this method takes in an SBML model and csv data to establish important
        attributes
        """
        r = te.loada(model_file)
        data = pd.read_csv(data_file)

        # sorting the data
        enzymes = ['e_' + i for i in r.getReactionIds()]
        e = data[enzymes]
        x = data[r.getFloatingSpeciesIds()]
        # y = data[r.getBoundarySpeciesIds()]
        y = data[[i for i in r.getBoundarySpeciesIds() if i not in bd_exclude]]
        v = data[['v_' + i for i in r.getReactionIds()]]

        # normalizing the data
        ref_ind = data.idxmax()[desired_product]   
        e_star = e.iloc[ref_ind].values
        x_star = x.iloc[ref_ind].values
        y_star = y.iloc[ref_ind].values
        v_star = v.iloc[ref_ind].values

        e_star[e_star == 0] = 1e-6
        x_star[x_star == 0] = 1e-6
        y_star[y_star == 0] = 1e-6
        v_star[v_star == 0] = 1e-9

        # Normalize to reference values (and drop trivial measurement)
        en = e.divide(e_star)
        xn = x.divide(x_star)
        yn = y.divide(y_star)
        vn = v.divide(v_star)

        en.drop(en.index[ref_ind], inplace=True)
        xn.drop(xn.index[ref_ind], inplace=True)
        yn.drop(yn.index[ref_ind], inplace=True)
        vn.drop(vn.index[ref_ind], inplace=True)

        N = r.getFullStoichiometryMatrix()
        
        # Correct negative flux values at the reference state
        N[:, v_star < 0] = -1 * N[:, v_star < 0]
        v_star = np.abs(v_star)
    
        return N, v_star, en, xn, yn, vn
    
    def create_Visser_elasticity_matrix(model_file, Ex=True):
        """Create an elasticity matrix for internal metabolites given the model in model.

        E[j,i] represents the elasticity of reaction j for metabolite i.

        """
        r = te.loada(model_file)
        # convert from antimony to sbml
        doc = libsbml.readSBMLFromString(r.getSBML())
        model = doc.getModel()

        if Ex: 
            m_list = r.getFloatingSpeciesIds()
            n_metabolites = len(r.getFloatingSpeciesIds())
        else:
            m_list = r.getBoundarySpeciesIds()
            n_metabolites = len(r.getBoundarySpeciesIds())
        n_reactions = len(r.getReactionIds())
        array = np.zeros((n_reactions, n_metabolites), dtype=float)

        for n in range(n_reactions): 
            rxn = model.getReaction(n)
            for reactant in range(rxn.getNumReactants()):                 
                metabolite = rxn.getReactant(reactant).species
                stoich = rxn.getReactant(reactant).getStoichiometry()
                if metabolite in m_list: 
                    # Reversible reaction, assign all elements to -stoich
                    if rxn.getReversible(): 
                        array[n, m_list.index(metabolite)] = -np.sign(stoich)

                    # Irreversible in forward direction, only assign if met is reactant
                    elif ((not rxn.getReversible()) & (stoich < 0)):
                        array[n, m_list.index(metabolite)] = -np.sign(stoich)
            
                    # Irreversible in forward direction, only assign if met is reactant
                    elif ((not rxn.getReversible()) & (stoich > 0)): 
                        array[n, m_list.index(metabolite)] = -np.sign(stoich)
        
        return array
