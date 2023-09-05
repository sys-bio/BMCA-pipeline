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
import scipy as sp
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
        
        self.n_exp = self.en.shape[0]

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
        y = data[r.getBoundarySpeciesIds()]
        # y = data[[i for i in r.getBoundarySpeciesIds() if i not in bd_exclude]]
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

        assert np.isclose(np.all(np.matmul(N, v_star)), 0), "data does not describe steady state"
        
        yn[yn == 0] = 1E-6

        return N, v_star, en, xn, yn, vn
    
    def create_Visser_elasticity_matrix(model_file, Ex=True):
        """Create an elasticity matrix for metabolites given the model in model.
        E[j,i] represents the elasticity of reaction j for metabolite i.
        """
        r = te.loada(model_file)
        
        if Ex:
            array = -r.getFullStoichiometryMatrix().T
        else: 
            doc = libsbml.readSBMLFromString(r.getSBML())
            model = doc.getModel()
            
            rxns = r.getReactionIds()
            bd_sp = r.getBoundarySpeciesIds()
            
            array = np.zeros((len(rxns), len(bd_sp)))
            for n in range(len(rxns)): 
                rxn = model.getReaction(n)
                for reactant in range(rxn.getNumReactants()):                 
                    sp = rxn.getReactant(reactant).species
                    stoich = rxn.getReactant(reactant).getStoichiometry()
                    if sp in bd_sp: 
                        array[n, bd_sp.index(sp)] = -np.sign(stoich)
                for prod in range(rxn.getNumProducts()):
                    sp = rxn.getProduct(prod).species
                    stoich = rxn.getProduct(prod).getStoichiometry()
                    if sp in bd_sp: 
                        array[n, bd_sp.index(sp)] = -np.sign(stoich)
        return array

    def calculate_steady_state(self, Ea, Eb):
        # equation 5 of PSJ's paper
        v_e = np.diag(np.matmul(self.en, self.v_star.reshape((-1,1))))
        N_v_e = self.N * v_e
        A = np.matmul(N_v_e, Ea)
        
        inner_v = (np.ones((self.N.shape[1], self.n_exp)) + np.matmul(Eb, np.log(self.yn).T))
        B = -np.matmul(N_v_e, inner_v)

        A_pinv = sp.linalg.pinv(A)
        chi = np.matmul(A_pinv, B)
        
        v_hat = np.diag(np.matmul((np.ones((self.N.shape[1], self.n_exp)) +
                np.matmul(Ea, chi) +
                np.matmul(Eb, self.yn.T)), self.en))
        
        return chi, v_hat
