# handy-dandy
import os
import sys
from tqdm import tqdm
     
# arrays/dataframes
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

import pandas as pd

from csv import writer

# math/stats
import scipy as sp
import scipy.stats

# biochemical pathway simulators
import tellurium as te
import libsbml # https://synonym.caltech.edu/software/libsbml/5.18.0/docs/https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/annotated.html
import cobra

import emll

class antemll():

    def __init__(self, model_file, data_file, cobra_sbml, ref_ind=0, filler_v_star=None):
                
        self.model_file = model_file
        self.model = te.loada(model_file)

        self.N, self.v_star, self.en, self.xn, self.yn, self.vn, self.x_star = \
            antemll.load_model_data(model_file, data_file, ref_ind, filler_v_star)
                
        self.n_exp = self.en.shape[0]

        self.Ex = antemll.create_Visser_elasticity_matrix(model_file, cobra_sbml)
        self.Ey = antemll.create_Visser_elasticity_matrix(model_file, cobra_sbml, Ex=False)


    def load_model_data(model_file, data_file, ref_ind, filler_v_star):
        # a function, because a method takes in 'self'
        """
        this method takes in an SBML model and csv data to establish important
        attributes
        """
        r = te.loada(model_file)
        r.conservedMoietyAnalysis = True
        
        if isinstance(data_file, str):
            df = pd.read_csv(data_file)
        elif isinstance(data_file, pd.DataFrame):
            df = data_file

        # in case of omitted data
        available_fl_sp = [i for i in r.getFloatingSpeciesIds() if i in df.columns]
        available_bd_sp = [i for i in r.getBoundarySpeciesIds() if i in df.columns]
        
        enzymes = ['e_' + i for i in r.getReactionIds()]
        available_enz = [i for i in enzymes if i in df.columns]

        fluxes = ['v_' + i for i in r.getReactionIds()]
        available_fluxes = [i for i in fluxes if i in df.columns]
        
        # clean the data
        data = df.drop(df[df[available_fl_sp].lt(0).any(axis=1)].index)
        
        # sorting the data
        
        e = data[available_enz]
        x = data[available_fl_sp]
        y = data[available_bd_sp]
        v = data[available_fluxes]

        # normalizing the data
        #  'ref' should be the first row of data
        # ref_ind = data.idxmax()[desired_product]   
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
        # print(en)
        xn = x.divide(x_star)
        yn = y.divide(y_star)
        vn = v.divide(v_star)

        #en.drop(en.index[ref_ind], inplace=True)
        #xn.drop(xn.index[ref_ind], inplace=True)
        #yn.drop(yn.index[ref_ind], inplace=True)
        #vn.drop(vn.index[ref_ind], inplace=True)

        N = r.getFullStoichiometryMatrix()
        
        """filler_v_star = list(v_star).copy()
        filler_indices = set(range(len(fluxes))) - set([fluxes.index(i) for i in available_fluxes])
        for i in filler_indices:
            filler_v_star.insert(i, 1)
        filler_v_star = np.array(filler_v_star)"""

        if filler_v_star is not None:
            v_star = filler_v_star
        # Correct negative flux values at the reference state
        N[:, v_star < 0] = -1 * N[:, v_star < 0]
        v_star = np.abs(v_star)

        assert np.all(np.isclose(np.matmul(N, v_star), 0, atol=1e-07)), "data does not describe steady state"
        
        yn[yn == 0] = 1E-6

        xn = np.log(xn)
        yn = np.log(yn)

        return N, v_star, en, xn, yn, vn, x_star
    
    def create_Visser_elasticity_matrix(model_file, cobra_smbl, Ex=True):
        """Create an elasticity matrix for metabolites given the model in model.
        E[j,i] represents the elasticity of reaction j for metabolite i.
        stolen from emll.util.create_elasticity_matrix
        """        
        if Ex: # making the Ex matrix
            cobra_ant = model_file.split('.ant')[0]
            cobra_ant = cobra_ant.split('-')[0] + '_cobra.ant'
            print(cobra_ant)
                # check if cobra version exists 
                # if not, make it.
                # load the cobra model

            r = te.loada(cobra_ant)
            r.conservedMoietyAnalysis = True
            model = cobra.io.read_sbml_model(cobra_smbl)
            # model = cobra.io.read_sbml_model("../../models/sbml/Simplified_Teusink_yeast_cobra.xml")
            # model = cobra.io.read_sbml_model("../../../models/sbml/JSexample22_cobra.xml")

            n_metabolites = len(model.metabolites)
            n_reactions = len(model.reactions)
            array = np.zeros((n_reactions, n_metabolites), dtype=float)

            m_ind = model.metabolites.index
            r_ind = model.reactions.index

            for reaction in model.reactions:
                for metabolite, stoich in reaction.metabolites.items():

                    # Reversible reaction, assign all elements to -stoich
                    if reaction.reversibility:
                        array[r_ind(reaction), m_ind(metabolite)] = -np.sign(stoich)

                    # Irrevesible in forward direction, only assign if met is reactant
                    elif ((not reaction.reversibility) & 
                        (reaction.upper_bound > 0) &
                        (stoich < 0)):
                        array[r_ind(reaction), m_ind(metabolite)] = -np.sign(stoich)

                    # Irreversible in reverse direction, only assign if met is product
                    elif ((not reaction.reversibility) & 
                        (reaction.lower_bound < 0) &
                        (stoich > 0)):
                        array[r_ind(reaction), m_ind(metabolite)] = -np.sign(stoich)

            array =  emll.util.create_elasticity_matrix(model)
            array = pd.DataFrame(array, index=r.getReactionIds(), columns=[i.id for i in model.metabolites])
            array = array.loc[:,r.getFloatingSpeciesIds()]
            
        else: # making the Ey matrix
            r = te.loada(model_file)
            r.conservedMoietyAnalysis = True
            
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
            array = pd.DataFrame(array, index=r.getReactionIds(), \
                                 columns=r.getBoundarySpeciesIds())

        return array

    def calculate_Smallbone_ss(self, Ea, Eb):
        """
        As outlined in Smallbone et al. 2007, Equations 6 and 10

        Because compartment volumes usually have a value of 1, I have 
        simply opted to make np.ones matrices instead of calling on 
        roadrunner for the compartment volumes.

        """
        r = te.loada(self.model_file)
        r.conservedMoietyAnalysis = True
        a = r.getIndependentFloatingSpeciesIds()
        b = r.getFloatingSpeciesIds()
        squiggle_idx = [b.index(i) for i in a if i in b]
        squiggle_idx.sort()
        
        t1 = np.linalg.inv(np.diag(self.x_star))
        t2 = np.linalg.inv(np.diag(np.ones(len(self.x_star))))
        t3 = self.N
        t4 = np.linalg.pinv(self.N[squiggle_idx,:])
        t5 = np.diag(np.ones(len(squiggle_idx)))
        t6 = np.diag(self.x_star[squiggle_idx])
        L = t1 @ t2 @ t3 @ t4 @ t5 @ t6 # link matrix
        
        L_ = r.getLinkMatrix()

        # solving for equation 10, steady state internal metabolite concentrations
        t7 = self.N[squiggle_idx,:]
        t7_ = r.getReducedStoichiometryMatrix()
        t8 = np.diag(self.v_star)
        
        # d = t7_ @ t8 @ Ea @ L
        # print(np.linalg.matrix_rank(d)) 

        t100 = -np.linalg.inv(t7_ @ t8 @ Ea @ L_) 
        t101 = t7_ @ t8 @ Eb @ (self.yn).T # here, we may need to add in another dimension
        chi_star = t100 @ t101

        # chi_star.rename({1:'met_conc'}, axis=1, inplace=True)
        res = {squiggle_idx[i]: b[i] for i in range(len(squiggle_idx))}
        chi_star = chi_star.rename(index=res)

        return chi_star
    
    def calculate_PSJ_ss(self, Ea, Eb):
        # equation 5 of PSJ's paper
        v_e = np.diag(np.squeeze(self.en * self.v_star).to_numpy())
        N_v_e = self.N @ v_e # would use @ instead of * for only 1 perturbation
        A = np.matmul(N_v_e, Ea)
        
        print(self.yn)
        inner_v = (np.ones((self.N.shape[1], self.n_exp)) + np.matmul(Eb, (self.yn).T))
        B = -np.matmul(N_v_e, inner_v)

        A_pinv = sp.linalg.pinv(A)
        chi = np.matmul(A_pinv, B)
        
        return chi
