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