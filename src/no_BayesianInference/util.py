# from PSJ's emll/util.py

import numpy as np
import scipy as sp
import tellurium as te
import os
import re
import csv
import numpy as np


def calculate_steady_state(Ea, Eb, N, v_star, en, yn):
        # equation 5 of PSJ's paper
        v_e = np.diag(np.matmul(v_star, en))
        N_v_e = np.matmul(N, v_e)
        A = np.matmul(N_v_e, Ea)
        
        inner_v = (np.ones(N.shape[1]) + np.matmul(Eb, yn))
        B = -np.matmul(N_v_e, inner_v)

        A_pinv = sp.linalg.pinv(A)
        xn = np.matmul(A_pinv, B)
        
        vn = en * (np.ones(N.shape[1]) +
                   np.matmul(Ea, xn) +
                   np.matmul(Eb, yn))

        return en, vn


def generate_data(model_file, perturbation_levels, data_folder):
    """
    Takes in SBML model_file and creates perturbation data files. 
    Deposits perturbation data files into data_folder
    """
    model_name = model_file.split('.')[0]
    model_name = model_name.split('/')[-1]
    for pl in perturbation_levels:
        
        r = te.loada(model_file)
        
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
                
                # perturbed boundary species cases
                for params in exMet:
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
    friendly format by removing all boundary species. The returned 
    file is in SBML
    
    r = te.loada(antimony_file)
    bd_sp = r.getBoundarySpecies()
    # read the file 
    outfile = open(antimony_file,"r")
    data = outfile.readlines()

    # set up 

    linenum=[]
    title=[]
    for n, line in enumerate(data):
        if '// ' in line:
            linenum.append(n)
            title.append(line)

    compmt_idx = [title.index(a) for a in title if 'Compartments' in a][0]

    rxn_idx = [title.index(a) for a in title if 'Reactions' in a][0]
    """
    pass
    
    # identify all the boundary species
    # between lines linenum[compmt_idx] and linenum[compmt_idx + 1] 
    # for each line
    # split the string by comma
    # if an element in the returned list contains '$', then remove that element
    # write the new list as a string into the file. 

    # DELETE ALL THE RATE LAWS
    # between lines linenum[rxn_idx + 1]  and linenum[rxn_idx + 1] 
    # split the line by ';' 
    # write the first element and add '; ;' after it

    # DELETE ALL BOUNDARY SPECIES IN THE REACTIONS
    # between lines linenum[rxn_idx + 1]  and linenum[rxn_idx + 1] 
    # split by ' '
    # if element contains '$', delete
    # reevaluate list; if '+' is next to a punctuation mark, then delete
    # write the resulting line back into file

    # DELETE ALL BOUNDARY SPECIES INITIALIZATIONS
    # for lines linenum[rxn_idx + 1] until the end of the document
    # if line contains ay element in bd_sp, delete the line

    # return cobra_file