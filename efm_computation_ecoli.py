# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:28:29 2023

@author: guill
"""

import numpy as np
import cobra
from efmtool import calculate_efms
import time
from scipy.optimize import linprog
from mip import *
import tqdm, sys
from multiprocessing import Pool
#from cobra.util.array import create_stoichiometric_matrix
#from cobamp.wrappers import KShortestEFMEnumeratorWrapper
#from mip import CPLEX

###########################################################

digit_tol = 6
def supp(vector, tol=digit_tol):
    return (list(np.nonzero(np.round(vector, tol))[0]))

def irr_supp(vector, zero_tol=digit_tol):
    #return list(np.intersect1d(supp(vector),supp(self.irr,zero_tol)))
    return (np.intersect1d(supp(vector), supp(irr,zero_tol)))

def candidates(vector,efvs):
     return efvs[np.where(np.all((np.round(efvs[:,np.setdiff1d(irr_v, supp(vector))],5) == 0), axis=1))]    #return efvs[np.where(np.all((np.round(efvs[:, np.setdiff1d(supp(m,), irr_supp(vector))], 5) == 0), axis=1))]
     
def decompose(vector, candidates):
    c=np.zeros(len(candidates))
    A=candidates.T
    b=vector
    decomp=linprog(c, A_eq=A, b_eq=b)
    return (decomp)

def decompose_efm(efm):
    CAN = candidates(efm,bases)
    return MILP_Shortest_decomp(efm,CAN)

def MILP_Shortest_decomp(target_vector, CAN, tol=1e-6):
    m = mip.Model(solver_name="CPLEX")
  
    #Numeric tolerance:
    #m.infeas_tol = 0.0000001
    #m.integer_tol = 0.0000001
    m.verbose = 0
  
    a = [
        m.add_var(var_type=mip.BINARY)
        for i in range(len(CAN))
    ]
  
    x = [
        m.add_var(var_type=mip.CONTINUOUS)
        for i in range(len(CAN))
    ]
  
    M = 1000000 #Big M
    #Need to make sure if a[i]=0 --> x[i]=0
    for i in range(len(CAN)):
      #m.add_constr(a[i] <= x[i])
      m.add_constr(x[i] <= M*a[i])
      m.add_constr(0 <= x[i])
  
    for flux in range(len(target_vector)):
      m += mip.xsum(
          x[i]*CAN[i][flux]
          for i in range(len(CAN))) >= target_vector[flux]-tol
	
	
    for flux in range(len(target_vector)):
      m += mip.xsum(
          x[i]*CAN[i][flux]
          for i in range(len(CAN))) <= target_vector[flux]+tol
h  
  
    m += mip.xsum(
        a[i]
        for i in range(len(CAN))) >= 1
  
    m.objective = mip.minimize(mip.xsum(
        a[i]
        for i in range(len(CAN))
    ))
  
    # Start the solver.
    solver_status = m.optimize()
  
    decomp_supp = [a[i].x for i in range(len(CAN))]
    decomp_coeff = [x[i].x for i in range(len(CAN))]
      
    contributors=[]
    coeff=[]
    if str(solver_status) == 'OptimizationStatus.OPTIMAL':
        # List with the coefficients of the EFMs involved in the decomp
        #coeff = [i for i in decomp_coeff if abs(i) > 0]   #To avoid numerical issues with 0
        for i in range(len(decomp_coeff)):
            if decomp_supp[i] != 0:  # Append only those EFMs that have a contribution coeff !=0
              contributors.append(CAN[i])  # Append the candidate EFM i
              coeff.append(decomp_coeff[i])
              
    objective_value=m.objective_value
    m.clear()  # Clears the model variables to avoid possible issues when doing multiple MILPs consecutively
    return (coeff)


def efms_carbon(model_path, sources):
    model=cobra.io.sbml.read_sbml_model(model_path)
    r_id=[r.id for r in model.reactions]
    
    #Automatise environment creation
    model.reactions.get_by_id("EX_o2_e").bounds=(0, 1000)
    model.reactions.get_by_id("EX_glc__D_e").bounds=(0, 1000)
    for source in sources:
        model.reactions.get_by_id(r_id[source]).bounds=(-10, 1000)    
    
    #model_bounds={r.id: (r.lower_bound, r.upper_bound)
    #                     for r in model.reactions}

    # Get stoichiometry, reversibility, reaction names and metabolite names
    stoichiometry=cobra.util.array.create_stoichiometric_matrix(model)
    rev=np.array([rea.reversibility for rea in model.reactions]).astype(int)
    irr_v=supp((np.ones(len(rev)) - rev).astype(int))
    met_id=[met.id for met in model.metabolites]

    # Compute efms using efmtool
    t_start=time.time()
    efms_efmtool=calculate_efms(stoichiometry=stoichiometry,
                                  reversibilities=rev,
                                  reaction_names=r_id,
                                  metabolite_names=met_id)
    time_efmtool=time.time() - t_start
    print(f"EFMTool found {efms_efmtool.shape[1]} EFMs in {round(time_efmtool * 1000, 0)} ms.")
    
    return(efms_efmtool.T, irr_v)
###########################################################


"""
Load model and create nutrient environments
"""
if __name__ == '__main__':

    model_path=r"C:\Users\guill\Documents\UNIVERSIDAD\Bioinformatics FU\Thesis\scripts\e_coli_core.xml"

    #Automatise environment creation
    '''To change the pair of C-sources, change the numbers to appropiate reaction'''
    source_A = 19
    source_B = 27
    efms_source_A = efms_carbon(model_path, [source_A])[0]
    efms_A_supp = [supp(efm) for efm in efms_source_A]
    efms_source_B = efms_carbon(model_path, [source_B])[0]
    efms_B_supp = [supp(efm) for efm in efms_source_B]
    efms_both,irr_v = efms_carbon(model_path, [source_A,source_B])

    to_decomp=[]
    bases = np.r_[efms_source_A,efms_source_B]
    for efm in efms_both:
        if supp(efm) not in efms_A_supp and supp(efm) not in efms_B_supp:
            to_decomp.append(efm)
        
    print('There are', len(to_decomp),'to decompose.', len(bases), 'EFMs as base.')

    #This needs to be converted to NumPy arrays for the candidate function to work.
    np_decomp = np.array(to_decomp) #NumPy array just in case
    
    np.save('all_candidates', bases)
    np.save('efms', np_decomp)
    np.save('irr_v', irr_v)
    sys.exit()