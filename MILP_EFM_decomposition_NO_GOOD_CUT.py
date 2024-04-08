from mip import *
import numpy as np
from scipy.optimize import linprog

def supp(vector): #RETURNS INDECES!!!
    return np.nonzero(vector)
    #return np.where(vector != 0, 0, 1)

def irr_supp(vector): #RETURNS A BINARY VECTOR WHERE 1 IS IRR AND 0 REV
    return np.where(vector == 0, 1, 0)

def candidates(vector,efvs):
    return efvs[np.where(np.all((np.round(efvs[:,np.setdiff1d(irr_v, supp(vector))],5) == 0), axis=1))]

def decompose(vector,candidates):
    c = np.zeros(len(candidates))
    A = candidates.T
    b = vector
    decomp = linprog(c,A_eq = A ,b_eq = b)
    return(decomp)

def clean_supp(vector):
  support = str(supp(vector))
  support = support[7:]
  support = support[:-3]
  support = support.replace('\t','')
  support = support.replace('\n','')
  support = support.replace(' ','')
  return support

#Target vector is the vector you want to decompose. CAN is an array of all possible vectors that could be contributors
def MILP_Shortest_decomp(target_vector, CAN, past_sol):
    m = mip.Model(solver_name="CPLEX")
  
    #Numeric tolerance:
    m.infeas_tol = 0.0000001
    m.integer_tol = 0.0000001
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
          for i in range(len(CAN))) == target_vector[flux]

    #NO-GOOD CUT
    for psol in past_sol:
      m += mip.xsum(
         a[i]
         for i in past_sol[psol]) <= len(past_sol[psol])-1
    
    m += mip.xsum(
       a[i]
       for i in range(len(CAN))) >= 1
    
    m.objective = mip.minimize(mip.xsum(
       a[i]
       for i in range(len(CAN))))
    
    # Start the solver.
    solver_status = m.optimize()
    
    decomp_supp = [a[i].x for i in range(len(CAN))]
    decomp_coeff = [x[i].x for i in range(len(CAN))]
           
    contributors=[]
    coeff=[]
    if str(solver_status) == 'OptimizationStatus.OPTIMAL':
        ID = str(len(past_sol.keys()))
        past_sol[ID] = []
        # List with the coefficients of the EFMs involved in the decomp
        #coeff = [i for i in decomp_coeff if abs(i) > 0]   #To avoid numerical issues with 0
        for i in range(len(decomp_coeff)):
            if decomp_supp[i] != 0:  # Append only those EFMs that have a contribution coeff !=0
              contributors.append(CAN[i])  # Append the candidate EFM i
              coeff.append(decomp_coeff[i])
              past_sol[ID].append(i)
              
    objective_value=m.objective_value
    m.clear()  # Clears the model variables to avoid possible issues when doing multiple MILPs consecutively
    return (solver_status, contributors, coeff, objective_value)


def MILP_NO_GOOD_CUT(target_vector, CAN): #TV is vector to decompose, Candidates for decomposition (np arrays)
  past_sol = {}
  sol_count = 0
  M = MILP_Shortest_decomp(target_vector, CAN, past_sol)
  print(M)
  sol_count += 1
  i = 0
  if len(CAN) == M[3]: #To avoid computing infeasible problems
    print(sol_count, 'solutions until bigger support.')
    i=1

  while i==0:
    M2 = MILP_Shortest_decomp(target_vector, CAN, past_sol)
    print(M2)
    if len(M[1])<len(M2[1]):
      print(sol_count, 'solutions until bigger support.')
      i=1
    sol_count += 1

  return sol_count


if __name__ == "__main__":
    past_sol={}
    CAN = np.array([[ 1.,  0., -3.,  2.,  0.],
                    [ 1., -3.,  0.,  2.,  0.],
                    [-3.,  1.,  0.,  0.,  2.],
                    [-3.,  0.,  1.,  0.,  2.],
                    [ 0.,  0., -4.,  3.,  1.],
                    [-4.,  0.,  0.,  1.,  3.],
                    [ 0.,  1., -1.,  0.,  0.],
                    [ 0., -1.,  1.,  0.,  0.]])

    tv = np.array([ 0., -4.,  0.,  3.,  1.])
    
    print(MILP_NO_GOOD_CUT(tv, CAN))