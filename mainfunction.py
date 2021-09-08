from math import gcd, ceil
from scipy import sparse
import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
import Engproblem
import Optproblem
import Mathproblem

if __name__ =='__main__':
##Solve Engineering problem
    width, height, lc, lt = 10.0 ,10.0, 1, 1 
    vt = [[0, 0], [width, 0], [width, height], [0, height]]
    material = Engproblem.Material(lc, lt)
    domain = Engproblem.Domain(vt)
    if domain.verify():
        loads = Engproblem.Load([width/2, height], [1, 0], 1)
        supports = [Engproblem.Support([0, 0], [0, 0]), Engproblem.Support([width, 0], [1, 0])]
    
##Solve Optimization problem
        fl, step = 1.42, 0.25
        groundstructure = Optproblem.Groundstructure(domain, loads, supports, material, step, step)
        prob = Optproblem.Optprob(groundstructure)
        SPV = prob.createSPV()
        LDV = prob.createLDV()
        MatrixB = prob.createMatrixB(SPV)
        LV = prob.createLV()
        ML = groundstructure.createmembers()
    
##Solve Mathematic problem
        mathprob = Mathproblem.Mathproblem(LV, MatrixB, LDV, SPV, lc, lt)
        vol, q, a = mathprob.solve()
        print("final vol: %f" % (vol))
        prob.plot(q, a, max(a) * 1e-3, "Finished")
    else:print("error")
