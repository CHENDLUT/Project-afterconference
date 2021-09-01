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
    width, height, lc, lt = 10 ,10, 1, 1 
    bp = [[0, 0], [width, 0], [width, height], [0, height]]
    material = Engproblem.Material(lc, lt)
    domain = Engproblem.Domain(bp)
    loads = Engproblem.Load([width/2, height], [1, 0], 1)
    supports = [Engproblem.Support([0, 0], [0, 0]), Engproblem.Support([width, 0], [1, 0])]
    
##Solve Optimization problem
    fl, step = 1.42, 1
    groundstructure = Optproblem.Groundstructure(domain, loads, supports, material, step, step, fl)
    prob = Optproblem.Optprob(groundstructure)
    
##Solve Mathematic problem
    vol = prob.solve()
    print("final vol:%f" % (vol))
