from Engproblem import *
from Optproblem import *
from Mathproblem import *
from shapely.geometry import Polygon
import numpy as np
if __name__ == '__main__':
    ## Solve engineering problem
    l, w, h = 10.0, 10.0, 10.0
    lc, lt ,step1, step2, step3, prearea, jc = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0 
    fl = 1.001*np.sqrt(step1**2+step2**2+step3**2)
    vt = [[0,0,0],[l,0,0],[l,w,0],[0,w,0],[0,0,h],[l,0,h],[l,w,h],[0,w,h]] #set vertex coordinates
    material = Material(lc, lt)
    domain = Domain(vt) # create the design domain
    domain.polyhedron.verify() # verify the convexity of the design domain
    loads = [Load([l/2,w/2,0],[0,0,-1],1),Load([5.2,3.3,3.1],[0,0,-1],1)]
    supports = [Support([0,0,0],[0,0,0]),Support([l,0,0],[1,1,0]),Support([l,w,0],[1,1,0]),Support([0,w,0],[0,1,0])]
    Pre_member = [Member(Node(-1,[0.0,0.0,0.0]),Node(-1,[3.0,1.0,1.0]),material),
                  Member(Node(-1,[4.0,4.0,4.0]),Node(-1,[8.0,5.0,6.0]),material),
                  Member(Node(-1,[3.0,1.0,1.0]),Node(-1,[5.2,3.3,3.1]),material)] # input the pre-existing members and the default initial node number is negative
    ## Solve optimization problem
    optprob = Optprob(domain, supports, loads, Pre_member, prearea, material, step1, step2, step3, fl, jc) # create a optimization problem
    optprob.initial() # initialize the optimization problem including creating the groundstructure, creating nodes, creating members and potential members
    optprob.processBoundaryCondition() # link nodes to loads and supports 

    ## solve() will create constraints and variables for the Mathematical problem; 
    ## and then Mathematical problem will pass the necessary information to the solver.
    
    ## cvxpy, pyomo are useful interface libraries for mathematical programming problems.
    optprob.solve()
    optprob.updatedata() # pass the area value to each member
    optprob.plot("Finished", False)
    
