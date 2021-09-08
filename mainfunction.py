from Engproblem import *
from Optproblem import *
from Mathproblem import *
from shapely.geometry import Polygon
if __name__ == '__main__':
    ## Solve engineering problem
    w, h = 10.0, 10.0
    lc, lt ,step1, step2, fl = 1.0, 1.0, 1.0, 1.0, 1.42
    vt = [[0,0],[w,0],[w,h],[0,h]]
    material = Material(lc, lt)
    domain = Domain(vt)
    domain.verify()
    loads = [Load([w/2,h],[1,0],1)]
    supports = [Support([0,0],[0,0]),Support([w,0],[0,0])]
    ## Solve optimization problem
    optprob = Optprob(domain, supports, loads, material, step1, step2, fl)
    optprob.initial()
    optprob.processBoundaryCondition()
    ## Solve mathematical problem
    optprob.solve()
    optprob.updatedata()
    optprob.plot("Finished", False)
    
