from Engproblem import *
from Optproblem import *
from Mathproblem import *
from shapely.geometry import Polygon
if __name__ == '__main__':
    ## Solve engineering problem
    l, w, h = 20.0, 10.0, 10.0
    lc, lt ,step1, step2, step3, fl = 1.0, 1.0, 1.0, 1.0, 1.0, 1.74
    vt = [[0,0,0],[l,0,0],[l,w,0],[0,w,0],[0,0,h],[l,0,h],[l,w,h],[0,w,h]]
    area = 2*(w*h+h*l+w*l) # uncertain and can be replaced by vtk
    material = Material(lc, lt)
    domain = Domain(vt, step1)
    domain.verify(area)
    loads = [Load([l/2,w/2,0],[0,0,-1],1)]
    supports = [Support([0,0,0],[0,0,0]),Support([l,0,0],[1,1,0]),Support([l,w,0],[1,1,0]),Support([0,w,0],[0,1,0])]
    ## Solve optimization problem
    optprob = Optprob(domain, supports, loads, material, step1, step2, step3, fl)
    optprob.initial()
    optprob.processBoundaryCondition()
    ## Solve mathematical problem
    optprob.solve()
    optprob.updatedata()
    optprob.plot("Finished", False)
    
