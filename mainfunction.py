from Engprob import *
from Optprob import *
from Mathprob import *

##Solve engineering problem
w, h = 10.0, 10.0 
vt = [[0,0],[w,0],[w,h],[0,h]]
step1, step2, lc, lt = 0.5, 0.5, 1.0, 1.0
material = Material(lc, lt) # define the material
domain = Domain(vt, material) # define design domain
domain.verify() # determine the convexity
supports = [Support([0,0],[0,0]),Support([w,0], [0,0])] # set support position, support condition
loads = [Load([w/2,h],[1.0,0.0],1)] # set load position, load magnitude
preml = [Member(Node(0,[0.0,0.0]),Node(0,[5.0,1.5]),material), Member(Node(0,[10.0,0.0]),Node(0,[5.0,1.5]),material)]

##Solve optimization problem
prob = Optprob(domain, supports, loads, material, preml, step1, step2) # define optimization problem
prob.initial() # initialize optimization problem
prob.processBoundaryCondition() # link nodes to loads and supports

##Solve math problem
prob.solve() # solve the problem
prob.updatedata()
print("Optimized Volume is: %.2f"%(prob.mathprob.result))
prob.plot()

