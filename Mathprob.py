from Engprob import *
from Mathprob import *
import numpy as np
import cvxpy as cvx

class Mathproblem:
    def __init__(self, B):
        self.B = B
        self.lengthlist = []
        self.arealist = []
        self.qlist = []
        self.ulist = []
        self.result = 0
    def solveLP(self, NDL, ML, dof, f, lc,lt):
        self.lengthlist = [ml.length for ml in ML]
        a = cvx.Variable(len(ML))
        obj = cvx.Minimize(self.lengthlist * a)
        q, eqn, cons= [], [], [a>=0]
        for k, fk in enumerate(f):
            q.append(cvx.Variable(len(ML)))
            eqn.append(self.B * q[k] == fk * dof)
            cons.extend([eqn[k], q[k] >= -lc * a, q[k] <= lt * a])
        prob = cvx.Problem(obj, cons)
        self.result = prob.solve()
        self.qlist = [np.array(qi.value).flatten() for qi in q]
        self.arealist = np.array(a.value).flatten()
        self.ulist = [-np.array(eqnk.dual_value).flatten() for eqnk in eqn]
