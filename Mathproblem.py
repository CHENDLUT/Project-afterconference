import cvxpy as cvx
import numpy as np
class Mathprob:
    def __init__(self, B):
        self.B = B
        self.lengthlist = []
        self.qlist = []
        self.ulist = []
        self.arealist = []
        self.result = 0
    def solveLP(self, NDL, TML, dof, f, lc, lt):
        self.lengthlist = [ml.length for ml in TML]
        a = cvx.Variable(len(TML))
        obj = cvx.Minimize(self.lengthlist * a)
        q, eqn, cons= [], [], [a>=0]
        for k, fk in enumerate(f):
            q.append(cvx.Variable(len(TML)))
            eqn.append(self.B * q[k] == fk * dof)
            cons.extend([eqn[k], q[k] >= -lc * a, q[k] <= lt * a])
        prob = cvx.Problem(obj, cons)
        self.result = prob.solve()
        self.qlist = [np.array(qi.value).flatten() for qi in q]
        self.arealist = np.array(a.value).flatten()
        self.ulist = [-np.array(eqnk.dual_value).flatten() for eqnk in eqn]
