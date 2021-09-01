import numpy as np
import cvxpy as cvx
import Engproblem
import Optproblem

class Mathproblem:
    def __init__(self, LV, B, LDV, SPV, lc, lt):
        self.LV = LV
        self.B = B
        self.LDV = LDV
        self.SPV = SPV
        self.lc = lc
        self.lt = lt
    def solveLP(self):
        a = cvx.Variable(len(self.LV))
        obj = cvx.Minimize(np.transpose(self.LV) * a)
        LDV = [self.LDV]
        q, eqn, cons = [], [], [a>=0]
        for k, fk in enumerate(LDV):
            q.append(cvx.Variable(len(self.LV)))
            eqn.append(self.B * q[k] == fk * self.SPV)
            cons.extend([eqn[k], q[k] >= -self.lc * a, q[k] <= self.lt * a])
        prob = cvx.Problem(obj, cons)
        vol = prob.solve()
        q = [np.array(qi.value).flatten() for qi in q]
        a = np.array(a.value).flatten()
        u = [-np.array(eqnk.dual_value).flatten() for eqnk in eqn]
        return vol, a, q, u
