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
    def solve(self):
        a = cvx.Variable(len(self.LV))
        obj = cvx.Minimize(self.LV * a)
        q, eqn, cons = [], [], [a>=0]
        LDV = [self.LDV]
        for k, fk in enumerate(LDV):
            q.append(cvx.Variable(len(self.LV)))
            eqn.append(self.B * q[k] == fk * self.SPV)
            cons.extend([eqn[k], q[k] >= -self.lc * a, q[k] <= self.lt * a])
        prob = cvx.Problem(obj, cons)
        vol = prob.solve()
        q = [np.array(qi.value).flatten() for qi in q]
        a = np.array(a.value).flatten()
        return vol, q, a
