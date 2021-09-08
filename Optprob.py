from Engprob import *
from Mathprob import *
import numpy as np
from math import gcd
import itertools
from scipy import sparse
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString

class Node:
    def __init__(self, number, position):
        self.number = number
        self.position = position
        self.x = position[0]
        self.y = position[1]
        self.load = []
        self.support = Support(position, [1,1])
    def updatesupport(self, supports):
        for sp in supports:
            if (self.position[0] == sp.position[0] and self.position[1] == sp.position[1]):
                self.support = sp  
    def updateload(self, loads):
        for ld in loads:
            if (self.position[0] == ld.position[0] and self.position[1] == ld.position[1]):
                self.load.append(ld) 
    
class Member:
    def __init__(self, start, end, material):
        self.start = start
        self.end = end
        self.material = material
        self.area = 0
        self.length = self.callength()
        self.component = self.calcomponent()
    def calcomponent(self):
        return (np.array(self.end.position)-np.array(self.start.position))
    def callength(self):
        return np.sqrt(np.sum([(np.array(self.end.position)-np.array(self.start.position))[i]**2 for i in range(2)]))
    def updatearea(self, a):
        self.area = a

class Groundstructure:
    def __init__(self, domain, supports, loads, step1, step2, material):
        self.domain = domain
        self.supports = supports
        self.loads = loads
        self.step1 = step1
        self.step2 = step2
        self.material = material
        self.NDL = []
        self.ML = []
    def createnodes(self):
        xv, yv = np.meshgrid(np.arange(0,self.domain.width+self.step1,self.step1), np.arange(0,self.domain.height+self.step2,self.step2))
        pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]
        for i in range(len(pts)):
            self.NDL.append(Node(i, [pts[i].x,pts[i].y]))
        for nd in self.NDL:
            nd.updatesupport(self.supports)
        for nd in self.NDL:
            nd.updateload(self.loads) 
    def createmembers(self):
        for i, j in itertools.combinations(range(len(self.NDL)), 2):
            member = Member(self.NDL[i],self.NDL[j],self.material)
            if gcd(int(member.component[0]/self.step1), int(member.component[1]/self.step2)) == 1:
                seg = [] if self.domain.convex else LineString([NDL[i].position, NDL[j].position])
                if self.domain.convex or self.domain.poly.contains(seg) or self.domain.poly.boundary.contains(seg):
                    self.ML.append(member)
    
            
class Optprob:
    def __init__(self, domain, supports, loads, material, step1, step2):
        self.groundstructure = Groundstructure(domain, supports, loads, step1, step2, material)
        self.mathprob = None
        self.B = None
        self.f = []
        self.dof = []
    def initial(self):
        self.groundstructure.createnodes()
        self.groundstructure.createmembers()
    def MatrixB(self):
        NDL, ML = self.groundstructure.NDL, self.groundstructure.ML
        m, n1, n2 = len(ML), np.array([int(mem.start.number) for mem in ML]), np.array([int(mem.end.number) for mem in ML])
        Nd = np.array([NDL[i].position for i in range(len(NDL))])
        l, X, Y = np.array([mem.length for mem in ML]), Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1]
        d0, d1, d2, d3 = self.dof[n1*2], self.dof[n1*2+1], self.dof[n2*2], self.dof[n2*2+1]
        s = np.concatenate((-X/l * d0, -Y/l * d1, X/l * d2, Y/l * d3))
        r = np.concatenate((n1*2, n1*2+1, n2*2, n2*2+1))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        self.B = sparse.coo_matrix((s, (r, c)), shape = (len(NDL)*2, m))
    def processBoundaryCondition(self):
        for nd in self.groundstructure.NDL:
            self.dof.append([nd.support.condition])
        ldcase = max([len(node.load) for node in self.groundstructure.NDL])
        for i in range(ldcase):
            self.f.append([nd.load[i].magnitude if i<=len(nd.load)-1 else [0,0] for nd in self.groundstructure.NDL])
            self.f[i] = np.array(self.f[i]).flatten()
        self.f, self.dof = np.array(self.f), np.array(self.dof).flatten()
    def solve(self):
        self.MatrixB()
        self.mathprob = Mathproblem(self.B)
        self.mathprob.solveLP(self.groundstructure.NDL, self.groundstructure.ML, self.dof, self.f, self.groundstructure.material.lc, self.groundstructure.material.lt)  
    def plot(self):
        mp = self.mathprob
        plt.ion(), plt.clf(); plt.axis('off'); plt.axis('equal'); plt.draw(); plt.title("Final result")
        tk = 5 / max(mp.arealist)
        for i in [i for i in range(len(mp.arealist)) if mp.arealist[i] >= max(mp.arealist) * 1e-3]:
            if all([mp.qlist[lc][i]>=0 for lc in range(len(mp.qlist))]): c = 'r'
            elif all([mp.qlist[lc][i]<=0 for lc in range(len(mp.qlist))]): c = 'b'
            else: c = 'tab:gray'
            plt.plot([self.groundstructure.ML[i].start.position[0],self.groundstructure.ML[i].end.position[0]], [self.groundstructure.ML[i].start.position[1],self.groundstructure.ML[i].end.position[1]], c, linewidth = mp.arealist[i] * tk)
        plt.show()
    def updatedata(self):
        for i in range(len(self.mathprob.arealist)):
            self.groundstructure.ML[i].updatearea(self.mathprob.arealist[i])
            
     
