import Engproblem
import numpy as np
from scipy import sparse
import itertools
from math import gcd, ceil
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon

class Node:
    def __init__(self, nd):
        self.nd = nd
        self.x = nd[1]
        self.y = nd[2]

class Member:
    def __init__(self, start, end, material):
        self.start = start
        self.end = end
        self.material = material
    def calcomponent(self):
        component = self.end - self.start
        return [component[1], component[2]]
    def callength(self):
        cp = self.calcomponent()
        return (np.sqrt(np.sum([cp[m]**2 for m in range(2)])))
    def direction(self):
        cp = self.end - self.start
        return [cp[1]/self.callength(), cp[2]/self.callength()]
    def crevector(self):
        return [self.end[0], self.start[0], selfself.callength()]

class Groundstructure:
    def __init__(self, domain, loads, supports, material, step1, step2):
        self.domain = domain
        self.loads = loads
        self.supports = supports
        self.material = material
        self.step1 = step1
        self.step2 = step2
    def createpoints(self):
        width, height = self.domain.size()
        step1, step2 = self.step1, self.step2
        xv, yv = np.meshgrid(np.arange(0, width+self.step1, self.step1), np.arange(0, height+self.step2, self.step2))
        pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]
        points = np.array([[i, pt.x, pt.y] for i,pt in enumerate(pts) if self.domain.poly.intersects(pt)])
        return points
    def createmembers(self):
        points = self.createpoints()
        ML = []
        for i,j in itertools.combinations(range(len(points)), 2):
            mem = Member(points[j], points[i], self.material)
            if gcd(int(abs(int(mem.calcomponent()[0]/self.step1))), int(abs(int(mem.calcomponent()[1]/self.step2)))) == 1:
                ML.append(mem)
        return ML
    
class Optprob:
    def __init__(self, groundstructure):
        self.groundstructure = groundstructure
        self.domain = groundstructure.domain
        self.loads = groundstructure.loads
        self.supports = groundstructure.supports
        self.material = groundstructure.material
        self.step1 = groundstructure.step1
        self.step2 = groundstructure.step2
    def createSPV(self):
        points = self.groundstructure.createpoints()
        SPV = np.ones(2*len(points))
        for sp in self.groundstructure.supports:
            for i,pt in enumerate(points):
                if (pt[1] == sp.position[0] and pt[2] == sp.position[1]):
                    SPV[2*i], SPV[2*i+1] = sp.condition[0], sp.condition[1]
        return SPV    
    def createLDV(self):
        points = self.groundstructure.createpoints()
        LDV = np.zeros(2*len(points))
        ld, ldp = self.groundstructure.loads, self.groundstructure.loads.position
        for i,pt in enumerate(points):
            if (pt[1] == ldp[0] and pt[2] == ldp[1]): LDV[2*i], LDV[2*i+1] = ld.magnitude[0], ld.magnitude[1]
        return LDV
    def createMatrixB(self, SPV):
        points, ML, SPV = self.groundstructure.createpoints(), self.groundstructure.createmembers(), self.createSPV()
        m, n1, n2 = len(ML), np.array([int(mem.start[0]) for mem in ML]), np.array([int(mem.end[0]) for mem in ML])
        l, X, Y = np.array([mem.callength() for mem in ML]), np.array([mem.calcomponent()[0] for mem in ML]), np.array([mem.calcomponent()[1] for mem in ML])
        d0, d1, d2, d3 = SPV[n1*2], SPV[n1*2+1], SPV[n2*2], SPV[n2*2+1]
        s = np.concatenate((-X/l * d0, -Y/l * d1, X/l * d2, Y/l * d3))
        r = np.concatenate((n1*2, n1*2+1, n2*2, n2*2+1))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        return sparse.coo_matrix((s, (r, c)), shape = (len(points)*2, m))
    def createLV(self):
        ML = self.groundstructure.createmembers()
        return np.array([mem.callength() for mem in ML])
    def plot(self, q, a, threshold, string):
        plt.ion(); plt.clf(); plt.axis('off'); plt.axis('equal');  plt.draw()
        plt.title(string)
        ML = self.groundstructure.createmembers()
        tk = 5 / max(a)
        for i in [i for i in range(len(a)) if a[i] >= threshold]:
            if all([q[lc][i]>=0 for lc in range(len(q))]): c = 'r'
            elif all([q[lc][i]<=0 for lc in range(len(q))]): c = 'b'
            else: c = 'tab:gray'
            plt.plot([ML[i].start[1], ML[i].end[1]], [ML[i].start[2], ML[i].end[2]], c, linewidth = a[i] * tk)
        plt.show()
        
        
