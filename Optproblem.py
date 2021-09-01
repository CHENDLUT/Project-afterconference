import Engproblem
import numpy as np
from scipy import sparse
import itertools
from math import gcd, ceil
from shapely.geometry import Point, LineString, Polygon
import matplotlib.pyplot as plt
import Mathproblem

class Node:
    def __init__(self, nd):
        self.nd = nd
        self.x = nd[1]
        self.y = nd[2]

class Member:
    def __init__(self, start, end, material, exist):
        self.start = start
        self.end = end
        self.material = material
        self.exist = exist
    def calcomponent(self):
        component = self.end - self.start
        return [component[1], component[2]]
    def callength(self):
        cp = self.calcomponent()
        return (np.sqrt(np.sum([cp[m]**2 for m in range(2)])))
    def direction(self):
        cp = self.calcomponent()
        return [cp[0]/self.callength(), cp[1]/self.callength()]
    def crevector(self):
        return [self.end[0], self.start[0], selfself.callength()]

class Groundstructure:
    def __init__(self, domain, loads, supports, material, step1, step2, fl):
        self.domain = domain
        self.loads = loads
        self.supports = supports
        self.material = material
        self.step1 = step1
        self.step2 = step2
        self.fl = fl
    def createpoints(self):
        width, height = self.domain.size()
        step1, step2 = self.step1, self.step2
        xv, yv = np.meshgrid(range(0, width+step1, step1), range(0, height+step2, step2))
        pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]
        points = np.array([[i, pt.x, pt.y] for i,pt in enumerate(pts) if self.domain.poly.intersects(pt)])
        return points
    def createmembers(self):
        points = self.createpoints()
        ML = []
        for i,j in itertools.combinations(range(len(points)), 2):
            mem = Member(points[j], points[i], self.material, False)
            if gcd(int(abs(mem.calcomponent()[0])), int(abs(mem.calcomponent()[1]))) == 1:
                ML.append(Member(points[j], points[i], self.material, False))
        return ML
    def createGSML(self):
        ML = self.createmembers()
        for mem in ML:
            if (mem.callength()<=self.fl): mem.exist = True
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
            if (pt[1] == ldp[0] and pt[2] == ldp[1]) :
                LDV[2*i], LDV[2*i+1] = ld.magnitude[0], ld.magnitude[1]
        return LDV
    def createMatrixB(self, points, TML, SPV):
        m, n1, n2 = len(TML), np.array([int(mem.start[0]) for mem in TML]), np.array([int(mem.end[0]) for mem in TML])
        l, X, Y = np.array([mem.callength() for mem in TML]), np.array([mem.calcomponent()[0] for mem in TML]), np.array([mem.calcomponent()[1] for mem in TML])
        d0, d1, d2, d3 = SPV[n1*2], SPV[n1*2+1], SPV[n2*2], SPV[n2*2+1]
        s = np.concatenate((-X/l * d0, -Y/l * d1, X/l * d2, Y/l * d3))
        r = np.concatenate((n1*2, n1*2+1, n2*2, n2*2+1))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        return sparse.coo_matrix((s, (r, c)), shape = (len(points)*2, m))
    def createLV(self):
        ML = self.groundstructure.createmembers()
        return np.array([mem.callength() for mem in ML])
    def createGSMLV(self):
        GSML, GSMLV = self.groundstructure.createGSML(), []
        for gsml in GSML:
            if (gsml.exist == True):GSMLV.append(gsml.callength())
        return GSMLV
    def solve(self):
        points, ML, SPV, LDV, LV = self.groundstructure.createpoints(), self.groundstructure.createGSML(), self.createSPV(), self.createLDV(), self.createGSMLV()
        for itr in range(1, 100):
            TML, lst, LV = [], [], []
            for i,mem in enumerate(ML):
                if (mem.exist == True):
                    TML.append(mem)
                    LV.append(mem.callength())
                else: lst.append(i)
            LV = np.array(LV)
            B = self.createMatrixB(points, TML, SPV)
            mathprob = Mathproblem.Mathproblem(LV, B, LDV, SPV, self.groundstructure.material.lc, self.groundstructure.material.lt)
            vol, a, q, u = mathprob.solveLP()
            self.plot(q, a, max(a) * 1e-3, TML, "Itr:" + str(itr))
            print("Itr: %d, vol: %f" % (itr, vol))
            bl, ML = self.stopviolation(u, ML, lst, points, SPV)
            if bl:break
        self.plot(q, a, max(a) * 1e-3, TML, "Finished", False)
        return vol
    def stopviolation(self, u, ML, lst, points, SPV):
        FML, lst = [], np.array(lst)
        for mem in ML:
            if(mem.exist == False):
                FML.append(mem)
        l = np.array([fml.callength() for fml in FML])
        B = self.createMatrixB(points, FML, SPV).tocsc()
        y = np.zeros(len(FML))
        for uk in u:
            yk = np.multiply(B.transpose().dot(uk) / l, np.array([[self.groundstructure.material.lt], [-self.groundstructure.material.lc]]))
            y += np.amax(yk, axis=0)
        vioFML = np.where(y>1.0001)[0]
        vioSort = np.flipud(np.argsort(y[vioFML]))
        num = ceil(min(len(vioSort), 0.05*max( [len(FML)*0.05, len(vioSort)])))
        for i in range(num):
            ML[lst[vioFML[vioSort[i]]]].exist = True
        return num == 0, ML
    def plot(self, q, a, threshold, TML, string, update = True):
        plt.ion() if update else plt.ioff()
        plt.clf(); plt.axis('off'); plt.axis('equal');  plt.draw()
        plt.title(string)
        tk = 5 / max(a)
        for i in [i for i in range(len(a)) if a[i] >= threshold]:
            if all([q[lc][i]>=0 for lc in range(len(q))]): c = 'r'
            elif all([q[lc][i]<=0 for lc in range(len(q))]): c = 'b'
            else: c = 'tab:gray'
            plt.plot([TML[i].start[1], TML[i].end[1]], [TML[i].start[2], TML[i].end[2]], c, linewidth = a[i] * tk)
        plt.pause(0.01) if update else plt.show()
