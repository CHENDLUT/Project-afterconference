from Engproblem import *
from Mathproblem import *
import numpy as np
import itertools
from shapely.geometry import Polygon, Point
from math import gcd,ceil
import matplotlib.pyplot as plt
from scipy import sparse

class Node:
    def __init__(self, number, position):
        self.number = number
        self.position = position
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]
        self.load = []
        self.support = Support(position, [1,1,1])
    def updatesupport(self, supports):
        for sp in supports:
            if (self.position[0] == sp.position[0] and self.position[1] == sp.position[1] and self.position[2] == sp.position[2]):
                self.support = sp 
    def updateload(self, loads):
        for ld in loads:
            if (self.position[0] == ld.position[0] and self.position[1] == ld.position[1] and self.position[2] == ld.position[2]):
                self.load.append(ld) 

class Member:
    def __init__(self, start, end, material):
        self.start = start
        self.end = end
        self.material = material
        self.area = 0.0
        self.length = self.callength()
        self.component = self.calcomponent()
        self.exist = False
    def calcomponent(self):
        return (np.array(self.end.position)-np.array(self.start.position))
    def callength(self):
        return np.sqrt(np.sum([(np.array(self.end.position)-np.array(self.start.position))[i]**2 for i in range(3)]))
    def updateexist(self, bl):
        self.exist = bl
    def updatearea(self, a):
        self.area = a 

class Groundstructure:
    def __init__(self, domain, supports, loads, step1, step2, step3, fl, material):
        self.domain = domain
        self.supports = supports
        self.loads = loads
        self.step1 = step1
        self.step2 = step2
        self.step3 = step3
        self.fl = fl
        self.material = material
        self.NDL = []
        self.ML = []
        self.TML = []
        self.PML = []
    def createnodes(self):# uncertain and can be replaced
        xv, yv, zv = np.meshgrid(np.arange(0,self.domain.length+self.step1,self.step1), np.arange(0,self.domain.width+self.step2,self.step2), np.arange(0,self.domain.height+self.step3,self.step3))
        pts = [Point(xv.flat[i], yv.flat[i], zv.flat[i]) for i in range(xv.size)]
        i = 0
        for pl in self.domain.polylist:
            poly = Polygon(pl)
            for pt in pts:
                if poly.intersects(pt):
                    self.NDL.append(Node(i,[pt.x,pt.y,pt.z]))
                    i+=1
        for nd in self.NDL:
            nd.updatesupport(self.supports)
        for nd in self.NDL:
            nd.updateload(self.loads)
    def createmembers(self):
        for i,j in itertools.combinations(range(len(self.NDL)),2):
            dx, dy, dz = abs(self.NDL[i].x - self.NDL[j].x), abs(self.NDL[i].y - self.NDL[j].y), abs(self.NDL[i].z-self.NDL[j].z)
            if gcd(int(dx/self.step1), int(dy/self.step2), int(dz/self.step3)) == 1:
                if self.domain.convex:
                    self.ML.append(Member(self.NDL[i], self.NDL[j], self.material))
        for mem in self.ML:
            if mem.length<=self.fl:mem.exist=True
    def createpotentials(self):
        self.TML,self.PML = [],[]
        for mem in self.ML:
            if mem.exist==True:
                self.TML.append(mem)
            else:
                self.PML.append(mem)
        
class Optprob:
    def __init__(self, domain, supports, loads, material, step1, step2, step3, fl):
        self.groundstructure = Groundstructure(domain, supports, loads, step1, step2, step3, fl, material)
        self.mathprob = None
        self.B = None
        self.judgeB = None
        self.f = []
        self.dof = []
        self.con = True
    def initial(self):
        self.groundstructure.createnodes()
        self.groundstructure.createmembers()
        self.groundstructure.createpotentials()
    def processBoundaryCondition(self):
        for nd in self.groundstructure.NDL:
            self.dof.append([nd.support.condition])
        ldcase = max([len(node.load) for node in self.groundstructure.NDL])
        for i in range(ldcase):
            self.f.append([nd.load[i].magnitude if i<=len(nd.load)-1 else [0,0,0] for nd in self.groundstructure.NDL])
            self.f[i] = np.array(self.f[i]).flatten()
        self.f, self.dof = np.array(self.f), np.array(self.dof).flatten()
    def MatrixB(self):
        NDL, TML = self.groundstructure.NDL, self.groundstructure.TML
        m, n1, n2 = len(TML), np.array([int(mem.start.number) for mem in TML]), np.array([int(mem.end.number) for mem in TML])
        Nd = np.array([NDL[i].position for i in range(len(NDL))])
        l, X, Y, Z = np.array([mem.length for mem in TML]), Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1], Nd[n2,2]-Nd[n1,2]
        d0, d1, d2, d3, d4, d5 = self.dof[n1*3], self.dof[n1*3+1], self.dof[n1*3+2], self.dof[n2*3], self.dof[n2*3+1], self.dof[n2*3+2]
        s = np.concatenate((-X/l * d0, -Y/l * d1, -Z/l*d2, X/l * d3, Y/l * d4, Z/l*d5))
        r = np.concatenate((n1*3, n1*3+1, n1*3+2, n2*3, n2*3+1, n2*3+2))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        self.B = sparse.coo_matrix((s, (r, c)), shape = (len(NDL)*3, m))
    def judgeMatrixB(self):
        NDL, PML = self.groundstructure.NDL, self.groundstructure.PML
        m, n1, n2 = len(PML), np.array([int(mem.start.number) for mem in PML]), np.array([int(mem.end.number) for mem in PML])
        Nd = np.array([NDL[i].position for i in range(len(NDL))])
        l, X, Y, Z = np.array([mem.length for mem in PML]), Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1], Nd[n2,2]-Nd[n1,2]
        d0, d1, d2, d3, d4, d5 = self.dof[n1*3], self.dof[n1*3+1], self.dof[n1*3+2], self.dof[n2*3], self.dof[n2*3+1], self.dof[n2*3+2]
        s = np.concatenate((-X/l * d0, -Y/l * d1, -Z/l*d2, X/l * d3, Y/l * d4, Z/l*d5))
        r = np.concatenate((n1*3, n1*3+1, n1*3+2, n2*3, n2*3+1, n2*3+2))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m),np.arange(m)))
        self.judgeB = (sparse.coo_matrix((s, (r, c)), shape = (len(NDL)*3, m))).tocsc()
    def solve(self):
        for itr in range(1,100):
            if self.con == True:
                self.MatrixB()
                self.mathprob = Mathprob(self.B)
                self.mathprob.solveLP(self.groundstructure.NDL, self.groundstructure.TML, self.dof, self.f, self.groundstructure.material.lc, self.groundstructure.material.lt)
                self.plot("Itr:" + str(itr),True)
                print("Itr: %d,vol: %f,mems:%d" % (itr, self.mathprob.result, len(self.groundstructure.TML)))
                self.judgeMatrixB()
                self.stopviolation()
                self.groundstructure.createpotentials()
            else:
                break
        print("Final Volume: %f" % (self.mathprob.result))       
    def stopviolation(self):
        l, lst = [], []
        for i,mem in enumerate(self.groundstructure.ML):
            if (mem.exist == False):
                l.append(mem.length)
                lst.append(i)
        l,lst = np.array(l), np.array(lst)
        y = np.zeros(len(self.groundstructure.PML)) 
        for uk in self.mathprob.ulist:
            yk = np.multiply(self.judgeB.transpose().dot(uk)/l, np.array([[self.groundstructure.material.lt], [-self.groundstructure.material.lc]]))
            y += np.amax(yk, axis=0)
        vioPML = np.where(y>1.0001)[0]
        vioSort = np.flipud(np.argsort(y[vioPML]))
        num = ceil(min(len(vioSort), 0.05*max([len(self.groundstructure.PML)*0.05, len(vioSort)])))
        for i in range(num):
            self.groundstructure.ML[lst[vioPML[vioSort[i]]]].exist = True
        self.con = (num != 0)
    def updatedata(self):
        for i in range(len(self.mathprob.arealist)):
            self.groundstructure.TML[i].updatearea(self.mathprob.arealist[i])
    def plot(self, str, update):
        mp = self.mathprob
        plt.ion() if update else plt.ioff()
        plt.clf(); plt.axis('off'); plt.draw()
        ax1 = plt.axes(projection='3d')
        plt.title(str)
        tk = 3/max(mp.arealist)
        for i in [i for i in range(len(mp.arealist)) if mp.arealist[i] >= max(mp.arealist) * 1e-3]:
            if all([mp.qlist[lc][i]>=0 for lc in range(len(mp.qlist))]): c = 'r'
            elif all([mp.qlist[lc][i]<=0 for lc in range(len(mp.qlist))]): c = 'b'
            else: c = 'tab:gray'
            plt.plot([self.groundstructure.TML[i].start.position[0],self.groundstructure.TML[i].end.position[0]],
                     [self.groundstructure.TML[i].start.position[1],self.groundstructure.TML[i].end.position[1]],
                     [self.groundstructure.TML[i].start.position[2],self.groundstructure.TML[i].end.position[2]],c, linewidth = mp.arealist[i] * tk)
        plt.show()    
