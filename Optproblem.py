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
    def updatenode(self, supports, loads): # update the condition of load and support of the node
        for sp in supports:
            if (np.sqrt(np.sum([(np.array(self.position)-np.array(sp.position))[i]**2 for i in range(3)]))<=1e-6):
                self.support = sp
        for ld in loads:
            if (np.sqrt(np.sum([(np.array(self.position)-np.array(ld.position))[i]**2 for i in range(3)]))<=1e-6):
                self.load.append(ld)   

class Member:
    def __init__(self, start, end, material):
        self.start = start
        self.end = end
        self.material = material
        self.area = 100.0
        self.length = self.callength()
        self.component = self.calcomponent()
        self.exist = False
    def calcomponent(self): # calculate the x,y,z components of the member
        return (np.array(self.end.position)-np.array(self.start.position))
    def callength(self): # calculate the length of the member
        return np.sqrt(np.sum([(np.array(self.end.position)-np.array(self.start.position))[i]**2 for i in range(3)]))
    def contain(self, nd): # verify whether this member contains the specific node
        if self.start.position == nd.position or self.end.position == nd.position: return False
        else:
            c1,c2 = np.array(nd.position)-np.array(self.start.position), np.array(self.end.position)-np.array(nd.position)
            L1,L2 = np.sqrt(np.sum([(np.array(self.start.position)-np.array(nd.position))[i]**2 for i in range(3)])), np.sqrt(np.sum([(np.array(self.end.position)-np.array(nd.position))[i]**2 for i in range(3)]))
            direction1, direction2 = [c1[0]/L1,c1[1]/L1,c1[2]/L1], [c2[0]/L2,c2[1]/L2,c2[2]/L2]
            return True if (direction1[0]-direction2[0]<=1e-6 and direction1[1]-direction2[1]<=1e-6 and direction1[2]-direction2[2]<=1e-6) else False
    def updatearea(self, a): # update the area of the member
        self.area = a

class Groundstructure:
    def __init__(self, domain, supports, loads, Pre_member, prearea, step1, step2, step3, fl, jc, material):
        self.domain = domain
        self.supports = supports
        self.loads = loads
        self.Pre_member = Pre_member
        self.prearea = prearea
        self.step1 = step1
        self.step2 = step2
        self.step3 = step3
        self.fl = fl
        self.jc = jc
        self.material = material
        self.NDL = []
        self.ML = []
        self.TML = []
        self.PML = []
        self.prenumber = []
    def createnodes(self): # create the nodes of grid
        xv, yv, zv = np.meshgrid(np.arange(0,self.domain.length+self.step1,self.step1),\
                                 np.arange(0,self.domain.width+self.step2,self.step2),\
                                 np.arange(0,self.domain.height+self.step3,self.step3))
        pts = [Point(xv.flat[i], yv.flat[i], zv.flat[i]) for i in range(xv.size)]
        for pt in pts:
            if self.domain.polyhedron.intersect(pt):
                self.NDL.append(Node(len(self.NDL), [pt.x, pt.y, pt.z]))
        for nd in self.NDL:
            nd.updatenode(self.supports, self.loads)
    def createmembers(self): # create the members based on the grid
        for i, j in itertools.combinations(range(len(self.NDL)),2):
            mem = Member(self.NDL[j], self.NDL[i], self.material)
            if gcd(int(mem.component[0]/self.step1), int(mem.component[1]/self.step2), int(mem.component[2]/self.step3)) == 1 or self.jc!=0:
                self.ML.append(mem)
        for mem in self.ML:
            if mem.length <= self.fl: mem.exist = True
    def addpremember(self): # add pre-existing nodes and members to the node list and member list
        for m,pm in enumerate(self.Pre_member):
            pm.exist = True
            pm.area = self.prearea
            for nd in self.NDL:
                if np.sqrt(np.sum([(np.array(pm.start.position)-np.array(nd.position))[i]**2 for i in range(3)]))<=1e-6:pm.start.number = nd.number
                elif np.sqrt(np.sum([(np.array(pm.end.position)-np.array(nd.position))[i]**2 for i in range(3)]))<=1e-6:pm.end.number = nd.number
            if pm.start.number<0:
                pm.start.number = len(self.NDL)
                self.NDL.append(pm.start)
                for nd in self.NDL:
                    mem = Member(pm.start, nd, self.material)
                    if mem.contain(nd) == False:self.ML.append(mem)
            elif pm.end.number<0:
                pm.end.number = len(self.NDL)
                self.NDL.append(pm.end)
                for nd in self.NDL:
                    mem = Member(pm.end, nd, self.material)
                    if mem.contain(nd) == False or self.jc!=0:self.ML.append(mem)
    def createpotentials(self): # update the exist of members in member list and sort out the initial members and potential members
        self.TML, self.PML = [], []
        for pm in self.Pre_member:
            for i,mem in enumerate(self.ML):
                if (pm.start.number == mem.start.number and pm.end.number == mem.end.number) or (pm.end.number == mem.start.number and pm.start.number == mem.end.number):self.ML[i] = pm
        for mem in self.ML:
            if mem.exist == True: self.TML.append(mem)
            else: self.PML.append(mem)
        self.prenumber = []
        for pm in self.Pre_member:
            for i, tml in enumerate(self.TML):
                if pm.start.number == tml.start.number and pm.end.number == tml.end.number:
                    self.prenumber.append(i)

class Optprob:
    def __init__(self, domain, supports, loads, Pre_member, prearea, material, step1, step2, step3, fl, jc):
        self.groundstructure = Groundstructure(domain, supports, loads, Pre_member, prearea, step1, step2, step3, fl, jc, material)
        self.mathprob = None
        self.B = None
        self.judgeB = None
        self.f = []
        self.dof = []
        self.con = True
    def initial(self):
        self.groundstructure.createnodes()
        self.groundstructure.createmembers()
        self.groundstructure.addpremember()
        self.groundstructure.createpotentials()
    def processBoundaryCondition(self):
        for nd in self.groundstructure.NDL:
            self.dof.append([nd.support.condition])
        ldcase = max([len(node.load) for node in self.groundstructure.NDL])
        for i in range(ldcase): # input all loadcase into nodes
            self.f.append([nd.load[i].magnitude if i<=len(nd.load)-1 else [0,0,0] for nd in self.groundstructure.NDL])
            self.f[i] = np.array(self.f[i]).flatten()
        self.f, self.dof = np.array(self.f), np.array(self.dof).flatten()
    def MatrixB(self): # this matrix serves as a coefficient matrix of groundstructure members
        NDL, TML = self.groundstructure.NDL, self.groundstructure.TML
        m, n1, n2 = len(TML), np.array([int(mem.start.number) for mem in TML]), np.array([int(mem.end.number) for mem in TML])
        Nd = np.array([NDL[i].position for i in range(len(NDL))])
        l, X, Y, Z = np.array([mem.length for mem in TML]), Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1], Nd[n2,2]-Nd[n1,2]
        d0, d1, d2, d3, d4, d5 = self.dof[n1*3], self.dof[n1*3+1], self.dof[n1*3+2], self.dof[n2*3], self.dof[n2*3+1], self.dof[n2*3+2]
        s = np.concatenate((-X/l * d0, -Y/l * d1, -Z/l*d2, X/l * d3, Y/l * d4, Z/l*d5))
        r = np.concatenate((n1*3, n1*3+1, n1*3+2, n2*3, n2*3+1, n2*3+2))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        self.B = sparse.coo_matrix((s, (r, c)), shape = (len(NDL)*3, m))
    def judgeMatrixB(self): # this matrix serves as a coefficient matrix of potential members 
        NDL, PML = self.groundstructure.NDL, self.groundstructure.PML
        m, n1, n2 = len(PML), np.array([int(mem.start.number) for mem in PML]), np.array([int(mem.end.number) for mem in PML])
        Nd = np.array([NDL[i].position for i in range(len(NDL))])
        l, X, Y, Z = np.array([mem.length for mem in PML]), Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1], Nd[n2,2]-Nd[n1,2]
        d0, d1, d2, d3, d4, d5 = self.dof[n1*3], self.dof[n1*3+1], self.dof[n1*3+2], self.dof[n2*3], self.dof[n2*3+1], self.dof[n2*3+2]
        s = np.concatenate((-X/l * d0, -Y/l * d1, -Z/l*d2, X/l * d3, Y/l * d4, Z/l*d5))
        r = np.concatenate((n1*3, n1*3+1, n1*3+2, n2*3, n2*3+1, n2*3+2))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m),np.arange(m)))
        self.judgeB = (sparse.coo_matrix((s, (r, c)), shape = (len(NDL)*3, m))).tocsc()
    def solve(self): # solve the optimization problem through iterative method
        for itr in range(1,100):
            if self.con == True:
                self.MatrixB()
                self.mathprob = Mathprob(self.B)
                self.mathprob.solveLP(self.groundstructure.NDL, self.groundstructure.TML, self.dof, self.f, self.groundstructure.material.lc, self.groundstructure.material.lt, self.groundstructure.jc)
                for i in self.groundstructure.prenumber:
                    self.mathprob.arealist[i] = self.groundstructure.prearea
                self.mathprob.result = np.dot(self.mathprob.lengthlist , np.transpose(self.mathprob.arealist))
                self.plot("Itr:" + str(itr),True)
                print("Itr: %d,vol: %f,mems:%d" % (itr, self.mathprob.result, len(self.groundstructure.TML)))
                self.judgeMatrixB()
                self.stopviolation()
                self.groundstructure.createpotentials()
            else:
                self.mathprob.result = np.dot(self.mathprob.lengthlist , np.transpose(self.mathprob.arealist))
                break
        print("Final Volume: %f" % (self.mathprob.result))       
    def stopviolation(self): # sort out the violating members
        l, lst = [], []
        for i,mem in enumerate(self.groundstructure.ML):
            if (mem.exist == False):
                l.append(mem.length+self.groundstructure.jc)
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
    def updatedata(self): # update the area of members
        for i in range(len(self.mathprob.arealist)):
            self.groundstructure.TML[i].updatearea(self.mathprob.arealist[i])
    def plot(self, str, update): # plot the structure
        mp = self.mathprob
        plt.ion() if update else plt.ioff()
        plt.clf(); plt.axis('off'); plt.draw()
        ax1 = plt.axes(projection='3d')
        plt.title(str)
        tk = 10/max(mp.arealist)
        for i in [i for i in range(len(mp.arealist)) if mp.arealist[i] >= max(mp.arealist) * 1e-3]:
            if all([mp.qlist[lc][i]>=0 for lc in range(len(mp.qlist))]): c = 'r'
            elif all([mp.qlist[lc][i]<=0 for lc in range(len(mp.qlist))]): c = 'b'
            else: c = 'tab:gray'
            plt.plot([self.groundstructure.TML[i].start.position[0],self.groundstructure.TML[i].end.position[0]],
                     [self.groundstructure.TML[i].start.position[1],self.groundstructure.TML[i].end.position[1]],
                     [self.groundstructure.TML[i].start.position[2],self.groundstructure.TML[i].end.position[2]],c, linewidth = mp.arealist[i] * tk)
        plt.show()    
