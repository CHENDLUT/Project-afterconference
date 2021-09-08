from scipy import sparse,spatial
import numpy as np

class Load:
    def __init__(self, position, magnitude, loadcase):
        self.position = position
        self.magnitude = magnitude
        self.loadcase = loadcase

class Support:
    def __init__(self, position, condition):
        self.position = position
        self.condition = condition
        
class Domain:
    def __init__(self, vertex, step1):
        self.vertex = vertex
        self.convex = True
        self.length = max(np.array(vertex)[:,0])
        self.width = max(np.array(vertex)[:,1])
        self.height = max(np.array(vertex)[:,2])
        self.polylist = [[(i, 0, 0),(i, self.width, 0),(i,self.width,self.height),(i, 0, self.height)] for i in np.arange(0,self.length+step1,step1)]
    def verify(self,area):
        hull = spatial.ConvexHull(np.array(self.vertex))
        self.convex = True if hull.area == area else False
        
class Material:
    def __init__(self,lc, lt):
        self.lc = lc
        self.lt = lt
