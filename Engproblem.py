from scipy import sparse,spatial
import numpy as np

class Polyhedron:
    def __init__(self, vertex):
        self.vertex = vertex
        self.hull = spatial.ConvexHull(np.array(vertex))
        self.convex = True
    def verify(self):
        self.convex = True if self.hull.vertices.size == len(self.vertex) else False
    def intersect(self, point): # verify whether the point is in the polyhedron but it can only solve the convex optimization temporarily
        if self.convex == True: 
            vt1 = np.array(self.vertex)+[point.x,point.y,point.z]
            pt_hull = spatial.ConvexHull(vt1)
            return True if self.hull.area == pt_hull.area else False

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
    def __init__(self, vertex):
        self.polyhedron = Polyhedron(vertex)
        self.length = max(np.array(vertex)[:,0])
        self.width = max(np.array(vertex)[:,1])
        self.height = max(np.array(vertex)[:,2])
        self.measure = self.calmeasure()
    def calmeasure(self):
        return np.sqrt(self.length**2+self.width**2+self.height**2)/1000

class Material:
    def __init__(self,lc, lt):
        self.lc = lc
        self.lt = lt
