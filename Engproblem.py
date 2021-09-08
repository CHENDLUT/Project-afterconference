from shapely.geometry import Polygon
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
    def __init__(self, vertex):
        self.vertex = vertex
        self.convex = True
        self.width = max(np.array(vertex)[:,0])
        self.height = max(np.array(vertex)[:,1])
        self.poly = Polygon(vertex)
    def verify(self):
        self.convex = True if self.poly.convex_hull.area == self.poly.area else False
        
class Material:
    def __init__(self,lc, lt):
        self.lc = lc
        self.lt = lt
