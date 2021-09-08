from shapely.geometry import Polygon
import numpy as np
class Support:
    def __init__(self, position, condition):
        self.position = position
        self.condition = condition

class Load:
    def __init__(self, position, magnitude, loadcase):
        self.position = position
        self.magnitude = magnitude

class Domain:
    def __init__(self, vertex, material):
        self.vertex = vertex
        self.material = material
        self.width = max(np.array(vertex)[:,0])
        self.height = max(np.array(vertex)[:,1])
        self.poly = Polygon(vertex)
        self.convex = self.verify()
    def verify(self):
        self.convex = True if (self.poly.convex_hull.area == self.poly.area) else False
    
class Material:
    def __init__(self, lc, lt):
        self.lt = lc
        self.lc = lt




 
        
