import numpy as np
import Optproblem
import Mathproblem
from shapely.geometry import Point, LineString, Polygon
class Domain:
    def __init__(self, vertex):
        self.vertex = vertex
        self.poly = Polygon(self.vertex)
    def verify(self):
        return self.poly.convex_hull.area == self.poly.area
    def size(self):
        vertex = np.array(self.vertex)
        width = max(vertex[:,0])
        height = max(vertex[:,1])
        return width, height

class Support:
    def __init__(self, position, condition):
        self.position = position
        self.condition = condition
        
class Load:
    def __init__(self, position, magnitude, loadcase):
        self.position = position
        self.magnitude = magnitude
        
class Material:
    def __init__(self, lc, lt):
        self.lc = lc
        self.lt = lt
