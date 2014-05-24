from numpy import *

class Point:
    """Encodes a point in space, and provides methods for displacing it."""
    
    def __init__(self, x,y,z):
        self.r= array([x,y,z]); # position r
        
    def rotate(self, axis, theta):
        """Rotates about axis by amount theta.
        axis is a numpy array giving direction of axis."""
        
    def translateAxis(self, axis, amount):
        """Translates along the axis given by the numpy array axis."""
        
    def translateXYZ(self, disp):
        """Translates in xyz-coordinates by amounts given in numpy array disp."""
        
class PointMass:
    """Encodes a point mass in space."""
    def __init__(self, point, mass):
        self.p= point; self.m= mass;
        
class RigidBody:
    """Encodes a rigid body, composed of discrete point masses."""