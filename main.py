from numpy import *
import matplotlib.pyplot as plt

# Mathematical support functions
def Rot1(t):
    """Returns a rotation matrix for axis 1 and finite angle t."""
    return array([[1.,0.,0.], [0.,cos(t),-sin(t)], [0.,sin(t),cos(t)]]);
def Rot2(t):
    """Returns a rotation matrix for axis 2 and finite angle t."""
    return array([[cos(t),0.,sin(t)], [0.,1.,0.], [-sin(t),0.,cos(t)]]);
def Rot3(t):
    """Returns a rotation matrix for axis 3 and finite angle t."""
    return array([[cos(t),-sin(t),0.], [sin(t),cos(t),0.], [0.,0.,1.]]);
def RotS123(t):
    """Returns a rotation matrix for a space 1-2-3 rotation with finite
    angles t1, t2, t3"""
    return Rot3(t[2]).dot(Rot2(t[1])).dot(Rot1(t[0]));
def Rotd(t):
    """Returns a rotation matrix for array of infinitesimal angles t1, t2, t3
    about axes 1, 2, 3, respectively."""
    return array([[1.,-t[2],t[1]], [t[2],1.,-t[0]], [-t[1],t[0],1.]]);

def unit(v):
    """Returns a unit vector as an array with direction given by array v."""
    return v/linalg.norm(v)

# Read in values from data input file
fname= raw_input("Input file> ");
f= open(fname, 'r')

mA= None; IA= None; mB= None; IB= None; thA= None; thB= None; p= None; q= None; r= None
if f.readline() == "3D double pendulum\n":
    mA= float(f.readline().translate(None, '\n'))
    IA= reshape([float(k) for k in f.readline().translate(None, '\n').split(' ')], (3,3))
    mB= float(f.readline().translate(None, '\n'))
    IB= reshape([float(k) for k in f.readline().translate(None, '\n').split(' ')], (3,3))
    
    thA= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    thB= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    p= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    q= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    r= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
else:
    print "Error: file not recognized!"

# Rotate inertia matricies according to input file
CA= RotS123(thA); CB= RotS123(thB);
IA= CA.dot(IA).dot(CA.T); IB= CB.dot(IB).dot(CB.T);

# Generate two additional generalized coordinates
# vA is a unit vector perpendicular to p, and vB is a unit vector perpendicular to r
vA= unit(array([0., 1., -p[1]/p[2]]))
vB= unit(array([0., 1., -r[1]/r[2]]))

# Solver
def solve():
    """Returns an array of alpha's, Ax, Ay, Az, Bx, By, Bz."""
    
    return True

#plt.plot(time, energy)
#plt.show()