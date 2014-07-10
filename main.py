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

x= array([1., 0., 0.]); y= array([0., 1., 0.]); z= array([0., 0., 1.])
g= 9.81

# Read in values from data input file
fname= raw_input("Input file> ");
f= open(fname, 'r')

mA= None; IA= None; mB= None; IB= None; thA= None; thB= None;
wA= None; wB= None; p= None; q= None; r= None
if f.readline() == "3D double pendulum\n":
    mA= float(f.readline().translate(None, '\n'))
    IA= reshape([float(k) for k in f.readline().translate(None, '\n').split(' ')], (3,3))
    mB= float(f.readline().translate(None, '\n'))
    IB= reshape([float(k) for k in f.readline().translate(None, '\n').split(' ')], (3,3))
    
    thA= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    thB= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    wA= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    wB= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
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

    bA= cross(wA, IA.dot(wA))-mA*cross(g*z+cross(wA, cross(wA, p)), p)-mB*cross(g*z+cross(wA, cross(wA, q))+cross(wB, cross(wB, r)), q)
    bB= cross(wB, IB.dot(wB))-mB*cross(cross(wA, cross(wA, q))+cross(wB, cross(wB, r))+g*z, r)
    b= array([bA.dot(x), bA.dot(y), bA.dot(z), bB.dot(x), bB.dot(y), bB.dot(z)])
    
    px, py, pz= p; qx, qy, qz= q; rx, ry, rz= r
    A= zeros((6,6)) # A will be the coefficient matrix for [A] alphas = b
    A[0]= [-mA*(pz**2+py**2)-mB*qy**2-IA[0][0], mA*px*py+mB*qx*qy-IA[0][1], mA*px*pz+mB*qx*qz-IA[0][2], -mB*(rz*qz+ry*qy), mB*rx*qy, mB*rx*qz]
    A[1]= [mA*py*px+mB*qy*qx-IA[1][0], -mA*(px**2+pz**2)-mB*(qx**2+qz**2)-IA[1][1], mA*py*pz+mB*qy*qz-IA[1][2], mB*ry*qx, -mB*(rx*qx+rz*qz), mB*ry*qz]
    A[2]= [mA*pz*px+mB*qz*qx-IA[2][0], mA*pz*py+mB*qz*qy-IA[2][1], -mA*(py**2+px**2)-mB*(qy**2+qx**2)-IA[2][2], mB*rz*qx, mB*rz*qy, -mB*(ry*qy+rx*qx)]
    
    A[3]= [-mB*(qz*rz+qy*ry), mB*qx*ry, mB*qx*rz, -mB*(rz**2+ry**2)-IB[0][0], mB*rx*ry-IB[0][1], mB*rx*rz-IB[0][2]]
    A[4]= [mB*qy*rx, -mB*(qx*rx+qz*rz), mB*qy*rz, mB*ry*rx-IB[1][0], -mB*(rx**2+rz**2)-IB[1][1], mB*ry*rz-IB[1][2]]
    A[5]= [mB*qz*rx, mB*qz*ry, -mB*(qy*ry+qx*rx), mB*rz*rx-IB[2][0], mB*rz*ry-IB[2][1], -mB*(ry**2+rx**2)-IB[2][2]]
    
    return linalg.solve(A, b)

print solve()
#plt.plot(time, energy)
#plt.show()
f.close()