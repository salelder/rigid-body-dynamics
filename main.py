from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import ode

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
def solverot(v1, v2, w1, w2):
    """Returns a rotation vector for the transformation of v1 to v2
    or w1 to w2 (in case rotation is along the axis of v1 or w1)."""
    if v1 == None:
        return array([0.,0.,0.])
    cr= cross(v1, v2)
    if linalg.norm(cr) == 0:
        return solverot(w1, w2, None, None)
    else:
        return arcsin(linalg.norm(cr)/(linalg.norm(v1)*linalg.norm(v2))) * unit(cr)

x= array([1., 0., 0.]); y= array([0., 1., 0.]); z= array([0., 0., 1.])

# Read in values from data input file
fname= raw_input("Input file> ");
f= open(fname, 'r')

mA= None; IA= None; mB= None; IB= None; thA= None; thB= None
wA0= None; wB0= None; p0= None; q0= None; r0= None
g= 0.; dt= 0.; T= 0.
if f.readline() == "3D double pendulum\n":
    mA= float(f.readline().translate(None, '\n'))
    IA= reshape([float(k) for k in f.readline().translate(None, '\n').split(' ')], (3,3))
    mB= float(f.readline().translate(None, '\n'))
    IB= reshape([float(k) for k in f.readline().translate(None, '\n').split(' ')], (3,3))
    
    thA= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    thB= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    wA0= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    wB0= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    p0= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    q0= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    r0= array([float(k) for k in f.readline().translate(None, '\n').split(' ')])
    
    g= float(f.readline().translate(None, '\n'))
    dt= float(f.readline().translate(None, '\n'))
    T= float(f.readline().translate(None, '\n'))    
else:
    print "Error: file not recognized!"
f.close()

# Rotate inertia matricies according to input file
CA= RotS123(thA); CB= RotS123(thB)
IA= CA.dot(IA).dot(CA.T); IB= CB.dot(IB).dot(CB.T)

# Generate two additional generalized coordinates,
# needed to fully define orientation in all cases
# vA is a unit vector perpendicular to p, and vB is a unit vector perpendicular to r
vA0= unit(array([0., 1., -p0[1]/p0[2]]))
vB0= unit(array([0., 1., -r0[1]/r0[2]]))

# Compute initial "velocity terms," with 'd' for 'dot'
pd0= cross(wA0, p0); qd0= cross(wA0, q0); vAd0= cross(wA0, vA0)
rd0= qd0 + cross(wB0, r0); vBd0= qd0 + cross(wB0, vB0)

# Solver
def solve(t, state):
    """t is current time.
    state is a list of: p, q, r, vA, vB, pd, qd, rd, vAd, vBd, wA, wB ('d' for 'dot').
    Returns the time-derivative of state."""
    
    p= state[0]; q= state[1]; r= state[2]; vA= state[3]; vB= state[4]
    pd= state[5]; qd= state[6]; rd= state[7]; vAd= state[8]; vBd= state[9]
    wA= state[10]; wB= state[11]    
    
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
    
    alphas= linalg.solve(A, b)
    aA= alphas[:3] # alpha_A
    aB= alphas[3:]
    
    pdd= cross(aA, p) + cross(wA, cross(wA, p))
    qdd= cross(aA, q) + cross(wA, cross(wA, q))
    vAdd= cross(aA, vA) + cross(wA, cross(wA, vA))
    rdd= qdd + cross(aB, r) + cross(wB, cross(wB, r))
    vBdd= qdd + cross(aB, vB) + cross(wB, cross(wB, vB))
    
    return [pd, qd, rd, vAd, vBd, pdd, qdd, rdd, vAdd, vBdd, aA, aB]

# Integration scheme
#r= ode(solve).set_integrator('dopri5')
#r.set_initial_value([p0,q0,r0,vA0,vB0,pd0,qd0,rd0,vAd0,vBd0,wA0,wB0], 0.)

# Euler's method. Easy.
time= 0.
s= [p0,q0,r0,vA0,vB0,pd0,qd0,rd0,vAd0,vBd0,wA0,wB0]
#print solve(0., [p0,q0,r0,vA0,vB0,pd0,qd0,rd0,vAd0,vBd0,wA0,wB0])

out= open('output.txt', 'w')
out.write('[')
while time < T:
    #print s
    
    p1= s[0]; r1= s[2] 
    vA1= s[3]; vB1= s[4] # will use in a moment to update inertia matricies

    sdot= solve(0., s)
    #print sdot
    for k in range(12):
        s[k]+= sdot[k]*dt
    time += dt
    
    # update inertia matricies
    p2= s[0]; r2= s[2] 
    vA2= s[3]; vB2= s[4]
    rotA= solverot(p1,p2,vA1,vA2); rotB= solverot(r1,r2,vB1,vB2)
    CA= RotS123(rotA); CB= RotS123(rotB)
    IA= CA.dot(IA).dot(CA.T); IB= CB.dot(IB).dot(CB.T)

    q= s[1]    
    
    out= open('output.txt', 'a')
    out.write('[['+str(p2[0])+','+str(p2[1])+','+str(p2[2])+'],['+str(q[0])+','+str(q[1])+','+str(q[2])+'],['+str(q[0]+r2[0])+','+str(q[1]+r2[1])+','+str(q[2]+r2[2])+']],')
    out.close()
# Runge-Kutta integration; had some trouble
#while r.t < T:
#    print r.t
#    # r.y is state variable
#    p1= r.y[0]; r1= r.y[2] 
#    vA1= r.y[3]; vB1= r.y[4] # will use in a moment to update inertia matricies
#    
#    r.integrate(r.t + dt) # integrates state variable!
#       
#    # update inertia matricies
#    p2= r.y[0]; r2= r.y[2] 
#    vA2= r.y[3]; vB2= r.y[4]
#    #rotA= solverot(p1,p2,vA1,vA2); rotB= solverot(r1,r2,vB1,vB2)
#    #CA= RotS123(rotA); CB= RotS123(rotB)
#    #IA= CA.dot(IA).dot(CA.T); IB= CB.dot(IB).dot(CB.T)