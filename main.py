from numpy import *
import matplotlib.pyplot as plt

def Rot1(t):
    """Returns a rotation matrix for axis 1 and finite angle t."""
    return array([[1.,0.,0.], [0.,cos(t),-sin(t)], [0.,sin(t),cos(t)]]);
def Rot2(t):
    """Returns a rotation matrix for axis 2 and finite angle t."""
    return array([[cos(t),0.,sin(t)], [0.,1.,0.], [-sin(t),0.,cos(t)]]);
def Rot3(t):
    """Returns a rotation matrix for axis 3 and finite angle t."""
    return array([[cos(t),-sin(t),0.], [sin(t),cos(t),0.], [0.,0.,1.]]);
def Rotd(t):
    """Returns a rotation matrix for array of infinitesimal angles t1, t2, t3
    about axes 1, 2, 3, respectively."""
    t1= t[0]; t2= t[1]; t3= t[2];    
    return array([[1.,-t3,t2], [t3,1.,-t1], [-t2,t1,1.]]);

def unit(v):
    """Returns a unit vector as an array with direction given by array v."""
    return v/linalg.norm(v)
        
class RigidBody:
    """Encodes a rigid body and provides methods for applying moments to it.
    At this time, only planar rotations about axis 3 are supported!"""
    def __init__(self, I11, I22, I33, a1, a2, a3, w):
        """a1...3 are the principal axes of the object. They define its orientation.
        w is the object's angular velocity. I11...33 define the inertia tensor
        for the object about the same point that moments will be computed about.
        E.g., In the case of a pendulum, this is most likely the pivot point."""
        self.I11= I11; self.I22= I22; self.I33= I33; self.a1= a1; self.a2= a2; self.a3= a3;
        self.I= array([[I11,0.,0.], [0.,I22,0.], [0.,0.,I33]]);
        self.w= w;
        
    def rotate(self, R):
        """Rotates the rigid object according to the rotation matrix provided."""
        self.a1= unit(dot(R, self.a1))
        self.a2= unit(dot(R, self.a2))
        self.a3= unit(dot(R, self.a3)) # rescale to unit vectors to help correct for errors
        
        self.I= dot(dot(R, self.I), R.transpose());
        
    def alpha(self, M):
        """Computes the angular acceleration of the object given the total
        external moments M."""
        #In general, M= I alpha + w cross I w, where I is the inertia matrix
        return array([0.,0.,M[2]/self.I33]); # For planar rotations only!

L= 1.
m= 1.
t0= pi/6
g= 9.81

rod= RigidBody(0.,(m*L**2)/3,1.,
    array([sin(t0),-cos(t0),0.]),array([cos(t0),sin(t0),0.]),array([0.,0.,1.]),
    array([0.,0.,0.]))
def Moment():
    """Return the total external moment on the rod defined above."""
    return array([0.,0.,-.5*L*m*g*rod.a1[0]]);

time= []
ypos= [] # of the tip of the pendulum    
dt=.001
k= 0
while k < 5/dt:
    time.append(k*dt)
    ypos.append(rod.a1[1]) # y-component of principal axis 1    
    
    a= rod.alpha(Moment())
    rod.rotate(Rotd(rod.w*dt))
    rod.w+= a*dt;
    #print(rod.a1[0]) # x-component of the pendulum
    k+= 1

plt.plot(time, ypos)
T= 2.*pi*sqrt(2.*L/(3.*g)) # the mathematically derived period
plt.plot([T,T],[-1,0],'r--')
plt.show()