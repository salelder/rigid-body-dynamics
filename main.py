from numpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def Rot1(t):
    """Returns a rotation matrix for axis 1 and finite angle t."""
    return array([[1,0,0], [0,cos(t),-sin(t)], [0,sin(t),cos(t)]]);
def Rot2(t):
    """Returns a rotation matrix for axis 2 and finite angle t."""
    return array([[cos(t),0,sin(t)], [0,1,0], [-sin(t),0,cos(t)]]);
def Rot3(t):
    """Returns a rotation matrix for axis 3 and finite angle t."""
    return array([[cos(t),-sin(t),0], [sin(t),cos(t),0], [0,0,1]]);
def Rotd(t1,t2,t3):
    """Returns a rotation matrix for infinitesimal angles t1, t2, t3
    about axes 1, 2, 3, respectively."""
    return array([[1,-t3,t2], [t3,1,-t1], [-t2,t1,1]]);

def RotI(I, t1, t2, t3):
    """."""
        
class RigidBody:
    """Encodes a rigid body and provides methods for applying moments to it.
    At this time, only planar rotations about axis 3 are supported!"""
    def __init__(self, I11, I22, I33, a1, a2, a3, w):
        """a1...3 are the principal axes of the object. They define its orientation.
        w is the object's angular velocity. I11...33 define the inertia tensor
        for the object about the same point that moments will be computed about.
        E.g., In the case of a pendulum, this is most likely the pivot point."""
        self.I11= I11; self.I22= I22; self.I33= I33; self.a1= a1; self.a2= a2; self.a3= a3;
        self.I= array([[I11,0,0], [0,I22,0], [0,0,I33]]);
        self.w= w;
        
    def Rotate(self, R):
        """Rotates the rigid object according to the rotation matrix provided."""
        self.a1= dot(R, self.a1); self.a2= dot(R, self.a2); self.a3= dot(R, self.a3);
        self.I= dot(R, self.I, transpose(R));
        
    def alpha(self, M):
        """Computes the angular acceleration of the object given the total
        external moments M."""
        #In general, M= I alpha + w cross I w, where I is the inertia matrix
        return array([0.,0.,M[2]/self.I33]); # For planar rotations only!

L= 1.
m= 1.
t0= pi/4
g= 9.81
a1= array([sin(t0),-cos(t0),0.]); a2= array([cos(t0),sin(t0),0.]); a3= array([0.,0.,1.])
rod= RigidBody(0,m*L**2/3,m*L**2/3, a1,a2,a3, array([0.,0.,0.]))
def Moment():
    """Return the total external moment on the rod defined above."""
    return array([0.,0.,L*m*g*rod.a1[1]/rod.a2[0]/2.]);
dt=.05; nframes= 500
a1= []; a2= []; # these will store the values of a1 and a2 for time k
k= 0
while k < nframes:
    a1.append(rod.a1[0]); a2.append(rod.a2[0]);
    a= rod.alpha(Moment());
    rod.w+= a*dt;
    rod.a1= dot(Rot3(rod.w[2]), rod.a1);
    k+= 1;
    
fig= plt.figure()
ax= fig.add_subplot(111, xlim=(-1.5, 1.5), ylim=(-1.5,1.5))
line, = ax.plot([],[],'b-')

def animate(k):
    line.set_data(array([0,a1[k][0]]), array([0,a1[k][1]]))  # update the data
    return line,

def init():
    line.set_data(np.array([0,1]), np.array([0,0]))
    return line,
# must assign to something    
a= animation.FuncAnimation(fig, animate, frames=arange(0,nframes), init_func=init, interval=1, blit=True)
plt.show()