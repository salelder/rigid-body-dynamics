from numpy import *
from scipy.integrate import ode
import matplotlib.pyplot as plt

y0= [1., 0.] # x, v
t0= 0.
def f(t, y):
    return [y[1], -y[0]]
    
r= ode(f).set_integrator('dopri5')
r.set_initial_value(y0, t0)

t1= 2*pi
dt= .1

t= []
vel= []
x= []
while r.t < t1:
    r.integrate(r.t + dt)
    t.append(r.t)
    vel.append(r.y[1])
    x.append(r.y[0])


plt.plot(t, x)
plt.plot(t, vel)
plt.show()

print x