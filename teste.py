import numpy as np
import matplotlib.pyplot as plt
import random

def lorenz_system(x, y, z, r, b=10, s=6):
    x_dot = b * (y - x)
    y_dot = r * x - y - x * z
    z_dot = x * y - s * z
    return x_dot, y_dot, z_dot
def quantumPendulum(x,y,px,py, gam):
    f = flutuation(x,y,gam)
    x_dot  = px
    y_dot  = py
    px_dot  = -gam*px - x*(1+y*y) + f*(random.random())
    py_dot  = -gam*py - y*(1+x*x) + f*(random.random())
    return x_dot ,y_dot ,px_dot ,py_dot 

def classicalPendulum(x,y,px,py,k):
    x_dot = px
    y_dot  = py
    px_dot  = - x - x*k*y*y
    py_dot  = - y - y*k*x*x
    return x_dot ,y_dot ,px_dot ,py_dot 

def flutuation(x,y,gam):
    return np.sqrt(
            gam*(0.62832912000*(x**2.0 + y**2.0) +
             0.01205153100* x**2.0 * y**2.0  +
             0.56437351570*(x**4.0 + y**4.0) +
             0.13998728990*(x**6.0 + y**6.0) +
             0.06202680930*(x**2.0*y**4.0 + y**2.0*x**4.0) +
             0.03555913955* x**4.0 * y**4.0  +
             0.01777956970*(x**8.0 + y**8.0) + 1.155436517)/
             (1.0 + 0.3345167463*(x**2.0 + y**2.0) +
              0.09871707060*(x**4.0 + y**4.0))**2.0)
            
def rk4(function, x,y,z,w,k,step):
    k1 = np.zeros(4)
    k2 = np.zeros(4)
    k3 = np.zeros(4)
    k4 = np.zeros(4)
    k1[0] = function(x,y,z,w,k)[0]
    k1[1] = function(x,y,z,w,k)[1]
    k1[2] = function(x,y,z,w,k)[2]
    k1[3] = function(x,y,z,w,k)[3]
    
    k2[0] = function(x+step*k1[0]/2,y+step*k1[1]/2,z+step*k1[2]/2,w+step*k1[3]/2,k)[0]
    k2[1] = function(x+step*k1[0]/2,y+step*k1[1]/2,z+step*k1[2]/2,w+step*k1[3]/2,k)[1]
    k2[2] = function(x+step*k1[0]/2,y+step*k1[1]/2,z+step*k1[2]/2,w+step*k1[3]/2,k)[2]
    k2[3] = function(x+step*k1[0]/2,y+step*k1[1]/2,z+step*k1[2]/2,w+step*k1[3]/2,k)[3]

    k3[0] = function(x+step*k2[0]/2,y+step*k2[1]/2,z+step*k2[2]/2,w+step*k2[3]/2,k)[0]
    k3[1] = function(x+step*k2[0]/2,y+step*k2[1]/2,z+step*k2[2]/2,w+step*k2[3]/2,k)[1]
    k3[2] = function(x+step*k2[0]/2,y+step*k2[1]/2,z+step*k2[2]/2,w+step*k2[3]/2,k)[2]
    k3[3] = function(x+step*k2[0]/2,y+step*k2[1]/2,z+step*k2[2]/2,w+step*k2[3]/2,k)[3]

    k4[0] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],w+step*k3[3],k)[0]
    k4[1] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],w+step*k3[3],k)[1]
    k4[2] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],w+step*k3[3],k)[2]
    k4[3] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],w+step*k3[3],k)[3]
    
    x = x + (1/6) * step * (k1[0]+2*k2[0]+2*k3[0]+k4[0])
    y = y + (1/6) * step * (k1[1]+2*k2[1]+2*k3[1]+k4[1])
    z = z + (1/6) * step * (k1[2]+2*k2[2]+2*k3[2]+k4[2])
    w = w + (1/6) * step * (k1[3]+2*k2[3]+2*k3[3]+k4[3])
    return x,y,z,w

dr = 0.0001  # parameter step size
r = np.arange(0, 2, dr)  # parameter range
dt = 0.01  # time step
t = np.arange(0, 10, dt)  # time range

# initialize solution arrays
xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
pxs = np.empty(len(t) + 1)
pys = np.empty(len(t) + 1)

# initial values x0,y0,z0 for the system
xs[0], ys[0], pxs[0], pys[0] = (6, 1, 0, 1)


# Save the plot points coordinates and plot the with a single call to plt.plot
# instead of plotting them one at a time, as it's much more efficient
r_maxes = []
z_maxes = []
r_mins = []
z_mins = []


for R in r:
    # Print something to show everything is running
    #print(f"{R=:.2f}")
    for i in range(len(t)):
        # approximate numerical solutions to system
       xs[i + 1], ys[i + 1], pxs[i + 1],pys[i + 1] = rk4(classicalPendulum,xs[i], ys[i], pxs[i], pys[i], R, dt)
       
    # calculate and save the peak values of the z solution
    for i in range(1, len(xs) - 1):
        # save the local maxima
        if xs[i - 1] < xs[i] and xs[i] > xs[i + 1]:
            r_maxes.append(R)
            z_maxes.append(xs[i])
        # save the local minima
        elif xs[i - 1] > xs[i] and xs[i] < xs[i + 1]:
            r_mins.append(R)
            z_mins.append(xs[i])

    # "use final values from one run as initial conditions for the next to stay near the attractor"
    xs[0], ys[0], pxs[0], pys[0] = xs[i], ys[i], pxs[i], pys[i]


plt.scatter(r_maxes, z_maxes, color="black", s=0.5, alpha=0.2)
plt.scatter(r_mins, z_mins, color="red", s=0.5, alpha=0.2)
plt.show()