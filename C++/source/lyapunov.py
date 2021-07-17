import numpy as np
import matplotlib.pyplot as plt
import random

from numpy.core.multiarray import dot
def lorenzJacobian(x,y,z, xd,yd,zd,r):
    jac = np.zeros(3)
    jac[0] = 16.0*yd - 16.0*xd
    jac[1] = r*xd -z*xd- yd - x*zd
    jac[2] = y*xd +x*yd - (4.0)*zd
    return jac
def lorenz(x, y, z, xd,yd,zd,r):
    x_dot = 16.0* y - 16.0*x
    y_dot = r* x - y - x * z
    z_dot = x * y - 4.0 * z
    return x_dot, y_dot, z_dot
            
def rk4(function, x,y,z, xd,yd,zd,r,step):
    k1 = np.zeros(3,float)
    k2 = np.zeros(3,float)
    k3 = np.zeros(3,float)
    k4 = np.zeros(3,float)
    k1[0] = function(x,y,z,xd,yd,zd,r)[0]
    k1[1] =  function(x,y,z,xd,yd,zd,r)[1]
    k1[2] =  function(x,y,z,xd,yd,zd,r)[2]
    # k1[3] =  function(x,y,z,xd,yd,zd)[3]
    
    k2[0] = function(x+step*k1[0]/2.0,y+step*k1[1]/2.0,z+step*k1[2]/2.0,xd,yd,zd,r)[0]
    k2[1] = function(x+step*k1[0]/2.0,y+step*k1[1]/2.0,z+step*k1[2]/2.0,xd,yd,zd,r)[1]
    k2[2] = function(x+step*k1[0]/2.0,y+step*k1[1]/2.0,z+step*k1[2]/2.0,xd,yd,zd,r)[2]
    # k2[3] = function(x+step*k1[0]/2.0,y+step*k1[1]/2.0,z+step*k1[2]/2.0,xd,yd,zd)[3]

    k3[0] = function(x+step*k2[0]/2.0,y+step*k2[1]/2.0,z+step*k2[2]/2.0,xd,yd,zd,r )[0]
    k3[1] = function(x+step*k2[0]/2.0,y+step*k2[1]/2.0,z+step*k2[2]/2.0,xd,yd,zd,r )[1]
    k3[2] = function(x+step*k2[0]/2.0,y+step*k2[1]/2.0,z+step*k2[2]/2.0,xd,yd,zd,r )[2]
    # k3[3] = function(x+step*k2[0]/2.0,y+step*k2[1]/2.0,z+step*k2[2]/2.0,xd,yd,z)[3]

    k4[0] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],xd,yd,zd,r)[0]
    k4[1] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],xd,yd,zd,r)[1]
    k4[2] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],xd,yd,zd,r)[2]
    # k4[3] = function(x+step*k3[0],y+step*k3[1],z+step*k3[2],xd,yd,z)[3]
    
    x = x + (1.0/6.0) * step * (k1[0]+2*k2[0]+2*k3[0]+k4[0])
    y = y + (1.0/6.0) * step * (k1[1]+2*k2[1]+2*k3[1]+k4[1])
    z = z + (1.0/6.0) * step * (k1[2]+2*k2[2]+2*k3[2]+k4[2])
    coord = np.zeros(3, float)
    coord[0] = x
    coord[1] = y
    coord[2] = z
    return coord, k1

def dotProduct(vector1, vector2):
    result = 0.0
    for i in range(len(vector1)):
        result += vector1[i]*vector2[i]
    return result

def projection(vector1, vector2):
    return dotProduct(vector1,vector2)/dotProduct(vector2,vector2)
def normalize(vector):
    return vector/np.linalg.norm(vector)

x,y,z = (10.0,1.0,0.0)
pertubation = np.zeros(3)
pertubation[0] = 2.0*random.random() - 1.0
pertubation[1] = 2.0*random.random() - 1.0
pertubation[2] = 2.0*random.random() - 1.0
# print(pertubation)
pertubation = normalize(pertubation) 
# print(pertubation)

xd, yd, zd = 0,0,0
# print(xd,yd,zd)
z1 =  np.zeros(3, float)
ux =  np.zeros(3, float)
step = 0.01
iterations = (int)(10/step )
cord = np.zeros(3, float)
# plotazao = np.zeros((3,1))
# xs  = np.zeros(1)
# ys  = np.zeros(1)
# zs  = np.zeros(1)
# xs[0],ys[0],zs[0]= x,y,z

# t  = np.zeros(iterations+1)

# lyapunov[0] = projection(ux, z1)
lya = 0.0
# R = np.zeros(200)
# for i in range(200):
#     R[i] = i*0.1
# print("(", z1[0],", ", z1[1], ", ", z1[2],")" )
# print(pertubation)
lyapunov = []
# for r in R:
r=45.95
dummy = np.zeros(3)
coord, dummy = rk4(lorenz,x,y,z,xd,yd,zd,r, step)
pertubation, ux = rk4(lorenzJacobian, coord[0],coord[1],coord[2],pertubation[0],pertubation[1],pertubation[2],r, step)
# print(ux)
# ux = lorenzJacobian(coord[0],coord[1],coord[2],pertubation[0],pertubation[1],pertubation[2],r)
pertubation = normalize(pertubation)
i = 1.0
lya += np.dot(ux, pertubation)/np.dot(pertubation,pertubation)
aux = 0.0
pertubation = normalize(pertubation)

# print(aux-lya)
while(i<iterations):

    pertubation = normalize(pertubation)

    coord, dummy = rk4(lorenz,coord[0],coord[1],coord[2],xd,yd,zd,r, step)
    # ux = lorenzJacobian(coord[0],coord[1],coord[2],pertubation[0],pertubation[1],pertubation[2],r)
    pertubation, ux = rk4(lorenzJacobian, coord[0],coord[1],coord[2],pertubation[0],pertubation[1],pertubation[2],r, step)
    aux = lya
    pertubation = normalize(pertubation)

    # z1 = rk4(lorenzJacobian,coord[0],coord[1],coord[2],xd,yd,zd, r, step)
    lya += np.dot(ux, pertubation)
    # xs = np.append(xs,ux[0])
    # ys = np.append(ys, ux[1])
    # zs = np.append(zs, ux[2])

    print("lyapunov exp = ", lya)
    
    i = i+1.0
    # xs[i+1],ys[i+1],zs[i+1] = rk4(lorenz,  xs[i],  ys[i],  zs[i],xd, yd, zd, step)

    # lyapunov[i+1] = lyapunov[i] + projection(ux, z1)
    # t[i+1] =t[i] + step 
    # print( z1[0]," ", z1[1], " ", z1[2] )
# # lyapunov.append(lya)
# print()
print((lya/((iterations))))



# print(lyapunov)
# ax = plt.figure().add_subplot(projection='3d')

# ax.plot(xs, ys ,zs, lw=0.5)
# # print(lya)
# # plt.plot(R,lyapunov)
# plt.show() 
