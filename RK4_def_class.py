#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from Accelerations import *
from Schwarschild_function import schwarzschild_2


# In[2]:


a = 1
b = 0.05
Mbh = 1e-1


# In[3]:


class MultipleStars:
    def __init__(self, number_stars, a, b):
        self.number_stars = number_stars
        self.a = a
        self.b = b
        self.x, self.y, self.z, self.vx, self.vy, self.vz = schwarzschild_2(self.a, self.b, self.number_stars)

    def RK4 (self, h, N):
        
        position3D = np.zeros((3, self.number_stars, int(N//h)+1))
        velocity3D = np.zeros((3, self.number_stars, int(N//h)+1))
        
        for j in range(0, self.number_stars):
            
            i = 0

            x = self.x[j]
            y = self.y[j]
            z = self.z[j]
            vx = self.vx[j]
            vy = self.vy[j]
            vz = self.vz[j]

            for k in np.arange(0,N,h):

                #order1
                kx1 = vx*h
                ky1 = vy*h
                kz1 = vz*h
                ku1 = h*acceleration_x( x, y, z)
                kv1 = h*acceleration_y( x, y, z)
                kw1 = h*acceleration_z( x, y, z)

                #order2
                kx2 = (vx+0.5*ku1)*h
                ky2 = (vy+0.5*kv1)*h
                kz2 = (vz+0.5*kw1)*h
                ku2 = h*acceleration_x( x+kx1/2, y+ky1/2, z+kz1/2)
                kv2 = h*acceleration_y( x+kx1/2, y+ky1/2, z+kz1/2)
                kw2 = h*acceleration_z( x+kx1/2, y+ky1/2, z+kz1/2)

                #order3
                kx3 = (vx+0.5*ku2)*h
                ky3 = (vy+0.5*kv2)*h
                kz3 = (vz+0.5*kw2)*h
                ku3 = h*acceleration_x( x+kx2/2, y+ky2/2, z+kz2/2)
                kv3 = h*acceleration_y( x+kx2/2, y+ky2/2, z+kz2/2)
                kw3 = h*acceleration_z( x+kx2/2, y+ky2/2, z+kz2/2)

                #order4
                kx4 = (vx+ku3)*h
                ky4 = (vy+kv3)*h
                kz4 = (vz+kw3)*h
                ku4 = h*acceleration_x( x+kx3, y+ky3, z+kz3)
                kv4 = h*acceleration_y( x+kx3, y+ky3, z+kz3)
                kw4 = h*acceleration_z( x+kx3, y+ky3, z+kz3)

                #calculation of the positions
                x += (kx1+2*kx2+2*kx3+kx4)/6.
                y += (ky1+2*ky2+2*ky3+ky4)/6.
                z += (kz1+2*kz2+2*kz3+kz4)/6.

                #calculation of the velocities
                vx += (ku1+2*ku2+2*ku3+ku4)/6.
                vy += (kv1+2*kv2+2*kv3+kv4)/6.
                vz += (kw1+2*kw2+2*kw3+kw4)/6.

                #adding to the list to get the informations
                position3D[0,j,i]=x
                position3D[1,j,i]=y
                position3D[2,j,i]=z
                velocity3D[0,j,i]=vx
                velocity3D[1,j,i]=vy
                velocity3D[2,j,i]=vz

                i=i+1           

        return position3D, velocity3D
    
    def RK4_perturbation (self, h, N, X, Y, Z):
        
        position3D = np.zeros((3, self.number_stars, int(N//h)+1))
        velocity3D = np.zeros((3, self.number_stars, int(N//h)+1))
        
        
        
        for j in range(0, self.number_stars):
            
            i = 0

            x = self.x[j]
            y = self.y[j]
            z = self.z[j]
            vx = self.vx[j]
            vy = self.vy[j]
            vz = self.vz[j]
            x_bh = X[0]
            y_bh = Y[0]
            z_bh = Z[0] 


            for k in np.arange(0,N,h):

                #order1
                kx1 = vx*h
                ky1 = vy*h
                kz1 = vz*h
                ku1 = h*acceleration_x_star( x, y, z, X[i], Y[i], Z[i])
                kv1 = h*acceleration_y_star( x, y, z, X[i], Y[i], Z[i])
                kw1 = h*acceleration_z_star( x, y, z, X[i], Y[i], Z[i])

                #order2
                kx2 = (vx+0.5*ku1)*h
                ky2 = (vy+0.5*kv1)*h
                kz2 = (vz+0.5*kw1)*h
                ku2 = h*acceleration_x_star( x+kx1/2, y+ky1/2, z+kz1/2, X[i], Y[i], Z[i])
                kv2 = h*acceleration_y_star( x+kx1/2, y+ky1/2, z+kz1/2, X[i], Y[i], Z[i])
                kw2 = h*acceleration_z_star( x+kx1/2, y+ky1/2, z+kz1/2, X[i], Y[i], Z[i])

                #order3
                kx3 = (vx+0.5*ku2)*h
                ky3 = (vy+0.5*kv2)*h
                kz3 = (vz+0.5*kw2)*h
                ku3 = h*acceleration_x_star( x+kx2/2, y+ky2/2, z+kz2/2, X[i], Y[i], Z[i])
                kv3 = h*acceleration_y_star( x+kx2/2, y+ky2/2, z+kz2/2, X[i], Y[i], Z[i])
                kw3 = h*acceleration_z_star( x+kx2/2, y+ky2/2, z+kz2/2, X[i], Y[i], Z[i])

                #order4
                kx4 = (vx+ku3)*h
                ky4 = (vy+kv3)*h
                kz4 = (vz+kw3)*h
                ku4 = h*acceleration_x_star( x+kx3, y+ky3, z+kz3, X[i], Y[i], Z[i])
                kv4 = h*acceleration_y_star( x+kx3, y+ky3, z+kz3, X[i], Y[i], Z[i])
                kw4 = h*acceleration_z_star( x+kx3, y+ky3, z+kz3, X[i], Y[i], Z[i])

                #calculation of the positions
                x += (kx1+2*kx2+2*kx3+kx4)/6.
                y += (ky1+2*ky2+2*ky3+ky4)/6.
                z += (kz1+2*kz2+2*kz3+kz4)/6.

                #calculation of the velocities
                vx += (ku1+2*ku2+2*ku3+ku4)/6.
                vy += (kv1+2*kv2+2*kv3+kv4)/6.
                vz += (kw1+2*kw2+2*kw3+kw4)/6.

                #adding to the list to get the informations
                position3D[0,j,i]=x
                position3D[1,j,i]=y
                position3D[2,j,i]=z
                velocity3D[0,j,i]=vx
                velocity3D[1,j,i]=vy
                velocity3D[2,j,i]=vz

                i=i+1           

        return position3D, velocity3D
 


# In[4]:


def RK4_BH (x, y, z, vx, vy, vz, h, N):
    
    X = []
    Y = []
    Z = []
    Vx = []
    Vy = []
    Vz = []
    
    for k in np.arange (0,N,h):

        #order1
        kx1 = vx*h
        ky1 = vy*h
        kz1 = vz*h
        ku1 = h*acceleration_x_BH( x, y, z)
        kv1 = h*acceleration_y_BH( x, y, z)
        kw1 = h*acceleration_z_BH( x, y, z)

        #order2
        kx2 = (vx+0.5*ku1)*h
        ky2 = (vy+0.5*kv1)*h
        kz2 = (vz+0.5*kw1)*h
        ku2 = h*acceleration_x_BH( x+kx1/2, y+ky1/2, z+kz1/2)
        kv2 = h*acceleration_y_BH( x+kx1/2, y+ky1/2, z+kz1/2)
        kw2 = h*acceleration_z_BH( x+kx1/2, y+ky1/2, z+kz1/2)

        #order3
        kx3 = (vx+0.5*ku2)*h
        ky3 = (vy+0.5*kv2)*h
        kz3 = (vz+0.5*kw2)*h
        ku3 = h*acceleration_x_BH( x+kx2/2, y+ky2/2, z+kz2/2)
        kv3 = h*acceleration_y_BH( x+kx2/2, y+ky2/2, z+kz2/2)
        kw3 = h*acceleration_z_BH( x+kx2/2, y+ky2/2, z+kz2/2)

        #order4
        kx4 = (vx+ku3)*h
        ky4 = (vy+kv3)*h
        kz4 = (vz+kw3)*h
        ku4 = h*acceleration_x_BH( x+kx3, y+ky3, z+kz3)
        kv4 = h*acceleration_y_BH( x+kx3, y+ky3, z+kz3)
        kw4 = h*acceleration_z_BH( x+kx3, y+ky3, z+kz3)

        #calcul of the positions
        x += (kx1+2*kx2+2*kx3+kx4)/6.
        y += (ky1+2*ky2+2*ky3+ky4)/6.
        z += (kz1+2*kz2+2*kz3+kz4)/6.

        #calcul of the celerities
        vx += (ku1+2*ku2+2*ku3+ku4)/6.
        vy += (kv1+2*kv2+2*kv3+kv4)/6.
        vz += (kw1+2*kw2+2*kw3+kw4)/6.

        #adding to the list to get the informations
        X.append(x)
        Y.append(y)
        Z.append(z)
        Vx.append(vx)
        Vy.append(vy)
        Vz.append(vz)
        
    return X,Y,Z
    

