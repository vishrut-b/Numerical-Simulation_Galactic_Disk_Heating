#!/usr/bin/env python
# coding: utf-8


import numpy as np



def density(R, z, a, b):
    return (b*b/(4*np.pi))*((a*R*R+(a+3*np.sqrt(z*z+b*b))*(a+np.sqrt(z*z+b*b))**2)/(pow(R*R+(a+np.sqrt(z*z+b*b))**2,2.5)*pow(z*z+b*b,1.5)))


def schwarzschild_2(a, b, nb_stars):
    
    #-------------------- Computation of initial position -------------------- 
    
    Rmax = a                           # Typical Radius of a Galaxy
    rho_max = density(Rmax, 0, a, b)   # Mass density for the plane z=0, inside the maximum radius = a.
    rho_0 = density(0, 0, a, b)        # Mass density at the center of the galaxy 

    list_x = []
    list_y = []
    list_z = np.zeros(nb_stars)        # For each star, the initial position is on the plane z=0
    
    count = 0
    while count<nb_stars:
        
        # We pick a position for a star from a uniform distribution in the range -a and a :
        x = np.random.uniform(-a, a) # for x axis
        y = np.random.uniform(-a, a) # for y axis 
        
        # We compute the radii of stars 
        R = np.sqrt(x*x+y*y)   
        
        # We pick the density for this star from a uniform distribution 
        rho = np.random.uniform(rho_0, rho_max)
        
        # We verify if the picked star is inside the distribution
        if rho < density(R, 0, a, b):
            list_x.append(x)
            list_y.append(y)
            count+=1
    
    #-------------------- Computation of initial velocity --------------------
    

    v_R = 0.1 #radial velocity                
    sigma_R = v_R/10
    sigma_z = np.sqrt(b*b*density(a, 0, a, b)) 

    # Generate random values from a Gaussian distribution :
    list_vz = np.random.normal(loc=0, scale=sigma_z, size=nb_stars)    # velocity along z
    list_vR = np.random.normal(loc=0, scale=sigma_R, size=nb_stars)    # velocity ellipsoid
    
    # We compute now the velocity elipsoid along x and y :
    list_theta = [np.arctan(list_y[i]/list_x[i]) for i in range(0, nb_stars)]
    list_vx = [list_vR[i]*np.cos(list_theta[i]) for i in range(0, nb_stars)]
    list_vy = [list_vR[i]*np.sin(list_theta[i]) for i in range(0, nb_stars)]
    
    
    # We compute the velocity along x and y to have a circular orbit and we add the velocity ellipsoid
    
    velocity_x = []
    velocity_y = []
    for i in range(0, nb_stars):
        
        R = np.sqrt(list_x[i]**2+list_y[i]**2)    
        S = np.sqrt(R*R+(a+np.sqrt(list_z[i]**2+b**2))**2)
        omega = pow(S,-1.5)
        v = R*omega

        vx = list_y[i]*omega
        vy = -list_x[i]*omega
        vz = 0

        velocity_x.append(vx+list_vx[i])
        velocity_y.append(vy+list_vy[i])
    
    return list_x, list_y, list_z, velocity_x, velocity_y, list_vz

