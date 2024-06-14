#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np

a = 1.0
b = 0.05
Mbh = 1e-1


# In[ ]:


def acceleration_x(x,y,z):
    return -x/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5)

def acceleration_y(x,y,z):
    return -y/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5)

def acceleration_z(x,y,z):
    return -z*(1+a/np.sqrt(z*z+b*b))/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,3/2)


# In[1]:


def acceleration_x_BH(x,y,z):
    return -x/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5)

def acceleration_y_BH(x,y,z):
    return -y/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5)

def acceleration_z_BH(x,y,z):
    return -z*(1+a/np.sqrt(z*z+b*b))/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,3/2)


# In[2]:


def acceleration_x_star( x, y, z, x_bh, y_bh, z_bh):
    return -x/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5) +  Mbh*(x_bh-x)/pow((x_bh-x)**2 + (y_bh-y)**2 + (z_bh-z)**2,1.5) 

def acceleration_y_star( x, y, z, x_bh, y_bh, z_bh):
    return -y/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5) + Mbh*(y_bh-y)/pow((x_bh-x)**2 + (y_bh-y)**2 + (z_bh-z)**2,1.5)

def acceleration_z_star( x, y, z, x_bh, y_bh, z_bh):
    return -z*(1+a/np.sqrt(z*z+b*b))/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,3/2) + Mbh*(z_bh-z)/pow((x_bh-x)**2 + (y_bh-y)**2 + (z_bh-z)**2,1.5)

