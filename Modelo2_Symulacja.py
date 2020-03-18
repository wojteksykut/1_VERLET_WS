#!/usr/bin/env python
# coding: utf-8

# In[4]:


from matplotlib import pyplot as plt
import numpy as np


# In[5]:


class atom:
    def __init__(self, x, v, m):
        self.x = x
        self.v = v
        self.m = m
        self.x_prev = None
        self.v_prev = None
    
    def F(self):
        G = 0.01
        M = 100
        r2 = self.x[0]**2 + self.x[1]**2
        wF = G * (M * self.m / r2)
        F = [-(self.x[0] * wF)/ r2**(1/2), -(self.x[1] * wF)/ r2**(1/2)]
        
        return F
                      
    def E(self):
        E_k = self.m * (self.v[0]**2 + self.v[1]**2) / 2
        E_p = (self.F()[0]**2 + self.F()[1]**2)**(1/2)  * (self.x[0]**2 + self.x[1]**2) ** (1/2)
        E_c = E_p + E_k
        
        return (E_k, E_p, E_c)
        
    def Euler_next(self, dt):
        
        F = self.F()
        
        x0 = self.x[0] + self.v[0] * dt + (F[0] * dt ** 2) / (2 * self.m)
        x1 = self.x[1] + self.v[1] * dt + (F[1] * dt ** 2) / (2 * self.m)
        
        v0 = self.v[0] + (F[0] * dt) / self.m
        v1 = self.v[1] + (F[1] * dt) / self.m
        
        self.x_prev = self.x
        self.v_prev = self.v
        
        self.x = [x0,x1]
        self.v = [v0,v1]
        
        return (self.x, self.v)
    
    def Verlet_next(self, dt):
        
        F = self.F()
        
        if self.x_prev == None:
            #self.x_prev = [2 * self.v[0] * dt - self.x[0], 2 * self.v[1] * dt - self.x[1]]
            #print('x_prev = ', self.x_prev )
            self.Euler_next(dt)
        
        x0 = 2 * self.x[0] - self.x_prev[0] + (F[0] * dt ** 2) / self.m
        x1 = 2 * self.x[1] - self.x_prev[1] + (F[1] * dt ** 2) / self.m
        
        v0 = (x0 - self.x_prev[0]) / 2 * dt
        v1 = (x1 - self.x_prev[1]) / 2 * dt
        
        self.x_prev = self.x
        self.v_prev = self.v
        
        self.x = [x0,x1]
        self.v = [v0,v1]
        
        return (self.x, self.v)
    
    
    def LF_next(self, dt):
        
        F = self.F()
        
        if self.x_prev == None:
            self.Euler_next(dt)
        
        v_less = [(self.x[0] - self.x_prev[0])/ dt , (self.x[1] - self.x_prev[1])/ dt]
        v_more = [v_less[0] + F[0] * dt / self.m, v_less[1] + F[1] * dt / self.m]
        
        x = [self.x[0] + v_more[0] * dt, self.x[1] + v_more[1] * dt]
        v = [(v_less[0] + v_more[0])/2 , (v_less[1] + v_more[1])/2]
        
        self.x_prev = self.x
        self.v_prev = self.v
        
        self.x = x
        self.v = v
        
        return (self.x, self.v)


# In[14]:


class simulation:
    def __init__(self, steps, dt):
        self.steps = steps
        self.dt = dt
        
    def do_simulation(self, atom, type = 'V'):
        
        if type == 'E':
            for i in range(self.steps):
                next_x = atom.Euler_next(self.dt)
                if (i%100 == 0):
                    plt.scatter(next_x[0][0], next_x[0][1])
        
        elif type == 'LF':
            for i in range(self.steps):
                next_x = atom.LF_next(self.dt)
                if (i%100 == 0):
                    plt.scatter(next_x[0][0], next_x[0][1])
        
        elif type == 'V':
            for i in range(self.steps):
                next_x = atom.Verlet_next(self.dt)
                if (i%50 == 0):
                    plt.scatter(next_x[0][0], next_x[0][1],marker=r'$\clubsuit$')
                    path = 'C:/Users/wojte/Desktop/obrazki/' + str(i) +'abc.png'
                    plt.savefig(path)
        return True


# In[15]:


A = atom([0,1], [1,0], 0.1)
S = simulation(10000, 0.001)
S.do_simulation(A, 'V')


# In[150]:





# In[ ]:




