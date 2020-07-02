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
    def __init__(self, atom, steps, dt):
        self.steps = steps
        self.dt = dt
        self.atom = atom
        
    def do_simulation(self, type = 'V'):
        X = np.empty((self.steps,2))
        Energy = np.empty((self.steps,3))
        for step in range(self.steps):

            if type == 'E':
                for i in range(self.steps):
                    next_x = self.atom.Euler_next(self.dt)
                    E = self.atom.E()
                    if (i%100 == 0):
                        X[step][0] = next_x[0][0]
                        X[step][1] = next_x[0][1]
                        Energy[step][0] = E[0]
                        Energy[step][1] = E[1]
                        Energy[step][2] = E[2]
            
            elif type == 'LF':
                for i in range(self.steps):
                    next_x = self.atom.LF_next(self.dt)
                    E = self.atom.E()
                    if (i%100 == 0):
                        X[step][0] = next_x[0][0]
                        X[step][1] = next_x[0][1]
                        Energy[step][0] = E[0]
                        Energy[step][1] = E[1]
                        Energy[step][2] = E[2]
            
            elif type == 'V':
                for i in range(self.steps):
                    next_x = self.atom.Verlet_next(self.dt)
                    E = self.atom.E()
                    if (i%100 == 0):
                        X[step][0] = next_x[0][0]
                        X[step][1] = next_x[0][1]
                        Energy[step][0] = E[0]
                        Energy[step][1] = E[1]
                        Energy[step][2] = E[2]
        
        return X, Energy


def save_frames(simulation, steps, prefix, skip=10):
    X = simulation.do_simulation()[0]
    k = 0
    for i in range(steps):
        plt.axis([-2, 2, -2, 2])
        plt.scatter(x=0,y=0,c='black')
        plt.scatter(x=X[i,0],y=X[i,1],c='blue')
        plt.savefig('./frames/'+prefix+'_{:05d}.png'.format(k), bbox_inches='tight', pad_inches=0., dpi=300)
        plt.close() 
        k += 1

def save_energy(simulation, steps, path):
    E = simulation.do_simulation()[1]

    fig, axs = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    axs[0].plot(E[:][0])
    axs[0].set_title('Energia kinetyczna')
    axs[1].plot(E[:][1])
    axs[1].set_title('Energia potencjalna')
    axs[2].plot(E[:][2])
    axs[2].set_title('Energia calkowita')
    
    plt.show()
    plt.savefig(path, bbox_inches='tight', pad_inches=0., dpi=300)
    plt.close() 

A = atom([0,1], [1,0], 0.1)
S = simulation(A, 100, 0.01)
S.do_simulation('V')

save_frames(S,100,'V')
save_energy(S,10000,'energiaV.png')

S.do_simulation('E')

save_frames(S,100,'E')
save_energy(S,10000,'energiaE.png')

S.do_simulation('LF')

save_frames(S,100,'LF')
save_energy(S,10000,'energiaLF.png')
