import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

kb = 1
epsilon = 1
sigma = 1


class Cell:
    
    def __init__(self, r_atom, m_atom, size, N, T, X=None, V=None):
        
        self.r_atom = r_atom
        self.m_atom = m_atom
        self.size = size
        self.N = N 
        self.T = T
        
        self.V = V #macierz Nx2 - predkosci atomow
        self.V_prev = None
        self.V_next = None
        self.Vp = V #macierz Nx2 - predkosci atomow bez sily oporu 

        if self.V == None:
            self.V = self.V_create()

        self.X = X #macierz Nx2 - polozenia atomow
        self.X_next = None
        if self.X == None:
            self.X = self.X_create()
    
    def V_create(self):
        
        V = np.random.random((self.N,2))
        V = V - np.sum(V, axis=0)/self.N
        Ek = 0
        for i in range(self.N):
            Ek += sum(V[i]**2)

        Ek /= self.N
        f = sqrt(self.T/Ek)
        self.Vp = self.V = V*f
        
        return self.V

    def X_create(self):
        self.X = np.random.random((self.N, 2)) * self.size
        overlap = False
        k = 0
        while True:
            for i in range(self.N-1):
                for j in range(i+1,self.N):
                    if self.count_r(i,j)[0] < 2*self.r_atom:
                        overlap = True
                        break
            if overlap:
                self.X = np.random.random((self.N, 2)) * self.size
                overlap = False
            else:
                break
            k+=1
            if k==100000: raise Exception('Time limit exceeded')
        return self.X

    def count_r(self, a1, a2):

        r_vector = np.array([0.0,0.0])

        x = self.X[a2,0]-self.X[a1,0]
        y = self.X[a2,1]-self.X[a1,1]
        x1 = min(abs(x), self.size - abs(x))
        y1 = min(abs(y), self.size - abs(y))
        
        r_value = sqrt(x1**2 + y1**2)
        
        
        if x1 == abs(x): 
            r_vector[0] = x
        else: 
            r_vector[0] = x - np.sign(x)*self.size
            r_vector[0] = x
        
        if y1 == abs(y): 
            r_vector[1] = y
        else: 
            r_vector[1] = y - np.sign(y)*self.size
        
        return r_value, -r_vector

    def count_forces(self, r_cut_off):

        F = np.zeros((self.N,self.N, 2))
        f_cut_off = -24 * epsilon/sigma * ((2 * (sigma/r_cut_off)**13) - (sigma/r_cut_off)**7)
        
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                r_value, r_vector = self.count_r(i, j)
                if r_value < r_cut_off:
                    f = (-24 * epsilon/sigma * ((2 * (sigma/r_value)**13) - (sigma/r_value)**7) - f_cut_off) * r_vector/r_value
                else: 
                    f = 0
                
                F[i,j] = f
                F[j,i] = -f
        
        return np.sum(F, axis=0)

    def count_potential(self, r_cut_off):
        N = self.N
        U = np.zeros((N,N))
        u_cut_off = 4 * epsilon * ((sigma/r_cut_off)**12 - (sigma/r_cut_off)**6)
        
        for i in range(N-1):
            for j in range(i+1, N):
                r_value, r_vector = self.count_r(i, j)
                if r_value < r_cut_off:
                    u = 4*epsilon*((sigma/r_value)**12 - (sigma/r_value)**6) - u_cut_off
                else: u = 0
                U[i,j] = u
                U[j,i] = u
        
        return np.sum(U, axis=0)

    def update_X(self, X):
        self.X = X%self.size

class Simulation:

    def __init__(self, Cell, dt, r_cut_off):
        
        self.Cell = Cell
        self.dt = dt
        self.r_cut_off = r_cut_off
    
    def steps(self, n, energy=False):
        
        C = self.Cell
        
        X = np.zeros((C.N, n*2))
        V = np.zeros((C.N, n*2))
        
        Ek = np.empty((n, C.N))
        Ep = np.empty((n, C.N))
        
        F = C.count_forces(self.r_cut_off)
        
        C.V_prev = C.V - F/(2*C.m_atom) * self.dt
        C.update_X(C.X)

        for i in range(n): 

            X[:, :2] = C.X
            V[:, :2] = C.V
            
            F = C.count_forces(self.r_cut_off)
            C.V_next = C.V_prev + F/C.m_atom * self.dt
            C.update_X(C.X + C.V_next*self.dt)
            C.V = (C.V_prev+C.V_next)/2
            C.V_prev = C.V_next
            
            X = np.roll(X, -2, axis=1)
            V = np.roll(V, -2, axis=1)

            if energy:
                Ek[i] = 0.5 * C.m_atom * np.sum(C.V**2, axis=1) 
                Ep[i] = C.count_potential(self.r_cut_off)

        if energy: 
            return X,V,Ek,Ep
        
        return X,V



def save_frames(S, n, path):
    
    k = 0
    X,V = S.steps(n)
    
    for i in range(0, n, 25):
        plt.axis([0, 8, 0, 8])
        plt.scatter(X[:,i*2],X[:,i*2+1], c='blue', s=1000, alpha=0.6)
        plt.savefig(path + '/' + '_{:05d}.png'.format(k), bbox_inches='tight', pad_inches=0., dpi=300)
        plt.close()
        k += 1

def plot_energy(S, n ,path):
    
    X,V,Ek,Ep = S.steps(n, energy=True)
    fig, axs = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    Ek = np.sum(Ek, axis=1)
    Ep = np.sum(Ep, axis=1)

    axs[0].plot(Ek)
    axs[1].plot(Ep)
    axs[2].plot(Ek+Ep)
    axs[0].set_title('Energia Kinetyczna')
    axs[1].set_title('Energia Potencjalna')
    axs[2].set_title('Energia Calkowita')
    
    plt.savefig(path, bbox_inches='tight', pad_inches=0., dpi=300)
    plt.close() 

C = Cell(0.5, 1, 8, 16, 3)
S = Simulation(C, 0.0001, 2.5)
save_frames(S, 5000, 'frames')
plot_energy(S, 1000, 'Wykres_energii.png')
