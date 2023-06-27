#---------------------------------------
# Imports
#---------------------------------------
import numpy as np
import pandas as pd

#---------------------------------------
# Variables setup
#---------------------------------------
# Geometry on X
Lx = 1
Nx = 51
dx = Lx/(Nx - 1)

# Geometry on Y
Ly = 1
Ny = 51
dy = Ly/(Ny - 1)

Npoints = Nx*Ny

# Mesh
Xv = np.linspace(0, Lx, Nx)
Yv = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(Xv, Yv)

X = np.reshape(X, Npoints)
Y = np.reshape(Y, Npoints)

# Time
Lt = 0.1
Nt = 51
Time = np.linspace(0, Lt, Nt)
dt = Time[1] - Time[0]

# Material is Cordierite 
k = 1    #W/(m*K)
cv = 1   #J/(kg*C)
rho = 1 #kg/m3

alpha = k/(rho*cv)

#---------------------------------------
# Problem Setup
#---------------------------------------
# Definition of the Boundary Conditions points
bc_bottom = np.arange(0, Nx, 1)
bc_left   = np.arange(Nx, Npoints - Nx, Nx)
bc_top    = bc_bottom + Npoints - Nx
bc_right  = bc_left + Nx - 1

bc = np.hstack([bc_bottom, bc_left, bc_top, bc_right])
core = np.hstack([np.arange(Nx + 1, 2*Nx - 1, 1) + Nx*j for j in range(0,Ny - 2)])

Q = 100*rho*cv

A = np.eye(Npoints)
b = np.zeros(Npoints, dtype= 'float64')

for i in core:
    A[i,      i] = 1 + 2*alpha*dt*(1/dx**2 + 1/dy**2)
    A[i,  i - 1] = -alpha*dt/dy**2
    A[i,  i + 1] = -alpha*dt/dy**2
    A[i, i - Nx] = -alpha*dt/dx**2
    A[i, i + Nx] = -alpha*dt/dx**2

    b[i] = Q*dt/(rho*cv)

#---------------------------------------
# Problem solution
#---------------------------------------
T_num_vec = np.zeros((Nt, Npoints))
Ainv = np.linalg.inv(A)

for k in range(1, Nt):
    T_num_vec[k] = Ainv@(T_num_vec[k - 1] + b.T)

#---------------------------------------
# Validation
#---------------------------------------
num_solution = pd.read_csv('transcal_num_solution_indexed.csv',
                            sep=',', 
                            names=["idx_x", "x", "idx_y", "y", "idx_t", "t", "f(x, y, t)"])

num_solution['idx_x'] = num_solution['idx_x'] - 1
num_solution['idx_y'] = num_solution['idx_y'] - 1
num_solution['idx_t'] = num_solution['idx_t'] - 1

num_solution['t'] = (num_solution['t']*1000)//1/1000

num_solution['idx_xy'] = num_solution['idx_x'] + Nx*num_solution['idx_y']
num_solution['calculated f'] = T_num_vec[num_solution['idx_t'], num_solution['idx_xy']]
num_solution['error'] = abs(num_solution['calculated f'] - num_solution['f(x, y, t)'] )

SRME = np.sqrt(sum(num_solution['error']**2)/(Nt*Npoints))
MAE = sum(num_solution['error'])/(Nt*Npoints)
Max  = max(num_solution['error'])

print('O erro quadrático médio é: %.3e' % (SRME))
print('O erro absoluto médio é: %.3e' % (MAE))
print('O erro máximo é: %.3e' % (Max))