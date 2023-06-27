#---------------------------------------
# Imports
#---------------------------------------
import numpy as np

#---------------------------------------
# Variables setup
#---------------------------------------
# Geometry on X
Lx = 1
Nx = 50
dx = Lx/(Nx - 1)

# Geometry on Y
Ly = 1
Ny = 50
dy = Ly/(Ny - 1)

Npoints = Nx*Ny

# Mesh
Xv = np.linspace(0, Lx, Nx)
Yv = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(Xv, Yv)

X = np.reshape(X, Npoints)
Y = np.reshape(Y, Npoints)

# Time
dt = 0.001
Nt = 10_000
Time = np.arange(0, Nt*dt, dt)

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

# Definition of manufactured functions
T_manu = lambda i, j:  (X[j]**2 + Y[j]**2)*np.exp(-Time[i])
Q_manu = lambda i, j: -(X[j]**2 + Y[j]**2 + 2*X[j] + 2*Y[j])*np.exp(-Time[i])
T_bc = np.array([[T_manu(i, j) for j in bc] for i in np.arange(0, Nt, 1)])

# Definition of linear system elements
A = np.eye(Npoints)
b = np.zeros(Npoints, dtype= 'float64')

for i in core:
    A[i,      i] = 1 + 2*alpha*dt*(1/dx**2 + 1/dy**2)
    A[i,  i - 1] = -alpha*dt/dy**2
    A[i,  i + 1] = -alpha*dt/dy**2
    A[i, i - Nx] = -alpha*dt/dx**2
    A[i, i + Nx] = -alpha*dt/dx**2

for i, idx in enumerate(bc):
    b[idx] = T_bc[0][i]

for i, idx in enumerate(core):
    b[idx] = Q_manu(0, i)*dt/(rho*cv)

#---------------------------------------
# Problem solution
#---------------------------------------
T_vec = np.zeros((Nt, Npoints))

# Initializing Temperature
for idx in np.hstack([bc, core]):
    T_vec[0, idx] = T_manu(0, idx)

Ainv = np.linalg.inv(A)

for k in range(1, Nt):
    for i, idx in enumerate(bc):
        b[idx] = T_manu(k, idx) - T_manu(k - 1, idx)

    for i, idx in enumerate(core):
        b[idx] = Q_manu(k, idx)*dt/(rho*cv)

    T_vec[k] = Ainv@(T_vec[k - 1] + b.T)

Manu_solution_vec = np.array([[T_manu(i,j)for j in np.arange(0, Npoints, 1)] for i in np.arange(0, Nt, 1)])

Error_vec = abs(T_vec - Manu_solution_vec)

SRME = np.sqrt(sum(sum(Error_vec**2))/(Nt*Npoints))
MAE = sum(sum(Error_vec))/(Nt*Npoints)
Max  = max(np.reshape(Error_vec, (Nt*Npoints)))

print('O erro quadrático médio é: %.3e' % (SRME))
print('O erro absoluto médio é: %.3e' % (MAE))
print('O erro máximo é: %.3e' % (Max))