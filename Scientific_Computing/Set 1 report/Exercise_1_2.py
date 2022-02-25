# import the packages
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

"""
Numerical solver for time-stepping method, default x range: [0, 1], y range: [0, 1] default t range: [0, 10]
Boundary conditions: c(x, y = 1; t) = 1 and c(x, y = 0; t) = 0; c(x = 0, y; t) = c(x = 1, y; t), initial
conditions: c(x, y; t = 0) = 0 for 0 ≤ x ≤ 1, 0 ≤ y < 1.
"""

delta_x = 1e-1
delta_t = 1e-3
mesh_x = np.linspace(0, 1, int(1 / delta_x))
mesh_t = np.linspace(0, 10, int(1 / delta_t))
n_i = len(mesh_x)
n_k = len(mesh_t)

# create empty matrix to store solution of PDE
c = np.zeros((n_i, n_i))
c[:, -1].fill(1)

c_next = np.zeros((n_i, n_i))
c_next[:, -1].fill(1)


def diffusion_solver(j, D=1, dx=delta_x, dt=delta_t):

    cx = np.zeros(n_i)

    # boundary conditions
    for i in range(n_i):
        if i == 0:
            c_ = c[i, j] + dt*D/dx**2 * (c[i+1, j] + c[-(i+1), j] + c[i, j+1] + c[i, j-1] - 4*c[i, j])
        elif i == n_i-1:
            c_ = c[i, j] + dt*D/dx**2 * (c[0, j] + c[i-1, j] + c[i, j+1] + c[i, j-1] - 4 * c[i, j])
        else:
            c_ = c[i, j] + dt*D/dx**2 * (c[i+1, j] + c[i-1, j] + c[i, j+1] + c[i, j-1] - 4*c[i, j])
        
        cx[i] = c_
    c_next[:, j] = cx

    return cx


def parallel_solving_pde():

    results = np.zeros((n_i, n_i))
    results[:, -1].fill(1)

    for j in range(n_i-2):
        j += 1
        p = mp.Process(target=diffusion_solver, args=(j, ))
        p.start()

    return results


# Create matrix to store solution at every 1000 time points
c_store = np.zeros((int(n_k/100), n_i, n_i))

if __name__ == '__main__':

    for k in range(n_k):

        count = 0

        parallel_solving_pde()
        c = c_next

        if k % 100 == 0:
            c_store[count, :, :] = c
            print(k, ' time step achieved!')
