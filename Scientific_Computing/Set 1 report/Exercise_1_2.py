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


def diffusion_solver(i, j, D=1, dx=delta_x, dt=delta_t):

    # boundary conditions
    if i == 0:
        c_ = c[i, j] + dt*D/dx**2 * (c[i+1, j] + c[-(i+1), j] + c[i, j+1] + c[i, j-1] - 4*c[i, j])
    elif i == n_i-1:
        c_ = c[i, j] + dt*D/dx**2 * (c[0, j] + c[i-1, j] + c[i, j+1] + c[i, j-1] - 4 * c[i, j])
    else:
        c_ = c[i, j] + dt*D/dx**2 * (c[i+1, j] + c[i-1, j] + c[i, j+1] + c[i, j-1] - 4*c[i, j])

    return c_


def parallel_solving_pde(j):

    results = []
    c_new = np.zeros(n_i)

    # Define the number of processors used
    pool = mp.Pool(mp.cpu_count())

    for i in range(n_i):
        result = pool.apply_async(diffusion_solver, args=(i, j+1))
        results.append(result)

    pool.close()
    pool.join()

    i = 0
    for r in results:
        c_new[i] = r.get()
        i += 1

    return c_new


# Create matrix to store solution at every 1000 time points
c_store = np.zeros((int(n_k/100), n_i, n_i))

if __name__ == '__main__':

    c_k_1 = c.copy()

    for k in range(n_k):

        count = 0
        c = c_k_1

        for j in range(n_i-2):
            c_next = parallel_solving_pde(j)
            c_k_1[:, j] = c_next

        if k % 100 == 0:
            c_store[count, :, :] = c_k_1
            print(k, ' time step achieved!')
