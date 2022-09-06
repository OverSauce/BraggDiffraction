import numpy as np
import matplotlib.pyplot as plt

M = 100

eps0 = 8.854e-12
mu0 = 4*np.pi*1e-7
c = 1/np.sqrt(eps0*mu0)

Ez = np.zeros((M+1))
Hy = np.zeros((M))

start, end = M//2 - 25, M//2 + 25
x = np.linspace(-250e-3, 250e-3, M+1)

dx = x[1] - x[0]
dt = dx / c ## Stability Condition
dT = c * dt
Z = np.sqrt(mu0/eps0)
dtmu = dt/dx/mu0
dteps = dt/dx/eps0

C_grid = 3*eps0*dx
L_grid = mu0*dx/4

timesteps = 400

ABSBC = (c*dt - dx) / (c*dt + dx)

plt.ion()
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

for l in range(timesteps):

    Hy_temp = np.zeros((M))
    Ez_temp = np.zeros((M+1))

    for m in range(0, M):
        Hy_temp[m] = Hy[m] + (Ez[m+1] - Ez[m]) / Z * dT/dx
    Hy = Hy_temp

    for m in range(1, M):
        Ez_temp[m] = Ez[m] + (Hy[m] - Hy[m-1]) * Z * dT/dx
    # Ez_temp[0] = Ez[1] + ABSBC * (Ez_temp[1] - Ez[0])
    # Ez_temp[-1] = Ez[-2] + ABSBC * (Ez_temp[-2] - Ez[-1])
    Ez = Ez_temp

    if l < 20:
        Ez[M//2-1:M//2+1] = np.sin(2*np.pi*l/20) * np.exp(-0.5*((l - 40)/20)**2)
    
    ax1.cla()
    ax2.cla()
    ax1.plot(x, Ez)
    ax2.plot(x[1:], Hy)
    ax1.set_title(f"C = {C_grid*1e15:.2f}fF, L = {L_grid*1e12:.2f}pH, $f_c$ = {1/20/dt*1e-9:.0f}GHz")
    ax1.set_xlabel("x (m)")
    ax1.set_ylim([-2, 2])
    ax2.set_ylim([-0.005, 0.005])
    plt.pause(0.01)

plt.ioff()