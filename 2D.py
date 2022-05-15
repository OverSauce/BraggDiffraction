import numpy as np 
import matplotlib.pyplot as plt
import GIF

# sqruare root of 1/2 = 0.7071067811865476
"""
Bragg's Law:
    d = lambda * cos(theta) / 2
    lambda = 7.071
    theta = 45°
    d = 7.071 * cos(45°) / 2 = 7.071 / 2 = 3.535 -> 4
"""
M, N = 100, 100 # M, N = number of rows, columns // careful: first argument of an array is y axis and second is x axis

eps0 = 8.854e-12
mu0 = 4*np.pi*1e-7

c = 1/np.sqrt(eps0*mu0)

Ex = np.zeros((N+1, M+1))
Ey = np.zeros((N+1, M+1))
Hz = np.zeros((N+1, M+1))

eps = eps0 * np.ones((N+1, M+1))
mu = mu0 * np.ones((N+1, M+1))

x = np.linspace(-1e-6, 1e-6, M+1)
y = np.linspace(-1e-6, 1e-6, N+1)

dx = x[1] - x[0]
dy = y[1] - y[0]
dt = 1 / c / np.sqrt((1/dx)**2 + (1/dy)**2)

T = 20
T_real = T * dt

Capacitance = 3*eps0*dx
Inductance = mu0*dy/4

"""
Hard Source

sourcetop = N
sourcesize = 10
sourceposition = ((sourcetop-sourcesize, sourcetop), (0, 0))



"""

# srcsize = np.abs(y[N] - y[N-10]) # 10um
# wavelength = 20 #um

yct = N//2
xct = M//2
λ = 20
n_d = 1 # diffraction order
d = int(n_d*λ/2)
print(d)

# eps[yct+4*d+1:, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct+4*d-1:yct+4*d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct+3*d-1:yct+3*d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct+2*d-1:yct+2*d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct+d-1:yct+d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct-1:yct+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct-d-1:yct-d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct-2*d-1:yct-2*d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct-3*d-1:yct-3*d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[yct-4*d-1:yct-4*d+1, xct-λ//2:xct+λ//2] = 100 * eps0
# eps[:yct-4*d-1, xct-λ//2:xct+λ//2] = 100 * eps0

for n in range(N):
    eps[n, xct-λ//2:xct+λ//2] = eps0 + (1 + np.cos(2*np.pi*n/d)) * eps0

# dT = c*dt
# Z0 = np.sqrt(mu0/eps0)
# Z = np.sqrt(mu/eps)

timesteps = 800

print(f"Sampling Period = {dt}s\nSampling Frequency = {1/dt}Hz")

ABSBC = (c*dt - dx) / (c*dt + dx)

sigma = 10

plt.ion()
# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
plt.figure()

for l in range(timesteps):

    Ex_temp = np.zeros((N+1, M+1))
    Ey_temp = np.zeros((N+1, M+1))
    Hz_temp = np.zeros((N+1, M+1))

    for n in range(0, N):
        for m in range(0, M):
            Hz_temp[n, m] = Hz[n, m] + (Ey[n+1, m] - Ey[n, m])*dt/dx/mu[n, m] - (Ex[n, m+1] - Ex[n, m])*dt/dy/mu[n, m]
    Hz_temp[-1, :] = Hz[-2, :] + ABSBC * (Hz_temp[-2, :] - Hz[-1, :])
    Hz_temp[:, -1] = Hz[:, -2] + ABSBC * (Hz_temp[:, -2] - Hz[:, -1])
    Hz = Hz_temp.copy()

    for n in range(1, N+1):
        for m in range(1, M+1):
            Ex_temp[n, m] = Ex[n, m] - (Hz[n, m] - Hz[n, m-1])*dt/dy/eps[n, m]
            Ey_temp[n, m] = Ey[n, m] + (Hz[n, m] - Hz[n-1, m])*dt/dx/eps[n, m]
    Ex_temp[0, :] = Ex[1, :] + ABSBC * (Ex_temp[1, :] - Ex[0, :])
    Ex_temp[:, 0] = Ex[:, 1] + ABSBC * (Ex_temp[:, 1] - Ex[:, 0])
    Ey_temp[0, :] = Ey[1, :] + ABSBC * (Ey_temp[1, :] - Ey[0, :])
    Ey_temp[N:N-20, 0] = Ey[N:N-20, 1] + ABSBC * (Ey_temp[N:N-20, 1] - Ey[N:N-20, 0])

    Ex = Ex_temp.copy()
    Ey = Ey_temp.copy()
    
    for n in range(N, N-20, -1):
        Ey[n, 1] = 10 * np.sin(2*np.pi*n/λ + 2*np.pi*l/T) * np.exp(-((n-(N-10))**2 / (2*sigma**2)))

    # plt.plot(y, Ey[:, M])
    # plt.ylim([-1, 1])
    plt.pcolormesh(y*1e6, x*1e6, Ey, cmap='inferno', vmin=-1, vmax=1)
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')
    plt.colorbar()
    plt.title(f"t = {l*dt*1e15:.3f}fs (It {l})\nT = {T_real*1e15:.3f}fs, f = {1/T_real*1e-15:.3f}PHz, λ = {λ*dx*1e9:.0f}nm, n = {n_d}")
    # plt.savefig(f"E{l}.png")
    plt.pause(0.0001)
    plt.clf()

plt.ioff()
GIF.GIF("E", timesteps)
