import numpy as np
import matplotlib.pyplot as plt

## Dimensions
M, N = 100, 100
timesteps = 400

## Constants
eps0 = 8.854e-12
mu0 = 4*np.pi*1e-7
c = 1/np.sqrt(eps0*mu0)

## Initializing Fields and permittivity
Ex = np.zeros((N+1, M+1))
Ey = np.zeros((N+1, M+1))
Hz = np.zeros((N+1, M+1))
Ex_temp = np.empty((N+1, M+1))
Ey_temp = np.empty((N+1, M+1))
Hz_temp = np.empty((N+1, M+1))
eps = eps0 * np.ones((N+1, M+1))

## Space
x = np.linspace(-5e-6, 5e-6, M+1) # 5e-6 => 10 x 10 um box size
y = np.linspace(-5e-6, 5e-6, N+1)

## Centers
yct = N//2
xct = M//2

## Step size for time and space
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = 1 / c / np.sqrt((1/dx)**2 + (1/dy)**2) ## Von Neumann stability condition

## Source parameters
f_real = 462e12
T_real = 1/f_real
T = T_real/dt 
λ_real = c/f_real
λ = λ_real/dx

src_size = 10
angle = 45 # Incident angle, degrees
λ_src = λ/np.cos(angle*np.pi/180)

m0 = src_size # Endpoints of source in x axis
slope = np.tan(angle*np.pi/180) # Slope of source
n0 = N # Endpoints of source in y axis

n_B = 1 # Diffraction order
Λ = n_B*λ/2/np.sin(angle*np.pi/180) # Diffraction wavelength
Λ_real = Λ*dx # Diffraction wavelength in real space
L = λ_src # Length of Bragg planes

print(f"f = {f_real/1e9:.2f}GHz\nλ = {λ_real/1e-9:.2f}nm\nΛ = {Λ_real/1e-9:.2f}nm")

for n in range(N):
    eps[n, int(xct-L/2):int(xct+L/2)] = eps0 + eps0 * (1 + np.cos(2*np.pi*n/Λ))

plt.figure()
plt.pcolormesh(x*1e6, y*1e6, eps/eps0, cmap='seismic')
plt.title(f"$\epsilon_0$ = {eps0}\nΛ = {Λ_real/1e-9:.2f}nm")
plt.xlabel("x [$\mu$m]")
plt.ylabel("y [$\mu$m]")
plt.colorbar()
plt.show()

ABSBC = (c*dt - dx) / (c*dt + dx)

count = 0

plt.ion()
plt.figure()

for l in range(1, timesteps):

    for n in range(0, N):
        for m in range(0, M):
            Hz_temp[n, m] = Hz[n, m] - (Ey[n, m+1] - Ey[n, m])*dt/dx/mu0 + (Ex[n+1, m] - Ex[n, m])*dt/dy/mu0 
    Hz_temp[-1, :] = Hz[-2, :] + ABSBC * (Hz_temp[-2, :] - Hz[-1, :])
    Hz_temp[:, 0] = Hz[:, 1] + ABSBC * (Hz_temp[:, 1] - Hz[:, 0])
    Hz_temp[:, -1] = Hz[:, -2] + ABSBC * (Hz_temp[:, -2] - Hz[:, -1])
    Hz_temp[0, :] = Hz[1, :] + ABSBC * (Hz_temp[1, :] - Hz[0, :])

    Hz = Hz_temp.copy()

    for n in range(1, N+1):        
        for m in range(1, M+1):
            Ex_temp[n, m] = Ex[n, m] + (Hz[n, m] - Hz[n-1, m])*dt/dy/eps[n, m]
            Ey_temp[n, m] = Ey[n, m] - (Hz[n, m] - Hz[n, m-1])*dt/dx/eps[n, m]

    Ex = Ex_temp.copy()
    Ey = Ey_temp.copy()

    for m in range(src_size):
        Ex[int(n0+(m-m0)*slope), m] = 2 * np.cos(2*np.pi*l/T) * (1/(1+np.exp(-0.22*(l-20))))

    plt.pcolormesh(y*1e6, x*1e6, Ey, cmap='seismic', vmin=-1, vmax=1)
    plt.xlabel('x [μm]')
    plt.ylabel('y [μm]')
    plt.colorbar()
    plt.title(f"It {l}, n = {n_B}\nλ = {λ_real*1e9:.0f}nm, Λ = {Λ_real*1e9:.0f}nm")
    plt.pause(0.0001)
    # if l % 10 == 0:
    #     plt.savefig(f"Final_{count}.png")
    #     count += 1
    plt.clf()

plt.ioff()
