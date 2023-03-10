import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma


def distance_between_particules(ri, rj):

    M = ri.shape[0]
    N = rj.shape[0]
    dr = []
    for ii in range(3):
        dr.append( ri[:,ii].reshape((M,1)) - rj[:,ii].reshape((N,1)).T)

    return dr

def density(r, pos, m, h):

    N = r.shape[0]
    dr = distance_between_particules(r, pos)
    
    rho = np.sum(m * kernel(dr, h), 1).reshape((N, 1))

    return rho

def pressure(rho, k , n):
    return k * rho **(1+1/n)

def acceleration(pos, vel, m, h, k, n, l, nu):
    N = pos.shape[0]
    rho = density(pos, pos, m, h)
    P = pressure(rho, k, n)
   
    dW = gradKer(distance_between_particules(pos, pos), h)
    a = []
    for ii in range(3):
        a.append( - np.sum( m*( P / rho**2 + P.T/rho.T**2 ) * dW[ii], 1 ).reshape((N,1)) )

    ac = np.hstack((a[0], a[1], a[2]))
    ac -= l * pos
    ac -= nu * vel

    return ac

def kernel(dr, h):
    
    r = np.sqrt(dr[0]**2 + dr[1]**2 + dr[2]**2)
    W =(1.0 / (h * np.sqrt(np.pi) ))**3 * np.exp(-r**2 / h**2)

    return  W

def gradKer(r, h):
        
    r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    n = -2 * np.exp(- r**2 / h**2 ) / (h**5 * (np.pi)**(3/2) )
    dW = [n * r[0], n * r[1], n * r[2]]

    return dW


#Parameter of system
N = 400    # Number of particles
t = 0      # current time of the simulation
tEnd = 12     # time at which simulation ends
dt = 0.04   # timestep
Nt = int(np.ceil(tEnd/dt))
M = 2      # star mass
R = 0.75   # star radius
h = 0.1    # smoothing length
k = 0.1    # equation of state constant
n = 1      # polytropic index
nu = 1      # damping
plotRealTime = True # 

np.random.seed(42)

l = 2*k*(1+n)*np.pi**(-3/(2*n)) * ( (M*gamma(5/2+n))/R**3/gamma(1+n))**(1/n) /R**2
m = M / N
pos = np.random.randn(N,3) #array de nx3
vel = np.zeros(pos.shape)

acc = acceleration(pos, vel, m, h, k, n, l, nu)

# prep figure
fig = plt.figure(figsize=(4,5), dpi=80)
grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
ax1 = plt.subplot(grid[0:2,0])
ax2 = plt.subplot(grid[2,0])
rr = np.zeros((100,3))
rlin = np.linspace(0,1,100)
rr[:,0] =rlin
rho_analytic = l / ( 4 * k ) * ( R**2 - rlin**2)

for kk in range(Nt):
    vel += acc * dt/2
    pos += vel * dt
    acc = acceleration( pos, vel, m, h, k, n, l, nu )	
    vel += acc * dt/2
    t += dt
    rho = density( pos, pos, m, h )

	# plot in real time
    if plotRealTime or (kk == Nt-1):
        plt.sca(ax1)
        plt.cla()
        cval = np.minimum((rho-3)/3,1).flatten()
        plt.scatter(pos[:,0],pos[:,1], c=cval, cmap=plt.cm.autumn, s=10, alpha=0.5)
        ax1.set(xlim=(-1.4, 1.4), ylim=(-1.2, 1.2))
        ax1.set_aspect('equal', 'box')
        ax1.set_xticks([-1,0,1])
        ax1.set_yticks([-1,0,1])
        ax1.set_facecolor('black')
        ax1.set_facecolor((.1,.1,.1))
        
        plt.sca(ax2)
        plt.cla()
        ax2.set(xlim=(0, 1), ylim=(0, 3))
        ax2.set_aspect(0.1)
        plt.plot(rlin, rho_analytic, color='gray', linewidth=2)
        rho_radial = density( rr, pos, m, h )
        plt.plot(rlin, rho_radial, color='blue')
        plt.pause(0.001)
        
plt.sca(ax2)
plt.xlabel('radius')
plt.ylabel('density')

# Save figure
plt.savefig('sph.png',dpi=240)
plt.show()
	