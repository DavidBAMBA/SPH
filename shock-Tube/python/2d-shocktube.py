import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from matplotlib.animation import FuncAnimation

# Parameters
N = 400
mass_value = 0.01
kappa = 2  # parameter to evaluate the h
nu = 1.4

mass = np.zeros(N) + mass_value

# Rectangle dimensions
x_len = 1.2  # total length in x-direction
y_len = 0.2  # total length in y-direction

# Calculate number of particles in each direction
n_x = int(np.sqrt(N * x_len/y_len))
n_y = int(N/n_x)

# Make sure N = n_x * n_y. Adjust if necessary
while n_x * n_y != N:
    n_x -= 1
    n_y = N // n_x

# Define the uniform grids in x and y directions
x_coords = np.linspace(-x_len/2, x_len/2, n_x, endpoint=False) + (x_len/n_x)/2
y_coords = np.linspace(-y_len/2, y_len/2, n_y, endpoint=False) + (y_len/n_y)/2

# Create a meshgrid to obtain the positions of the particles
X, Y = np.meshgrid(x_coords, y_coords)

# Flatten the arrays to get coordinates of each particle
x = X.flatten()
y = Y.flatten()

# Assign properties based on x-coordinates 
# (you can adjust this logic as per your requirements)
rho = np.where(x <= 0, 1, 0.25)
e = np.where(x <= 0, 2.5, 1.795)
p = np.where(x <= 0, 1, 0.1795)

# Initial velocities are the same, set at 0
v_x = np.zeros(N)
v_y = np.zeros(N)

# Create the initial state vector as a matrix
S = np.zeros((N, 7))
S[:, 0] = x
S[:, 1] = y
S[:, 2] = rho
S[:, 3] = v_x
S[:, 4] = v_y
S[:, 5] = e
S[:, 6] = p

NParams = len(S[1])

# Reshape matrix into a vector
S = S.reshape(NParams * N)

def h_len(mass, density):
    """ Calculates the smoothing length for all the particles for a given
    state vector. """
    return np.zeros(N) + nu*(mass/density)

def smoothingW(dx, h):
    """ Utilizes the relative distances of a pair to calculate the 
    smoothing function. The input relative distance has already been calculated
    with the smoothing factor h."""
    
    ad = (1/(h**2 * np.pi)) # Alpha-d factor for 2D
    R = np.linalg.norm(dx)/h
    
    if R >= 0 and R < 1:
        smoothW = ad*(2/3 - R**2 + 0.5*R**3)
        
    elif R >=1 and R < 2:
        smoothW = ad*(((2-R)**3)/6)
        
    else:
        smoothW = 0.0
    
    return smoothW

def smoothingdW(dx, h):
    """ Utilizes the relative distances of a pair to calculate the 
    derivative of the smoothing function. The input relative distance has 
    already been calculated with the smoothing factor h."""
    
    ad = (1/h**2) # Alpha-d factor for 2D
    R = np.linalg.norm(dx)/h
    q = dx/h

    if R >= 0 and R < 1:
        grad_W = -ad * (2*R - 1.5*q**2) * q / h
        
    elif R >= 1 and R < 2:
        grad_W = ad * (2 - R)**2 * q / (2*h)
        
    else:
        grad_W = np.array([0.0, 0.0])

    return grad_W    

def artificialv(x1, x2, rho1, rho2, v1, v2, energy1, energy2, h1, h2):
    """ Calculates the artificial viscosity for a pair. """
    gamma = 1.4 
    alpha = 1
    beta = 1
    # Relative quantities
    dx = x1 - x2
    dv = v1 - v2
    c1 = np.sqrt((gamma-1)*energy1)
    c2 = np.sqrt((gamma-1)*energy2)
    c_rel = (c1 + c2)/2
    rho_rel = (rho1 + rho2)/2
    h_pair = (h1 + h2)/2    
    theta = 0.1*h_pair
    phi_rel = (h_pair*np.dot(dv, dx))/(np.linalg.norm(dx)**2 + theta**2) 
   
    # Calculate viscosity
    visc = 0
    if np.dot(dv, dx) < 0:
        visc = (-alpha*c_rel*phi_rel + beta*phi_rel**2)/rho_rel
    
    else:
        visc = 0
    
    return visc

def pressure(rho, e):
    """ Calculates the pressure for the given particle. """
    gamma = 1.4
    return (gamma-1)*rho*e

def velocity(mass, density1, density2, pressure1, pressure2, smoothingdW, artvisc):
    """ Calculates the derivative of velocity for a given pair of particles. """
    return -mass*(pressure1/density1**2 + pressure2/density2**2 + artvisc)*smoothingdW


def energy(mass, density1, density2, pressure1, pressure2, vel1, vel2, smoothingdW, artvisc):
    return 0.5*mass*((pressure1/(density1**2)) + (pressure2/(density2**2)) + artvisc)*(vel1-vel2)*smoothingdW

def density(mass, vel1, vel2, smoothingW):
    return mass*(vel1-vel2)*smoothingW

def integrate(t, S):
    S = S.reshape(N, NParams)  # Reshape vector into matrix.

    h = h_len(mass, S[:, 2])  # Updated smoothing lengths for each integrated W

    # Calculate pair interactions
    pair_i = []
    pair_j = []
    dx_s = []
    dy_s = []
    for kk in range(N - 1):
        for j in range(i + 1, N):
            dx = S[i, 0] - S[j, 0]
            dy = S[i, 1] - S[j, 1]
            h_pair = (h[i] + h[j]) / 2
            if np.linalg.norm([dx, dy]) <= kappa * h_pair:
                # Store indexes and dx, dy values
                pair_i.append(i)
                pair_j.append(j)
                dx_s.append(dx)
                dy_s.append(dy)

    npairs = len(pair_i)  # Define number of pairs in the given loop

    # Initialize derivatives array
    dW = np.zeros(np.shape(S))

    # Calculate pairwise interactions
    for k in range(npairs):
        pi = pair_i[k]
        pj = pair_j[k]

        # Update densities
        rho_inc = density(mass[pi], S[pi, 3], S[pj, 3], smoothingW(np.array([dx_s[k], dy_s[k]]), h_pair))
        S[pi, 2] += rho_inc
        S[pj, 2] += rho_inc

        # Calculate pairwise pressures
        S[pi, 6] = pressure(S[pi, 2], S[pi, 5])
        S[pj, 6] = pressure(S[pj, 2], S[pj, 5])

        # Calculate artificial viscosity
        artvisc = artificialv(np.array([S[pi, 0], S[pi, 1]]), np.array([S[pj, 0], S[pj, 1]]),
                              S[pi, 2], S[pj, 2], np.array([S[pi, 3], S[pi, 4]]), np.array([S[pj, 3], S[pj, 4]]),
                              S[pi, 5], S[pj, 5], h[pi], h[pj])

        # Calculate velocity derivatives
        dV = velocity(mass[pj], S[pi, 2], S[pj, 2], S[pi, 6], S[pj, 6], smoothingdW(np.array([dx_s[k], dy_s[k]]), h_pair), artvisc)
        dW[pi, 3:5] += dV
        dW[pj, 3:5] -= dV

        # Calculate derivatives of the internal energy
        dE = energy(mass[pj], S[pi, 2], S[pj, 2], S[pi, 6], S[pj, 6], np.array([S[pi, 3], S[pi, 4]]),
                    np.array([S[pj, 3], S[pj, 4]]), smoothingdW(np.array([dx_s[k], dy_s[k]]), h_pair), artvisc)
        dW[pi, 5] += dE
        dW[pj, 5] -= dE

    # Derivatives of position
    dW[:, 0] = S[:, 3]  # Derivative of the x position is the input x velocity
    dW[:, 1] = S[:, 4]  # Derivative of the y position is the input y velocity

    dW = dW.reshape(N * NParams)
    return dW


#_________________________________________ INTEGRATION ____________________________________________#

# Integration parameters
tstep = 0.005
tmin = 0
tmax = tstep * 40
steps = np.arange(tmin, tmax, tstep)  # Starts from tmin
NSteps = len(steps)

# Integrator setup. Note that you initially set tmin and tend as the same value. 
# I've adjusted tend to be tmax for the entire integration period.
W_int = RK45(integrate, tmin, tmax, S, first_step=tstep, max_step=tstep)
W_i = np.zeros([NSteps, N, NParams])

# Loop for the integration
for i, t in enumerate(steps):  # Looping over steps directly
    W_i[i] = np.array(W_int.y).reshape(N, NParams)  # Select current state, reshape and store

    W_int.step()        
    print(i, t)  # Print step number and time value

#___________________________________________ PLOTTING _____________________________________________#

# Define the state vector for the last integrated timestep.
S_plot = W_i[38]

# Plot densities
plt.figure(figsize=[10,10])
plt.scatter(S_plot[:,0], S_plot[:,1], c=S_plot[:,2], cmap='viridis', marker='o')
plt.colorbar().set_label('$Density \, \, [kg/m^{3}]$')
plt.title('Density plot')
plt.xlabel('$x \, Position \, \,[m]$')
plt.ylabel('$y \, Position \, \,[m]$')
plt.axis('equal')

# Plot velocities (magnitude of velocities is taken)
velocity_magnitude = np.sqrt(S_plot[:,3]**2 + S_plot[:,4]**2)
plt.figure(figsize=[10,10])
plt.scatter(S_plot[:,0], S_plot[:,1], c=velocity_magnitude, cmap='viridis', marker='o')
plt.colorbar().set_label('$Velocity \, [m/s]$')
plt.title('Velocity plot')
plt.xlabel('$x \, Position \, \,[m]$')
plt.ylabel('$y \, Position \, \,[m]$')
plt.axis('equal')

# Plot energies
plt.figure(figsize=[10,10])
plt.scatter(S_plot[:,0], S_plot[:,1], c=S_plot[:,5], cmap='viridis', marker='o')
plt.colorbar().set_label('$Internal energy \, \, [J/kg]$')
plt.title('Energy plot')
plt.xlabel('$x \, Position \, \,[m]$')
plt.ylabel('$y \, Position \, \,[m]$')
plt.axis('equal')

# Plot pressures
plt.figure(figsize=[10,10])
plt.scatter(S_plot[:,0], S_plot[:,1], c=S_plot[:,6], cmap='viridis', marker='o')
plt.colorbar().set_label('$Pressure \, \, [N/m^{2}]$')
plt.title('Pressure plot')
plt.xlabel('$x \, Position \, \,[m]$')
plt.ylabel('$y \, Position \, \,[m]$')
plt.axis('equal')

plt.show()
