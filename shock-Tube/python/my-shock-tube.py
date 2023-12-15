import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation


def integration(S):
    # Integration parameters
    tstep  = 0.005
    tmax   = tstep*40
    steps  = np.arange(tstep, tmax, tstep)
    NSteps = len(steps)
    print('Num_Steps: ', NSteps)

    # Integrator setup
    S_int = RK45(Equations, 0.005, S, 40, first_step=0.005, max_step=0.005, rtol=1e-9, atol=1e-9)
    S_i   = np.zeros([NSteps, N, NParams])

    # Loop for the integration
    for ii in range(NSteps):
        S_i[ii] = np.array(S_int.y).reshape(N, NParams) 
        S_int.step()
        print('ii: ',ii )#' val: ', S_i[ii,:,4])   

    return S_i
      

def initialize_positions(N, x_left, x_right):
    """Initialize positions of particles based on wall positions."""
    pos = np.zeros(N)
    x_l = np.linspace(-0.6, 0, x_left, endpoint=False)
    x_r = np.linspace( 0, 0.6, x_right, endpoint=False)
    pos[:x_left]  = x_l
    pos[-x_right:] = x_r
    return pos


def initialize_variable(value_left, value_right, N, x_left, x_right):
    """Utility function to initialize particle properties based on wall positions."""
    variable = np.zeros(N)
    variable[:x_left] = value_left
    variable[-x_right:] = value_right
    return variable


def h_len(mass, density, nu, N):
    """ Calculates the smoothing length for all the particles for a given
    state vector. """
    return np.zeros(N) + nu*(mass/density)


def Kernel(dx, h):
    """
    Utilizes the relative distances of a pair to calculate the 
    smoothing function. The input relative distance has already been calculated
    with the smoothing factor h.
    """
    ad = 1 / h  # Alpha-d factor
    R = np.linalg.norm(dx)/h # Assume dx is already an array of differences, so just get the absolute value
    #print('len_R: ',len(R),'  R:  ', R)
    # Compute all possible values
    condition1_val = ad*(2/3 - R**2 + 0.5*R**3)
    condition2_val = ad*(((2-R)**3)/6)
    # Apply conditions and sum results
    q = np.where((0 <= R) & (R < 1), condition1_val, 
              np.where((1 <= R) & (R < 2), condition2_val, 
              0.0))

    return q


def Der_kernel(dx, h):
    """ Utilizes the relative distances of a pair to calculate the 
    derivative of the smoothing function. The input relative distanc has 
    already been calculated with the smoothing factor h."""
    
    ad = 1 / h
    R = np.abs(dx) / h

    # Compute all possible values
    condition1_val =  ad*(-2 + (3/2)*R)*(dx/h**2)
    condition2_val = -ad*(0.5*((2-R)**2))*(dx/(h*abs(dx)))

    # Apply conditions and sum results
    dq = np.where((0 <= R) & (R < 1), condition1_val, 
              np.where((1 <= R) & (R < 2), condition2_val, 
              0.0))

    return dq


def pressure(rho, e):
    """ Calculates the derivative pressure for the given particle. """
    gamma = 1.4
    return (gamma-1)*rho*e


def Euler(mass, rho_i, rho_j, P_i, P_j, dq, artVis):
    """ 
    Calculates the derivative of density for a given pair of particles.
    Assumes that all inputs can be numpy arrays.
    """
    return -mass * (P_i / np.power(rho_i, 2) + P_j / np.power(rho_j, 2) + artVis) * dq


def Energy(mass, rho_i, rho_j, P_i, P_j, dq, artVis, v_i, v_j):
    """ 
    Calculates the energy density for a given pair of particles. 
    Assumes that all inputs can be numpy arrays.
    """
    return 0.5 * mass * ((P_i / np.power(rho_i, 2) + P_j / np.power(rho_j, 2)) + artVis) * (v_i - v_j) * dq


def artVis(x_i, x_j, rho_i, rho_j, v_i, v_j, E_i, E_j, h_i, h_j):
    """ Calculates the artificial viscosity for a pair in a vectorized manner. """
    
    # Define parameters that go into the artificial viscosity
    gamma  = 1.4
    c_i    = np.sqrt((gamma-1) * E_i)
    c_j    = np.sqrt((gamma-1) * E_j)
    alpha  = 1.0
    beta   = 1.0

    # Relative quantities
    dx      = x_i - x_j
    dv      = v_i - v_j
    c_rel   = (c_i + c_j)/2
    rho_rel = (rho_i + rho_j)/2
    h_pair  = (h_i + h_j)/2    
    theta   = 0.1 * h_pair
    phi_rel = (h_pair*np.dot(dv, dx))/(np.linalg.norm(dx)**2 + theta**2) 
   
    # Calculate viscosity using array operations
    condition_mask = np.dot(dv, dx) < 0
    visc_positive  = 0
    visc_negative  = (-alpha * c_rel * phi_rel + beta * phi_rel**2) / rho_rel
    visc           = np.where(condition_mask, visc_negative, visc_positive)
    
    return visc


def Equations(t, S):

    S = S.reshape(N, NParams) # Reshape vector into matrix.

    h = h_len(mass, S[:,1], nu, N) 
    np.triu_indices(N, k=0)
    i, j = np.triu_indices(N, k=0)  # Get indices for unique combinations

    # Compute the differences and h_pair for all possible pairs
    dx = S[i, 0] - S[j, 0]
    h_pair = (h[i] + h[j]) / 2

    condition = np.linalg.norm(dx)<= kappa*h_pair

    pair_i = i[condition]
    pair_j = j[condition]
    dx_s   = dx[condition]
    
    q   = Kernel(dx_s, h_pair[condition])
    dq  = Der_kernel(dx_s, h_pair[condition])

    # Calculate density using the summation density
    S[:,1] = mass * (2 / (3 * h))  # Density self effect (for every particle)
    dS = np.zeros(np.shape(S))  # Empty array to store the derivatives

    # Update densities 
    S[pair_i, 1] += mass[pair_i] * q
    S[pair_j, 1] += mass[pair_j] * q

    
    S[:,4]  = pressure(S[:,1], S[:,3]) # Updates pressure.
    
    # Compute the artificial viscosity for all pairs
    artvisc = artVis(S[pair_i, 0], S[pair_j, 0], S[pair_i, 1], S[pair_j, 1],
                          S[pair_i, 3], S[pair_j, 3], S[pair_i, 2], S[pair_j, 2], h[pair_i], h[pair_j])

    # Compute (derivatives of) ,momentums for each particle from the pairs
    dS[pair_i, 3] += Euler(mass[pair_j], S[pair_i, 1], S[pair_j, 1], S[pair_i, 4], S[pair_j, 4], dq, artvisc)
    dS[pair_j, 3] -= Euler(mass[pair_i], S[pair_j, 1], S[pair_i, 1], S[pair_j, 4], S[pair_i, 4], dq, artvisc)

    # Derivatives of the internal energy
    dS[pair_i, 2] += Energy(mass[pair_j], S[pair_i, 1], S[pair_j, 1],
                        S[pair_i, 4], S[pair_j, 4], dq, artvisc, S[pair_i, 3], S[pair_j, 3])
    dS[pair_j, 2] -= Energy(mass[pair_i], S[pair_j, 1], S[pair_i, 1],
                       S[pair_j, 4], S[pair_i, 4], dq, artvisc, S[pair_j, 3], S[pair_i, 3])

    # Derivatives of density and pressure are 0     
    dS[:,1] = 0 
    dS[:,4] = 0
    dS[:,0] = S[:,3] # Derivative of the position is the input velocity
    
    dS = dS.reshape(N*NParams)
    
    return dS      


def plot_state(S_i):
    """
    Plot the state properties: density, velocity, energy, and pressure.
    :param S_plot: State matrix, each row is a particle's properties, columns are [position, density, velocity, energy, pressure]
    """
    S_plot = S_i[38] 

    # Plot densities
    plt.figure(figsize=[10,10])
    plt.scatter(S_plot[:,0], S_plot[:,1])
    plt.title('Density plot')
    plt.xlabel('$Position \, \,[m]$')
    plt.xlim([-0.4, 0.4])
    plt.ylabel('$Density \, \, [kg/m^{3}]$')
    
    # Plot energy
    plt.figure(figsize=[10,10])
    plt.scatter(S_plot[:,0], S_plot[:,2])
    plt.title('Energy plot')
    plt.xlabel('$Position \, \,[m]$')
    plt.xlim([-0.4, 0.4])
    plt.ylim([1.8, 2.6]) 
    plt.ylabel('$Energy \, [J/kg]$')
    
    # Plot velocities
    plt.figure(figsize=[10,10])
    plt.scatter(S_plot[:,0], S_plot[:,3])
    plt.title('Velocities plot')
    plt.xlabel('$Position \, \,[m]$')
    plt.xlim([-0.4, 0.4])
    plt.ylim([0, 1])
    plt.ylabel('$velocities \, \, [m/s]$')
    
    # Plot pressures
    plt.figure(figsize=[10,10])
    plt.scatter(S_plot[:,0], S_plot[:,4])
    plt.xlim([-0.4, 0.4])
    plt.ylim([0, 1.2])
    plt.title('Pressure plot')
    plt.xlabel('$Position \, \,[m]$')
    plt.ylabel('$Pressure \, \, [N/m^{2}]$')
    
    plt.show()

    fig, ax = plt.subplots(figsize=[10,10])
    line, = ax.plot([], [], 'o', lw=2)
    ax.set_xlim([-0.4, 0.4])
    ax.set_ylim(0, 2) # Este rango es el rango de las posiciones x de las partículas
    ax.set_title('Density')
    ax.set_xlabel('$positions \, \,[m]$')
    ax.set_ylabel('$density \, \,[m]$')

    def init():
        line.set_data([], [])
        return line,

    def update(i):
        S_plot = S_i[i]
        line.set_data(S_plot[:,0], S_plot[:,1]) # Usar las posiciones x para ambos ejes para que solo se muestren las posiciones x
        return line,
    ani = FuncAnimation(fig, update, frames=39, init_func=init, blit=True)
    plt.show()


#-----------------------------------------------------------------------------
#--------------------------------------------MAIN-----------------------------


# Constants
N       = 400
x_left  = 320
x_right = 80

# Initialize Parameters
mass_value = 0.001875
kappa      = 2.0  
nu         = 1.4 
mass       = np.full(N, mass_value) 

# Initialize Positions and Variables
x   = initialize_positions(N, x_left, x_right)
rho = initialize_variable(1, 0.25, N, x_left, x_right)
e   = initialize_variable(2.5, 1.795, N,  x_left, x_right)
v   = np.zeros(N)  
p   = initialize_variable(1, 0.1795, N, x_left, x_right)
S   = np.column_stack((x, rho, e, v, p))

NParams = S.shape[1]

# Reshape matrix into a vector
S = S.reshape(NParams * N)

# Integration
S_i = integration(S)

# Plotting
plot_state(S_i)
    
