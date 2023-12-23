import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from scipy.integrate import RK23
from scipy.optimize import newton
#_______________________________________ INITIAL CONDITIONS _______________________________________#
# Parameters
N = 400 # Number of particles 
mass_value = 0.001875 # Mass of the particles 
kappa = 2 # Used in the NNPS algorithm
nu = 1.4 # Scaling factor for the smoothing factor

mass = np.zeros(N) + mass_value

# Give the particles evenly spaced positions before and after the wall
x = np.zeros(N)
xs1 = np.linspace(-0.6, 0, 320, endpoint=False) # 320 left of the wall
xs2 = np.linspace(0, 0.6, 80, endpoint=False) # 80 right of the wall
x[:len(xs1)] = xs1
x[-len(xs2):] = xs2

# Initial density for the particles
rho = np.zeros(N)
rho[:len(xs1)] = 1 # Density of the particles in the left
rho[-len(xs2):] = 0.25 # Density of the particles in the right

# Initial energy
e = np.zeros(N)
e[:len(xs1)] = 2.5 # Density of the particles in the left
e[-len(xs2):] = 1.795 # Density of the particles in the right

# Initial velocities are the same, set at 0
v = np.zeros(N)

# Initial pressures-
p = np.zeros(N)
p[:len(xs1)] = 1 # Pressure of the particles in the left
p[-len(xs2):] = 0.1795 #  of the particles in the right

# Create the initial state vector as a matrix
W = np.zeros((N, 5)) # Change to W0.
W[:,0] = x
W[:,1] = rho
W[:,2] = v
W[:,3] = e
W[:,4] = p 

NParams = len(W[1])

# Reshape matrix into a vector
W = W.reshape(NParams*N)

#________________________________________ FUNCTIONS _______________________________________________#
def h_len(mass, density):
    """ Calculates the smoothing length for all the particles for a given
    state vector. """
    return np.zeros(N) + nu*(mass/density)

def smoothingW(dx, h):
    """ Utilizes the relative distances of a pair to calculate the 
    smoothing function. The input relative distanc has already been calculated
    with the smoothing factor h."""
    
    ad = (1/h) # Alpha-d factor
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
    derivative of the smoothing function. The input relative distanc has 
    already been calculated with the smoothing factor h."""
    
    ad = (1/h) # Alpha-d factor
    R = np.linalg.norm(dx)/h
    
    if R >= 0 and R < 1:
        smoothdW = ad*(-2 + (3/2)*R)*(dx/h**2)
        
    elif R >=1 and R < 2:
        smoothdW = -ad*(0.5*((2-R)**2))*(dx/(h*abs(dx)))
        
    else:
        smoothdW = 0.0
    
    return smoothdW    

def artificialv(x1, x2, rho1, rho2, v1, v2, energy1, energy2, h1, h2):
    """ Calculates the artificial viscosity for a pair. """
    
    # Define parameters that go into the artificial viscosity
    gamma = 1.4
    c1 = np.sqrt((gamma-1)*energy1)
    c2 = np.sqrt((gamma-1)*energy2)
    alpha = 1
    beta = 1

    # Relative quantities
    dx = x1 - x2
    dv = v1 - v2
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
    """ Calculates the derivative pressure for the given particle. """
    gamma = 1.4
    return (gamma-1)*rho*e

def velocity(mass, density1, density2, pressure1, pressure2, smoothingdW, artvisc):
    """ Calculates the derivative of density for a given pair of particles. """
    return -mass*(pressure1/density1**2 + pressure2/density2**2 + artvisc)*smoothingdW
    
def energy(mass, density1, density2, pressure1, pressure2, vel1, vel2, smoothingdW, artvisc):
    return 0.5*mass*((pressure1/(density1**2)) + (pressure2/(density2**2)) + artvisc)*(vel1-vel2)*smoothingdW

#def density(mass, vel1, vel2, smoothingdW):
#    return mass*(vel1-vel2)*smoothingdW

#______________________________________ FUNCTION TO INTEGRATE _____________________________________#
    
# Nearest neighbor search loop:
def integrate(t, W):
    W = W.reshape(N, NParams) # Reshape vector into matrix.

    # Initialize NNPS algorithm - create empty list for storing the results
    pair_i = []
    pair_j = []
    smoothW = [] 
    smoothdW = []
    h = h_len(mass, W[:,1]) #np.zeros(len(W)) + 0.0075 # Updated smoothing lengths for each integrated W
    dx_s = []
    
    for i in range(N - 1):    
        for j in range(i+1, N):
            dx = (W[i,0] - W[j,0])
            h_pair = (h[i] + h[j])/2
            if np.linalg.norm(dx) <= kappa*h_pair:
                # Store indexes 
                pair_i.append(i)
                pair_j.append(j)
                
                # Calculate smoothing functions for the pair
                smoothW.append(smoothingW(dx, h_pair))
                smoothdW.append(smoothingdW(dx, h_pair))
                
                dx_s.append(dx)
   
    npairs = len(pair_i) # Define number of pairs in the given loop

    # Calculate density using the summation density
    W[:,1] = mass*(2/(3*h)) # Density self effect (for every particle)
    dW = np.zeros(np.shape(W)) # Empty array to store the derivatives

    for k in range(npairs):         
        pi = pair_i[k]
        pj = pair_j[k]
        
        W[pi,1] += mass[pj]*smoothW[k]
        W[pj,1] += mass[pi]*smoothW[k]
    
    W[:,4]  = pressure(W[:,1], W[:,3]) # Updates pressure. Depends only on particle
    
    # Compute the derivatives         
    for k in range(npairs):        
        #W[:,0] = x
        #W[:,1] = rho
        #W[:,2] = v
        #W[:,3] = e
        #W[:,4] = p 
    
        # Get index for each particle of the pair
        pi = pair_i[k]
        pj = pair_j[k]
        
        # Calculate artificial viscosity for each pair
        #artificialv(dx, rho1, rho2, v1, v2, energy1, energy2, h1, h2):
        artvisc = artificialv(W[pi,0], W[pj,0], W[pi,1], W[pj,1],
                              W[pi,2], W[pj,2], W[pi,3], W[pj,3], h[pi], h[pj])
        
        # Compute (derivatives of) velocities for each particle from the pair
        dW[pi,2] += velocity(mass[pj], W[pi,1], W[pj,1], W[pi,4], W[pj,4], smoothdW[k], artvisc)
        dW[pj,2] -= velocity(mass[pi], W[pj,1], W[pi,1], W[pj,4], W[pi,4], smoothdW[k], artvisc)
        
        # Derivatives of the internal energy
        dW[pi,3] += energy(mass[pj], W[pi,1], W[pj,1],
          W[pi,4], W[pj,4], W[pi,2], W[pj,2], smoothdW[k], artvisc)
        dW[pj,3] -= energy(mass[pi], W[pj,1], W[pi,1],
          W[pj,4], W[pi,4], W[pj,2], W[pi,2], smoothdW[k], artvisc)
      
    # Derivatives of density and pressure are 0     
    dW[:,1] = 0 
    dW[:,4] = 0
    dW[:,0] = W[:,2] # Derivative of the position is the input velocity
    
    dW = dW.reshape(N*NParams)
    
    return dW        

#_________________________________________ INTEGRATION ____________________________________________#

# Integration parameters
tstep = 0.005
tmin = 0
tmax = tstep*40
steps = np.arange(tstep, tmax, tstep)
NSteps = len(steps)

# Integrator setup
W_int = RK45(integrate, 0.005, W, 40) #tmax, first_step = tstep, max_step = tstep, rtol=10e-9, atol=10-9)
W_i = np.zeros([NSteps, N, NParams])

# Loop for the integration
for i in range(NSteps):
    W_i[i] = np.array(W_int.y).reshape(N, NParams) # Select current state, reshape and store
    W_int.step()        
    print(i)
    
W_plot = W_i[38]    

# Function to find the roots of!
def f(P, pL, pR, cL, cR, gg):
    a = (gg-1)*(cR/cL)*(P-1) 
    b = np.sqrt( 2*gg*(2*gg + (gg+1)*(P-1) ) )
    return P - pL/pR*( 1 - a/b )**(2.*gg/(gg-1.))

# Solución analítica del problema del tubo de choque de Sod
def SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, x0, T, gg):
    # rL, uL, pL, rR, uR, pR : Initial conditions of the Reimann problem 
    # xs: position array (e.g. xs = [0,dx,2*dx,...,(Nx-1)*dx])
    # x0: THIS IS AN INDEX! the array index where the interface sits.
    # T: the desired solution time
    # gg: adiabatic constant 1.4=7/5 for a 3D diatomic gas
    dx = xs[1]
    Nx = len(xs)
    v_analytic = np.zeros((3,Nx),dtype='float64')

    # compute speed of sound
    cL = np.sqrt(gg*pL/rL)
    cR = np.sqrt(gg*pR/rR)
    # compute P
    P = newton(f, 0.5, args=(pL, pR, cL, cR, gg), tol=1e-12)

    # compute region positions right to left# region R
    c_shock = uR + cR*np.sqrt( (gg-1+P*(gg+1)) / (2*gg) )
    x_shock = x0 + int(np.floor(c_shock*T/dx))
    v_analytic[0,x_shock-1:] = rR
    v_analytic[1,x_shock-1:] = uR
    v_analytic[2,x_shock-1:] = pR
    
    # Region 2
    alpha = (gg+1)/(gg-1)
    c_contact = uL + 2*cL/(gg-1)*( 1-(P*pR/pL)**((gg-1.)/2/gg) )
    x_contact = x0 + int(np.floor(c_contact*T/dx))
    v_analytic[0,x_contact:x_shock-1] = (1 + alpha*P)/(alpha+P)*rR
    v_analytic[1,x_contact:x_shock-1] = c_contact
    v_analytic[2,x_contact:x_shock-1] = P*pR
    
    # Region 3
    r3 = rL*(P*pR/pL)**(1/gg)
    p3 = P*pR
    c_fanright = c_contact - np.sqrt(gg*p3/r3)
    x_fanright = x0 + int(np.ceil(c_fanright*T/dx))
    v_analytic[0,x_fanright:x_contact] = r3
    v_analytic[1,x_fanright:x_contact] = c_contact
    v_analytic[2,x_fanright:x_contact] = P*pR
    
    # Region 4
    c_fanleft = -cL
    x_fanleft = x0 + int(np.ceil(c_fanleft*T/dx))
    u4 = 2 / (gg+1) * (cL + (xs[x_fanleft:x_fanright]-xs[x0])/T )
    v_analytic[0,x_fanleft:x_fanright] = rL*(1 - (gg-1)/2.*u4/cL)**(2/(gg-1))
    v_analytic[1,x_fanleft:x_fanright] = u4
    v_analytic[2,x_fanleft:x_fanright] = pL*(1 - (gg-1)/2.*u4/cL)**(2*gg/(gg-1))

    # Region L
    v_analytic[0,:x_fanleft] = rL
    v_analytic[1,:x_fanleft] = uL
    v_analytic[2,:x_fanleft] = pL

    return v_analytic

# Parámetros de la solución analítica
gg = 1.4  # gamma = C_v / C_p
rL, uL, pL = 1.0, 0.0, 1.0
rR, uR, pR = 0.2, 0.0, 0.1795
Nx = 400
X = 1.2
dx = X / (Nx - 1)
xs = np.linspace(-0.6, 0.6, Nx)
x0 = Nx // 2
T = 0.2
analytic = SodShockAnalytic(rL, uL, pL, rR, uR, pR, xs, x0, T, gg)

# Gráficas de comparación
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

# Densidad
axs[0, 0].scatter(W_plot[:, 0], W_plot[:, 1], label='Numerical Density')
axs[0, 0].plot(xs, analytic[0], label='Analytical Density', color='red')
axs[0, 0].set_title('Density Comparison')
axs[0, 0].set_xlabel('Position [m]')
axs[0, 0].set_ylabel('Density [kg/m^3]')
axs[0, 0].legend()

# Velocidad
axs[0, 1].scatter(W_plot[:, 0], W_plot[:, 2], label='Numerical Velocity')
axs[0, 1].plot(xs, analytic[1], label='Analytical Velocity', color='red')
axs[0, 1].set_title('Velocity Comparison')
axs[0, 1].set_xlabel('Position [m]')
axs[0, 1].set_ylabel('Velocity [m/s]')
axs[0, 1].legend()

# Energía
axs[1, 0].scatter(W_plot[:, 0], W_plot[:, 3], label='Numerical Energy')
axs[1, 0].set_title('Energy Plot')
axs[1, 0].set_xlabel('Position [m]')
axs[1, 0].set_ylabel('Internal Energy [J/kg]')

# Presión
axs[1, 1].scatter(W_plot[:, 0], W_plot[:, 4], label='Numerical Pressure')
axs[1, 1].plot(xs, analytic[2], label='Analytical Pressure', color='red')
axs[1, 1].set_title('Pressure Comparison')
axs[1, 1].set_xlabel('Position [m]')
axs[1, 1].set_ylabel('Pressure [N/m^2]')
axs[1, 1].legend()

plt.tight_layout()
plt.show()