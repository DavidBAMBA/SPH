import numpy as np

def h_len(mass, S1, nu):
    return np.ones_like(mass)

def Kernel(dx_s, h_pair):
    return np.ones_like(dx_s)

def Der_kernel(dx_s, h_pair):
    return np.ones_like(dx_s)

def pressure(S1, S3):
    return np.ones_like(S1)

def artVis(*args):
    return np.ones_like(args[0])

def Euler(*args):
    return np.ones_like(args[0])

def Energy(*args):
    return np.ones_like(args[0])

def integrate(t, W):
    W = W.reshape(N, NParams) 

    pair_i = []
    pair_j = []
    smoothW = [] 
    smoothdW = []
    h = h_len(mass, W[:,1], nu) 
    dx_s = []
    
    for i in range(N - 1):    
        for j in range(i+1, N):
            dx = (W[i,0] - W[j,0])
            h_pair = (h[i] + h[j])/2
            if np.linalg.norm(dx.reshape(-1,1)) <= kappa*h_pair:
                pair_i.append(i)
                pair_j.append(j)
                
                smoothW.append(Kernel(dx, h_pair))
                smoothdW.append(Der_kernel(dx, h_pair))
                
                dx_s.append(dx)

    npairs = len(pair_i) 
    print('\nkernel: ', smoothdW)
    W[:,1] = mass * (2 / (3 * h))
    dW = np.zeros(np.shape(W))

    for k in range(npairs):         
        pi = pair_i[k]
        pj = pair_j[k]
        
        W[pi,1] += mass[pj] * smoothW[k]
        W[pj,1] += mass[pi] * smoothW[k]

    W[:,4] = pressure(W[:,1], W[:,3])
  
    for k in range(npairs):        
        pi = pair_i[k]
        pj = pair_j[k]

        artvisc = artVis(W[pi,0], W[pj,0], W[pi,1], W[pj,1],
                         W[pi,2], W[pj,2], W[pi,3], W[pj,3], h[pi], h[pj])
        
        dW[pi,3] += Euler(mass[pj], W[pi,1], W[pj,1], W[pi,4], W[pj,4], smoothdW[k], artvisc)
        dW[pj,3] -= Euler(mass[pi], W[pj,1], W[pi,1], W[pj,4], W[pi,4], smoothdW[k], artvisc)
        
        dW[pi,2] += Energy(mass[pj], W[pi,1], W[pj,1],
                           W[pi,4], W[pj,4], smoothdW[k], artvisc, W[pi,3], W[pj,3]) 
      
        dW[pj,2] -= Energy(mass[pi], W[pj,1], W[pi,1],
                            W[pj,4], W[pi,4], smoothdW[k], artvisc, W[pj,3], W[pi,3],) 
    dW[:,1] = 0 
    dW[:,4] = 0
    dW[:,0] = W[:,2] 
    
    dW = dW.reshape(N*NParams)
    
    return dW        


def Equations(t, S, NParams, mass, kappa, N, nu):

    S = S.reshape(N, NParams) # Reshape vector into matrix.

    h = h_len(mass, S[:,1], nu) 
    
    i_indx, j_indx = np.triu_indices(N, k=1)  # Get indices for unique combinations

    # Compute the differences and h_pair for all possible pairs
    dx = S[i_indx, 0] - S[j_indx, 0]
    h_pair = (h[i_indx] + h[j_indx]) / 2

    condition_mask = np.linalg.norm(dx)<= kappa*h_pair

    pair_i = i_indx[condition_mask]
    pair_j = j_indx[condition_mask]
    dx_s   = dx[condition_mask]
    
    q   = Kernel(dx_s, h_pair[condition_mask])
    dq  = Der_kernel(dx_s, h_pair[condition_mask])

    
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



# Generate random data
N = 4
NParams = 5
t = np.random.rand()
W = np.random.rand(N * NParams)
S = W.copy()  # Ensuring both functions get the exact same input

mass = np.random.rand(N)
kappa = np.random.rand()
nu = np.random.rand()

# Call both functions
result_integrate = integrate(t, W)
result_Equations = Equations(t, S, NParams, mass, kappa, N, nu)

# Compare results
#print('v1: ',result_integrate )
#print('\n v2: ', result_Equations)
assert np.allclose(result_integrate, result_Equations, atol=1e-1), "Results are not the same!"

print("Test passed! Both functions return the same result.")

