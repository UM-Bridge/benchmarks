#%%
import numpy as np
import cuqi


#%%
dim = 20  # Number of KL modes
N = 100   # Number of solution nodes
L = 1.0  # Length of the domain
dx = L/(N+1)   # Space step size
decay = 1.5 # Decay rate for the KL expansion
normalizer = 10 # Normalization constant for the KL expansion
grid = np.linspace(dx, L, N, endpoint=False)
small_noise_level = 0.001
large_noise_level = 0.05
obs_everywhere_indices = range(N)
obs_left_half_indices = range(int(N//2))
sigma = lambda noise_level, y_exact:\
    1.0/np.sqrt(N)* noise_level*np.linalg.norm(y_exact)
T = 0.01 # Final time
prior_cov = 1
prior_mean = np.zeros(dim)


def create_forward_PDE_model(obs_indices):
    # Set up test problem
    cfl = 5/11 # The cfl condition to have a stable solution
    dt_approx = cfl*dx**2 # Defining approximate time step size
    num_time_steps = int(T/dt_approx)+1 # Number of time steps

    # Observation grid for the heat problem
    grid_obs = grid[obs_indices]
    
    # Time steps
    time_steps = np.linspace(0, T, num_time_steps, endpoint=True)
    
    # FD diffusion operator
    Dxx = (np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1), -1)
           + np.diag(np.ones(N-1), 1))/dx**2  
    
    # PDE form (diff_op, IC, time_steps)
    def PDE_form(initial_condition, t): return (Dxx, np.zeros(N),
                                                initial_condition)
    
    PDE = cuqi.pde.TimeDependentLinearPDE(
        PDE_form, time_steps, grid_sol=grid, grid_obs=grid_obs)
    
    # Set up geometries for model
    domain_geometry = cuqi.geometry.KLExpansion(grid, decay_rate=decay,
                                                normalizer=normalizer, 
                                                num_modes=dim)

    range_geometry = cuqi.geometry.Continuous1D(grid_obs)
    
    # Prepare model
    model = cuqi.model.PDEModel(PDE, range_geometry, domain_geometry)
    return model

