# %% Based on CUQIpy version v0.3.0
from heat1D_problem import *
import numpy as np
import cuqi
import matplotlib.pyplot as plt

# Check cuqi version
assert cuqi.__version__ == "0.3.0"
print ("CUQIpy version: ", cuqi.__version__)

# %% Benchmark data generation



# Data for the small noise case
np.random.seed = 0
model = create_forward_PDE_model(obs_everywhere_indices)

# exact solution (x_exact)
x_exact =1/30*(1-np.cos(2*np.pi*(L-grid)/(L)))\
                        +1/30*np.exp(-2*(10*(grid-0.5))**2)+\
                         1/30*np.exp(-2*(10*(grid-0.8))**2)

x_exact = cuqi.array.CUQIarray(x_exact, is_par=False,
                                geometry=model.domain_geometry)

y_exact = model(x_exact)
y = cuqi.distribution.Gaussian(y_exact, sigma(small_noise_level, y_exact)**2*np.eye(model.range_dim))
data = y.sample()

# Save data code (uncomment to overwrite data file)
#np.savez("data_small_noise.npz", data=data,
#         x_exact=x_exact,
#         y_exact=y_exact)

# Data for the large noise case
np.random.seed = 0
model = create_forward_PDE_model(obs_left_half_indices)
y_exact = model(x_exact)
y = cuqi.distribution.Gaussian(y_exact, sigma(large_noise_level, y_exact)**2*np.eye(model.range_dim))
data = y.sample()

# Save data code (uncomment to overwrite data file)
#np.savez("data_large_noise.npz", data=data,
#            x_exact=x_exact,
#            y_exact=y_exact)



